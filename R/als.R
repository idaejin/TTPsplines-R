#' Gaussian TT-ALS with fixed or cGCV λ (R backend).
#' @keywords internal
tt_als_fit <- function(y, basis, ranks, lambda_spec, control, penalty_order = 2,
                       init_cores = NULL, offset = NULL) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  method <- lambda_spec$method
  lambda <- lambda_spec$values %||% lambda_spec$lambda0
  offset <- normalize_offset(offset, length(y))
  intercept <- mean(y - offset)
  yc <- y - offset - intercept
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p, ranks[k + 1L], penalty_order)
  })

  n_eval <- 0L
  n_sweeps <- 0L
  history <- list()
  prev_lam <- lambda
  prev_eta <- NULL
  t0 <- proc.time()[["elapsed"]]
  use_spec <- isTRUE(control$use_spectral_gcv)

  for (sw in seq_len(control$max_sweeps)) {
    for (k in seq_len(d)) {
      L <- left_interfaces(cores, basis)
      R <- right_interfaces(cores, basis)
      Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
      ws <- make_core_workspace(
        yc, Xk, penalties[[k]], lambda[k],
        control$lambda_bounds, control$tol_lambda,
        use_spectral = use_spec
      )
      upd <- update_lambda(method, ws)
      cores[[k]] <- array(upd$g, c(ranks[k], p, ranks[k + 1L]))
      lambda[k] <- upd$lambda
      n_eval <- n_eval + upd$n_eval
    }
    n_sweeps <- sw
    eta <- tt_eta(offset, intercept, cores, basis)
    rss <- sum((y - eta)^2)
    # Full Gaussian penalized objective (same as LBFGS)
    pen_val <- 0
    for (kk in seq_len(d)) {
      g <- as.numeric(cores[[kk]])
      pen_val <- pen_val + 0.5 * lambda[kk] * sum(g * as.numeric(penalties[[kk]] %*% g))
    }
    obj <- 0.5 * rss + pen_val
    d_eta <- if (sw == 1L) NA_real_ else sqrt(mean((eta - prev_eta)^2))
    history[[sw]] <- list(
      sweep = sw, rss = rss, objective = obj, penalty = pen_val,
      lambda = lambda, d_eta = d_eta
    )
    prev_eta <- eta
    if (control$trace) {
      cat(sprintf("  ALS sweep %2d | obj=%.6g | RSS=%.6g | λ=%s\n",
                  sw, obj, rss, paste(sprintf("%.3g", lambda), collapse = ",")))
    }
    if (identical(method, "cGCV")) {
      dlog <- max(abs(log(lambda) - log(pmax(prev_lam, 1e-12))))
      if (dlog < control$tol_lambda && sw > 1L) break
      prev_lam <- lambda
    } else if (sw > 2L) {
      prev_rss <- history[[sw - 1L]]$rss
      if (abs(prev_rss - rss) / max(1, abs(prev_rss)) < control$tol) break
    }
  }

  eta <- tt_eta(offset, intercept, cores, basis)
  list(
    cores = cores,
    intercept = intercept,
    lambda = lambda,
    ranks = ranks,
    eta = eta,
    mu = eta,
    deviance = sum((y - eta)^2),
    n_sweeps = n_sweeps,
    n_pirls = NA_integer_,
    n_criterion_evals = n_eval,
    history = history,
    penalties = penalties,
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = TRUE,
    method_lambda = method,
    optimizer = "ALS"
  )
}

#' Try Rcpp Gaussian cGCV / fixed-λ path.
#' @keywords internal
tt_als_fit_rcpp <- function(y, basis, ranks, lambda_spec, control, penalty_order = 2,
                            init_cores = NULL, offset = NULL) {
  offset <- normalize_offset(offset, length(y))
  if (any(offset != 0)) {
    # C++ path does not yet accept offsets
    out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                      init_cores = init_cores, offset = offset)
    out$backend <- "R"
    return(out)
  }
  d <- length(basis)
  p <- ncol(basis[[1]])
  lambda0 <- lambda_spec$values %||% lambda_spec$lambda0
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p, ranks[k + 1L], penalty_order)
  })
  t0 <- proc.time()[["elapsed"]]
  if (identical(lambda_spec$method, "cGCV") && exists("tt_cgcv_fit_cpp", mode = "function")) {
    fit <- tt_cgcv_fit_cpp(
      y = as.numeric(y),
      basis_list = basis,
      init_cores = cores,
      penalties_list = penalties,
      lambda_init = lambda0,
      sweeps = as.integer(control$max_sweeps),
      lambda_min = control$lambda_bounds[1],
      lambda_max = control$lambda_bounds[2],
      tol = control$tol_lambda,
      tol_lambda = control$tol_lambda,
      return_jacobian = FALSE
    )
    out_cores <- lapply(seq_len(d), function(k) {
      array(as.numeric(fit$cores[[k]]), c(ranks[k], p, ranks[k + 1L]))
    })
    return(list(
      cores = out_cores,
      intercept = fit$intercept,
      lambda = as.numeric(fit$lambda),
      ranks = ranks,
      eta = as.numeric(fit$mu),
      mu = as.numeric(fit$mu),
      deviance = sum((y - as.numeric(fit$mu))^2),
      n_sweeps = fit$n_sweeps,
      n_pirls = NA_integer_,
      n_criterion_evals = fit$n_criterion_evals,
      history = list(),
      penalties = penalties,
      elapsed = proc.time()[["elapsed"]] - t0,
      converged = TRUE,
      method_lambda = "cGCV",
      optimizer = "ALS",
      backend = "Rcpp"
    ))
  }
  if (identical(lambda_spec$method, "fixed") && exists("tt_fit_d_cpp", mode = "function")) {
    fit <- tt_fit_d_cpp(
      y = as.numeric(y),
      basis_list = basis,
      init_cores = cores,
      lambda = lambda0,
      penalties_list = penalties,
      sweeps = as.integer(control$max_sweeps),
      return_jacobian = FALSE
    )
    out_cores <- lapply(seq_len(d), function(k) {
      array(as.numeric(fit$cores[[k]]), c(ranks[k], p, ranks[k + 1L]))
    })
    return(list(
      cores = out_cores,
      intercept = fit$intercept,
      lambda = lambda0,
      ranks = ranks,
      eta = as.numeric(fit$mu),
      mu = as.numeric(fit$mu),
      deviance = sum((y - as.numeric(fit$mu))^2),
      n_sweeps = control$max_sweeps,
      n_pirls = NA_integer_,
      n_criterion_evals = 0L,
      history = list(),
      penalties = penalties,
      elapsed = proc.time()[["elapsed"]] - t0,
      converged = TRUE,
      method_lambda = "fixed",
      optimizer = "ALS",
      backend = "Rcpp"
    ))
  }
  # fallback
  out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                    init_cores = init_cores, offset = offset)
  out$backend <- "R"
  out
}
