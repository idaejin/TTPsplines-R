#' Gaussian TT-ALS with fixed or cGCV λ (R backend).
#' @keywords internal
tt_als_fit <- function(y, basis, ranks, lambda_spec, control, penalty_order = 2) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  method <- lambda_spec$method
  lambda <- lambda_spec$lambda0
  intercept <- mean(y)
  yc <- y - intercept
  cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  penalties <- lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p, ranks[k + 1L], penalty_order)
  })

  n_eval <- 0L
  n_sweeps <- 0L
  history <- list()
  prev_lam <- lambda
  t0 <- proc.time()[["elapsed"]]

  for (sw in seq_len(control$max_sweeps)) {
    for (k in seq_len(d)) {
      L <- left_interfaces(cores, basis)
      R <- right_interfaces(cores, basis)
      Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
      ws <- make_core_workspace(
        yc, Xk, penalties[[k]], lambda[k],
        control$lambda_bounds, control$tol_lambda
      )
      upd <- update_lambda(method, ws)
      cores[[k]] <- array(upd$g, c(ranks[k], p, ranks[k + 1L]))
      lambda[k] <- upd$lambda
      n_eval <- n_eval + upd$n_eval
    }
    n_sweeps <- sw
    eta <- intercept + tt_contraction(cores, basis)
    rss <- sum((y - eta)^2)
    history[[sw]] <- list(sweep = sw, rss = rss, lambda = lambda)
    if (control$trace) {
      cat(sprintf("  ALS sweep %2d | RSS=%.6g | λ=%s\n",
                  sw, rss, paste(sprintf("%.3g", lambda), collapse = ",")))
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

  eta <- intercept + tt_contraction(cores, basis)
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
    method_lambda = method
  )
}

#' Try Rcpp Gaussian cGCV / fixed-λ path.
#' @keywords internal
tt_als_fit_rcpp <- function(y, basis, ranks, lambda_spec, control, penalty_order = 2) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
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
      lambda_init = lambda_spec$lambda0,
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
      backend = "Rcpp"
    ))
  }
  if (identical(lambda_spec$method, "fixed") && exists("tt_fit_d_cpp", mode = "function")) {
    fit <- tt_fit_d_cpp(
      y = as.numeric(y),
      basis_list = basis,
      init_cores = cores,
      lambda = lambda_spec$lambda0,
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
      lambda = lambda_spec$lambda0,
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
      backend = "Rcpp"
    ))
  }
  # fallback
  out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order)
  out$backend <- "R"
  out
}
