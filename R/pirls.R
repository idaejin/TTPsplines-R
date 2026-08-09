#' GLM PIRLS + weighted TT-ALS (fixed or cGCV λ), R backend.
#' @keywords internal
tt_pirls_fit <- function(y, basis, family, ranks, lambda_spec, control,
                         penalty_order = 2, init_cores = NULL) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  method <- lambda_spec$method
  lambda <- lambda_spec$values %||% lambda_spec$lambda0
  fam <- normalize_family(family)
  intercept <- init_intercept(fam, y)
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
    for (k in seq_len(d)) cores[[k]] <- cores[[k]] * 0.05
  } else {
    cores <- init_cores
  }
  penalties <- lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p, ranks[k + 1L], penalty_order)
  })

  eta <- intercept + tt_contraction(cores, basis)
  mu <- invlink_eta(fam, eta)
  dev <- glm_deviance(fam, y, mu)
  history <- data.frame(pirls = integer(), deviance = numeric())
  n_eval <- 0L
  t0 <- proc.time()[["elapsed"]]
  n_pirls <- 0L
  use_spec <- isTRUE(control$use_spectral_gcv)

  for (it in seq_len(control$pirls_maxit)) {
    work <- glm_working(fam, y, eta)
    w <- work$weight
    z <- work$z

    for (sw in seq_len(control$als_sweeps_per_pirls)) {
      zc <- z - intercept
      for (k in seq_len(d)) {
        L <- left_interfaces(cores, basis)
        R <- right_interfaces(cores, basis)
        Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
        ws <- make_core_workspace(
          zc, Xk, penalties[[k]], lambda[k],
          control$lambda_bounds, control$tol_lambda,
          weight = w, use_spectral = use_spec
        )
        upd <- update_lambda(method, ws)
        cores[[k]] <- array(upd$g, c(ranks[k], p, ranks[k + 1L]))
        lambda[k] <- upd$lambda
        n_eval <- n_eval + upd$n_eval
      }
      f <- tt_contraction(cores, basis)
      intercept <- sum(w * (z - f)) / max(sum(w), 1e-12)
    }

    eta_cand <- intercept + tt_contraction(cores, basis)
    mu_cand <- invlink_eta(fam, eta_cand)
    dev_cand <- glm_deviance(fam, y, mu_cand)

    # Soft deviance guard for Bernoulli
    if (isTRUE(control$damping) && identical(family_key(fam), "bernoulli") &&
        is.finite(dev) && is.finite(dev_cand) && dev_cand > dev * 1.25 && it > 1L) {
      if (control$trace) {
        cat(sprintf("  PIRLS %2d | damped (dev %.4g → %.4g)\n", it, dev, dev_cand))
      }
      history <- rbind(history, data.frame(pirls = it, deviance = dev))
      n_pirls <- it
      break
    }

    eta <- eta_cand
    mu <- mu_cand
    n_pirls <- it
    history <- rbind(history, data.frame(pirls = it, deviance = dev_cand))
    if (control$trace) {
      cat(sprintf("  PIRLS %2d | dev=%.6g | λ=%s\n",
                  it, dev_cand, paste(sprintf("%.3g", lambda), collapse = ",")))
    }
    rel <- abs(dev - dev_cand) / max(1, abs(dev))
    if (!is.finite(dev_cand)) break
    if (it > 2L && rel < control$tol) {
      dev <- dev_cand
      break
    }
    dev <- dev_cand
  }

  list(
    cores = cores,
    intercept = intercept,
    lambda = lambda,
    ranks = ranks,
    eta = eta,
    mu = mu,
    deviance = dev,
    n_sweeps = NA_integer_,
    n_pirls = n_pirls,
    n_criterion_evals = n_eval,
    history = history,
    penalties = penalties,
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = is.finite(dev),
    method_lambda = method,
    optimizer = "ALS",
    backend = "R"
  )
}

#' GLM PIRLS via Rcpp when available.
#' @keywords internal
tt_pirls_fit_rcpp <- function(y, basis, family, ranks, lambda_spec, control,
                              penalty_order = 2, init_cores = NULL) {
  key <- family_key(family)
  if (!exists("tt_glm_pirls_cgcv_cpp", mode = "function")) {
    return(tt_pirls_fit(y, basis, family, ranks, lambda_spec, control,
                        penalty_order, init_cores = init_cores))
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
  fit <- tt_glm_pirls_cgcv_cpp(
    y = as.numeric(y),
    basis_list = basis,
    init_cores = cores,
    penalties_list = penalties,
    family = key,
    lambda_init = lambda0,
    pirls_iter = as.integer(control$pirls_maxit),
    als_sweeps = as.integer(control$als_sweeps_per_pirls),
    lambda_min = control$lambda_bounds[1],
    lambda_max = control$lambda_bounds[2],
    tol = control$tol_lambda,
    tol_dev = control$tol,
    select_lambda = identical(lambda_spec$method, "cGCV")
  )
  out_cores <- lapply(seq_len(d), function(k) {
    array(as.numeric(fit$cores[[k]]), c(ranks[k], p, ranks[k + 1L]))
  })
  list(
    cores = out_cores,
    intercept = fit$intercept,
    lambda = as.numeric(fit$lambda),
    ranks = ranks,
    eta = as.numeric(fit$eta),
    mu = as.numeric(fit$mu),
    deviance = fit$deviance,
    n_sweeps = NA_integer_,
    n_pirls = fit$n_pirls,
    n_criterion_evals = fit$n_criterion_evals,
    history = as.data.frame(fit$history),
    penalties = penalties,
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = is.finite(fit$deviance),
    method_lambda = lambda_spec$method,
    optimizer = "ALS",
    backend = "Rcpp"
  )
}
