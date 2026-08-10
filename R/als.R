#' Gaussian TT-ALS with fixed or cGCV λ (R backend).
#'
#' Always uses the classical global discrete P-spline penalty on Θ via
#' \(P_k^{\mathrm{full}}=A_k^\top S_{\boldsymbol\lambda} A_k\).
#' @keywords internal
tt_als_fit <- function(y, basis, ranks, lambda_spec, control, penalty_order = 2,
                       init_cores = NULL, offset = NULL, weights = NULL) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  method <- lambda_spec$method
  lambda <- lambda_spec$values %||% lambda_spec$lambda0
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  intercept <- sum(w * (y - offset)) / sum(w)
  yc <- y - offset - intercept
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- tt_core_penalties_from_basis(ranks, basis, penalty_order)
  cyclic <- attr(basis, "cyclic")
  # Per-core Q descent for fixed λ (P6)
  check_q_descent <- identical(method, "fixed")
  q_tol <- as.numeric(control$objective_tol %||% 1e-10)
  q_max_increase <- 0
  q_violations <- 0L

  n_eval <- 0L
  n_sweeps <- 0L
  history <- list()
  prev_lam <- lambda
  prev_eta <- NULL
  t0 <- proc.time()[["elapsed"]]
  use_spec <- isTRUE(control$use_spectral_gcv)

  for (sw in seq_len(control$max_sweeps)) {
    for (k in seq_len(d)) {
      if (check_q_descent) {
        q_old <- tt_gaussian_Q(
          y, cores, intercept, basis, lambda,
          offset = offset, weights = w,
          penalty_order = penalty_order, cyclic = cyclic
        )$value
      }
      L <- left_interfaces(cores, basis)
      R <- right_interfaces(cores, basis)
      Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
      pen_k <- tt_conditional_penalty_full(
        cores, k, lambda, penalty_order = penalty_order, cyclic = cyclic
      )
      Pk <- pen_k$P_own
      P0 <- pen_k$P_other
      penalties[[k]] <- Pk
      ws <- make_core_workspace(
        yc, Xk, Pk, lambda[k],
        control$lambda_bounds, control$tol_lambda,
        weight = w, use_spectral = use_spec, P0 = P0
      )
      upd <- update_lambda(method, ws)
      cores[[k]] <- array(upd$g, c(ranks[k], p, ranks[k + 1L]))
      lambda[k] <- upd$lambda
      n_eval <- n_eval + upd$n_eval
      if (check_q_descent) {
        q_new <- tt_gaussian_Q(
          y, cores, intercept, basis, lambda,
          offset = offset, weights = w,
          penalty_order = penalty_order, cyclic = cyclic
        )$value
        dq <- q_new - q_old
        if (is.finite(dq) && dq > q_max_increase) q_max_increase <- dq
        if (is.finite(dq) && dq > q_tol * max(1, abs(q_old))) {
          q_violations <- q_violations + 1L
        }
      }
    }
    n_sweeps <- sw
    eta <- tt_eta(offset, intercept, cores, basis)
    rss <- sum(w * (y - eta)^2)
    pen_val <- tt_global_penalty_value(
      cores, lambda, penalty_order = penalty_order, cyclic = cyclic
    )
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
    deviance = sum(w * (y - eta)^2),
    n_sweeps = n_sweeps,
    n_pirls = NA_integer_,
    n_criterion_evals = n_eval,
    history = history,
    penalties = penalties,
    penalty_mode = "global",
    q_descent = if (check_q_descent) {
      list(
        checked = TRUE,
        violations = q_violations,
        max_increase = q_max_increase,
        tol = q_tol
      )
    } else {
      list(checked = FALSE)
    },
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = TRUE,
    method_lambda = method,
    optimizer = "ALS"
  )
}

#' Rcpp entry for Gaussian ALS — redirects to the global-penalty R path.
#'
#' Legacy full-sweep Rcpp ALS used the own-margin surrogate and is no longer
#' a production fitter. Penalty helpers remain available via Rcpp.
#' @keywords internal
tt_als_fit_rcpp <- function(y, basis, ranks, lambda_spec, control, penalty_order = 2,
                            init_cores = NULL, offset = NULL, weights = NULL) {
  out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                    init_cores = init_cores, offset = offset, weights = weights)
  out$backend <- "R"
  out
}
