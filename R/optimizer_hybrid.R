#' Experimental hybrid: PIRLS+ALS warm-start then short L-BFGS polish.
#'
#' Uses the same Bernoulli (or GLM) penalized objective as ALS and L-BFGS.
#' Intended when ALS reaches poor-predicting stationary regions.
#'
#' @keywords internal
tt_hybrid_fit <- function(y, basis, family, ranks, lambda_spec, control,
                          penalty_order = 2, init_cores = NULL) {
  t0 <- proc.time()[["elapsed"]]
  fam <- normalize_family(family)
  key <- family_key(fam)

  als <- if (identical(key, "gaussian")) {
    tt_als_fit(
      y, basis, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores
    )
  } else {
    tt_pirls_fit(
      y, basis, fam, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores
    )
  }

  polish_ctrl <- control
  polish_ctrl$lbfgs_maxit <- as.integer(control$hybrid_lbfgs_maxit %||% 50L)

  lbf <- tt_lbfgs_fit(
    y, basis, ranks, lambda_spec, polish_ctrl, penalty_order,
    init_cores = als$cores,
    family = if (identical(key, "gaussian")) NULL else fam,
    intercept0 = als$intercept
  )

  conv_als <- als$convergence %||% list(
    overall = isTRUE(als$converged),
    pirls = NA, als = TRUE,
    reason = "ALS",
    n_pirls = als$n_pirls, n_als_sweeps = als$n_sweeps,
    n_step_halvings = 0L
  )
  conv_lbf <- lbf$convergence %||% list(
    overall = isTRUE(lbf$converged),
    reason = "L-BFGS polish"
  )

  list(
    cores = lbf$cores,
    intercept = lbf$intercept,
    lambda = als$lambda,
    ranks = ranks,
    eta = lbf$eta,
    mu = lbf$mu,
    deviance = lbf$deviance,
    n_sweeps = als$n_sweeps,
    n_pirls = als$n_pirls,
    n_criterion_evals = (als$n_criterion_evals %||% 0L) +
      (lbf$n_criterion_evals %||% 0L),
    n_opt_iter = lbf$n_opt_iter,
    n_outer = lbf$n_outer %||% 1L,
    history = als$history,
    penalties = als$penalties %||% lbf$penalties,
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = isTRUE(lbf$converged),
    convergence = list(
      overall = isTRUE(lbf$converged),
      pirls = conv_als$pirls,
      als = conv_als$als %||% TRUE,
      lbfgs = isTRUE(lbf$converged),
      reason = paste0(
        "hybrid: ALS[", conv_als$reason %||% "?", "] -> LBFGS[",
        conv_lbf$reason %||% "polish", "]"
      ),
      n_pirls = conv_als$n_pirls %||% als$n_pirls,
      n_als_sweeps = conv_als$n_als_sweeps %||% als$n_sweeps,
      n_step_halvings = conv_als$n_step_halvings %||% 0L,
      n_opt_iter = lbf$n_opt_iter,
      hybrid_lbfgs_maxit = polish_ctrl$lbfgs_maxit
    ),
    method_lambda = als$method_lambda %||% lambda_spec$method,
    optimizer = "hybrid",
    backend = "R"
  )
}
