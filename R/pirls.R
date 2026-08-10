#' GLM PIRLS + weighted TT-ALS (fixed or cGCV λ), R backend.
#'
#' Bernoulli uses consistent W/z working variance and optional true-objective
#' outer step-halving (see [tt_control()]).
#'
#' @keywords internal
tt_pirls_fit <- function(y, basis, family, ranks, lambda_spec, control,
                         penalty_order = 2, init_cores = NULL, offset = NULL,
                         weights = NULL) {
  method <- lambda_spec$method
  if (identical(method, "cGCV") &&
      identical(.cgcv_update_mode(control), "outer_simultaneous")) {
    return(tt_pirls_fit_cgcv_outer(
      y, basis, family, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores, offset = offset, weights = weights
    ))
  }
  tt_pirls_fit_sequential(
    y, basis, family, ranks, lambda_spec, control, penalty_order,
    init_cores = init_cores, offset = offset, weights = weights
  )
}

#' Sequential PIRLS + Gauss–Seidel cGCV (legacy dynamics).
#' @keywords internal
#' @noRd
tt_pirls_fit_sequential <- function(y, basis, family, ranks, lambda_spec, control,
                                    penalty_order = 2, init_cores = NULL,
                                    offset = NULL, weights = NULL) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  method <- lambda_spec$method
  lambda <- lambda_spec$values %||% lambda_spec$lambda0
  fam <- normalize_family(family)
  key <- family_key(fam)
  offset <- normalize_offset(offset, length(y))
  w_obs <- normalize_weights(weights, length(y))
  intercept <- init_intercept(fam, y, offset = offset, weights = w_obs)
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
    for (k in seq_len(d)) cores[[k]] <- cores[[k]] * 0.05
  } else {
    cores <- init_cores
  }
  penalties <- tt_core_penalties_from_basis(ranks, basis, penalty_order)
  penalty_mode <- "global"
  margin_order <- .cgcv_margin_order(control$cgcv_margin_order, d)
  rho <- control$cgcv_damping %||% 1
  delta <- control$cgcv_max_log10_step %||% Inf
  do_trace <- identical(method, "cGCV") && isTRUE(control$cgcv_trace %||% TRUE)
  cgcv_trace <- list()

  eta <- tt_eta(offset, intercept, cores, basis)
  mu <- invlink_eta(fam, eta)
  dev <- glm_deviance(fam, y, mu, weights = w_obs)
  obj <- tt_glm_penalized_objective(
    y, cores, intercept, basis, penalties, lambda, fam,
    offset = offset, weights = w_obs,
    penalty_mode = penalty_mode, penalty_order = penalty_order
  )

  hist_rows <- list()
  n_eval <- 0L
  n_pirls <- 0L
  n_als_sweeps_total <- 0L
  n_step_halvings_total <- 0L
  use_spec <- isTRUE(control$use_spectral_gcv)
  do_halving <- identical(key, "bernoulli") &&
    isTRUE(control$pirls_step_halving %||% control$damping %||% TRUE)
  step_factor <- control$step_factor %||% 0.5
  step_min <- control$step_min %||% (1 / 128)
  obj_tol <- control$objective_tol %||% 1e-10
  t0 <- proc.time()[["elapsed"]]

  pirls_ok <- TRUE
  reason <- "maxit"
  als_ok <- TRUE

  for (it in seq_len(control$pirls_maxit)) {
    cores_old <- cores
    intercept_old <- intercept
    eta_old <- eta
    obj_old <- obj
    work <- glm_working(fam, y, eta, control = control)
    w <- work$weight * w_obs
    z <- work$z

    # --- inner weighted ALS ---
    base_sw <- max(1L, as.integer(control$als_sweeps_per_pirls %||% 1L))
    max_sw <- if (isTRUE(control$als_sweeps_adaptive %||% TRUE)) {
      max(base_sw, as.integer(control$als_sweeps_per_pirls_max %||% 3L))
    } else {
      base_sw
    }
    for (sw in seq_len(max_sw)) {
      zc <- z - offset - intercept
      eta_before <- if (sw >= base_sw) tt_eta(offset, intercept, cores, basis) else NULL
      for (k in margin_order) {
        built <- .cgcv_core_workspace(
          cores, k, lambda, basis, zc, ranks, control,
          weight = w, penalty_order = penalty_order,
          use_spectral = identical(method, "cGCV") && isTRUE(control$use_spectral_gcv),
          compute_op_norms = do_trace
        )
        Pk <- built$P_own
        penalties[[k]] <- Pk
        ws <- built$workspace
        lambda_old_k <- lambda[k]
        gcv_old <- if (do_trace) .cgcv_eval_at(ws, lambda_old_k) else NULL
        upd <- update_lambda(method, ws)
        n_eval <- n_eval + upd$n_eval

        if (identical(method, "cGCV") &&
            (rho < 1 - 1e-15 || is.finite(delta))) {
          step <- .cgcv_damped_trust_update(
            lambda_old = lambda_old_k,
            lambda_tilde = upd$lambda,
            rho = rho,
            max_log10_step = delta,
            bounds = control$lambda_bounds
          )
          lam_new <- step$lambda_new[[1L]]
          fit_clip <- .cgcv_eval_at(ws, lam_new)
          g_use <- fit_clip$g
          ed_use <- fit_clip$ed
          gcv_use <- fit_clip$value
          rss_use <- fit_clip$rss
          tilde <- upd$lambda
        } else {
          lam_new <- upd$lambda
          g_use <- upd$g
          ed_use <- upd$ed
          gcv_use <- upd$value
          tilde <- upd$lambda
          rss_use <- if (identical(method, "cGCV")) {
            .cgcv_eval_at(ws, lam_new)$rss
          } else {
            NA_real_
          }
        }

        cores[[k]] <- array(g_use, c(ranks[k], p, ranks[k + 1L]))
        lambda[k] <- lam_new

        if (do_trace) {
          cgcv_trace[[length(cgcv_trace) + 1L]] <- data.frame(
            pirls = it,
            sweep = sw,
            margin = k,
            lambda_old = lambda_old_k,
            lambda_tilde = tilde,
            lambda_new = lam_new,
            log10_old = log10(lambda_old_k),
            log10_tilde = log10(tilde),
            log10_new = log10(lam_new),
            ed = ed_use,
            ed_old = if (is.null(gcv_old)) NA_real_ else gcv_old$ed,
            gcv = gcv_use,
            gcv_old = if (is.null(gcv_old)) NA_real_ else gcv_old$value,
            rss = rss_use,
            P_other_op = built$P_other_op,
            P_own_op = built$P_own_op,
            lambda_P_own_op = lam_new * built$P_own_op,
            boundary = .lambda_boundary_status(lam_new, control$lambda_bounds),
            mode = "sequential",
            stringsAsFactors = FALSE
          )
        }
      }
      f <- tt_contraction(cores, basis)
      intercept <- sum(w * (z - offset - f)) / max(sum(w), 1e-12)
      n_als_sweeps_total <- n_als_sweeps_total + 1L
      if (sw >= base_sw && !is.null(eta_before)) {
        eta_after <- tt_eta(offset, intercept, cores, basis)
        if (max(abs(eta_after - eta_before)) < control$tol) break
      }
    }

    cores_cand <- cores
    intercept_cand <- intercept
    eta_cand <- tt_eta(offset, intercept_cand, cores_cand, basis)
    obj_cand <- tt_glm_penalized_objective(
      y, cores_cand, intercept_cand, basis, penalties, lambda, fam,
      offset = offset, weights = w_obs,
      penalty_mode = penalty_mode, penalty_order = penalty_order
    )

    accepted_step <- 1
    n_halve_it <- 0L
    line_ok <- TRUE

    # Absolute + relative slack: near a plateau, ALS working-LS steps need not
    # be exact descent directions of the true Bernoulli objective.
    accept_tol <- max(obj_tol, control$tol * max(1, abs(obj_old$value)))

    if (do_halving) {
      accepted_step <- NA_real_
      alpha <- 1
      while (alpha + 1e-15 >= step_min) {
        blended <- tt_blend_params(
          cores_old, intercept_old, cores_cand, intercept_cand, alpha
        )
        obj_try <- tt_glm_penalized_objective(
          y, blended$cores, blended$intercept, basis, penalties, lambda, fam,
          offset = offset, weights = w_obs,
          penalty_mode = penalty_mode, penalty_order = penalty_order
        )
        if (is.finite(obj_try$value) &&
            obj_try$value <= obj_old$value + accept_tol) {
          cores <- blended$cores
          intercept <- blended$intercept
          eta <- obj_try$eta
          obj <- obj_try
          accepted_step <- alpha
          break
        }
        alpha <- alpha * step_factor
        n_halve_it <- n_halve_it + 1L
        n_step_halvings_total <- n_step_halvings_total + 1L
      }
      if (is.na(accepted_step)) {
        # restore last valid iterate
        cores <- cores_old
        intercept <- intercept_old
        eta <- eta_old
        obj <- obj_old
        line_ok <- FALSE
        n_pirls <- it
        mu <- invlink_eta(fam, eta)
        dev <- glm_deviance(fam, y, mu, weights = w_obs)
        cand_rel <- abs(obj_cand$value - obj_old$value) /
          max(1, abs(obj_old$value))
        # Late stall at a numerical stationary region of the PIRLS path:
        # keep last valid params; do not claim a healthy full step, but do
        # not treat plateau termination as catastrophic optimizer failure.
        if (it > 2L && is.finite(obj_old$value) &&
            (cand_rel < control$tol ||
               abs(obj_cand$value - obj_old$value) <= accept_tol)) {
          pirls_ok <- TRUE
          reason <- "stationary (no improving PIRLS step)"
        } else {
          pirls_ok <- FALSE
          reason <- "PIRLS line search failed"
        }
        hist_rows[[length(hist_rows) + 1L]] <- .pirls_hist_row(
          it, obj, y, eta, fam, work, control$als_sweeps_per_pirls,
          proposed_step = 1, accepted_step = 0, n_halve_it = n_halve_it,
          line_ok = FALSE, weights = w_obs
        )
        break
      }
    } else {
      # Gaussian / Poisson / Bernoulli without halving: accept candidate
      cores <- cores_cand
      intercept <- intercept_cand
      eta <- eta_cand
      obj <- obj_cand
      accepted_step <- 1
    }

    mu <- invlink_eta(fam, eta)
    dev <- glm_deviance(fam, y, mu, weights = w_obs)
    n_pirls <- it
    hist_rows[[length(hist_rows) + 1L]] <- .pirls_hist_row(
      it, obj, y, eta, fam, work, control$als_sweeps_per_pirls,
      proposed_step = 1, accepted_step = accepted_step, n_halve_it = n_halve_it,
      line_ok = TRUE, weights = w_obs
    )

    if (isTRUE(control$trace)) {
      cat(sprintf(
        "  PIRLS %2d | L=%.6g | dev=%.6g | step=%.4g | halvings=%d | max|eta|=%.3g\n",
        it, obj$value, dev, accepted_step, n_halve_it, max(abs(eta))
      ))
    }

    # Outer relative change on true objective (Bernoulli) or deviance (else)
    if (identical(key, "bernoulli")) {
      rel <- abs(obj_old$value - obj$value) / max(1, abs(obj_old$value))
    } else {
      rel <- abs(obj_old$value - obj$value) / max(1, abs(obj_old$value))
      # also track deviance for compatibility
      rel_dev <- abs(glm_deviance(fam, y, invlink_eta(fam, eta_old),
                                  weights = w_obs) - dev) /
        max(1, abs(dev))
      rel <- min(rel, rel_dev)
    }
    if (!is.finite(obj$value)) {
      pirls_ok <- FALSE
      reason <- "non-finite objective"
      break
    }
    if (it > 2L && rel < control$tol && accepted_step >= 1 - 1e-12) {
      reason <- "relative objective change"
      break
    }
    if (it == control$pirls_maxit) reason <- "maxit"
  }

  history <- if (length(hist_rows)) do.call(rbind, hist_rows) else
    data.frame(pirls = integer(), deviance = numeric())

  convergence <- list(
    overall = isTRUE(pirls_ok) && is.finite(obj$value),
    pirls = isTRUE(pirls_ok),
    als = isTRUE(als_ok),
    reason = reason,
    n_pirls = as.integer(n_pirls),
    n_als_sweeps = as.integer(n_als_sweeps_total),
    n_step_halvings = as.integer(n_step_halvings_total),
    gauge_note = paste(
      "Outer step-halving blends TT cores in parameter space without",
      "gauge alignment; R ALS does not re-orthogonalize between outer iterates."
    )
  )

  list(
    cores = cores,
    intercept = intercept,
    lambda = lambda,
    ranks = ranks,
    eta = eta,
    mu = invlink_eta(fam, eta),
    deviance = glm_deviance(fam, y, invlink_eta(fam, eta), weights = w_obs),
    n_sweeps = n_als_sweeps_total,
    n_pirls = n_pirls,
    n_criterion_evals = n_eval,
    history = history,
    penalties = penalties,
    penalty_mode = penalty_mode,
    cgcv = list(
      update = "sequential",
      parameterization = .cgcv_parameterization(control),
      damping = rho,
      max_log10_step = delta,
      margin_order = margin_order,
      trace = if (length(cgcv_trace)) do.call(rbind, cgcv_trace) else NULL,
      proposals = NULL
    ),
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = isTRUE(convergence$overall),
    convergence = convergence,
    method_lambda = method,
    optimizer = "ALS",
    backend = "R",
    objective = obj$value
  )
}

#' Outer simultaneous cGCV for GLM PIRLS.
#'
#' Fit all cores at fixed λ (PIRLS) → Jacobi proposals on last working
#' response → damped / trust-region λ update.
#' @keywords internal
#' @noRd
tt_pirls_fit_cgcv_outer <- function(y, basis, family, ranks, lambda_spec,
                                    control, penalty_order = 2,
                                    init_cores = NULL, offset = NULL,
                                    weights = NULL) {
  d <- length(basis)
  lambda <- as.numeric(lambda_spec$values %||% lambda_spec$lambda0)
  fam <- normalize_family(family)
  offset <- normalize_offset(offset, length(y))
  w_obs <- normalize_weights(weights, length(y))
  rho <- control$cgcv_damping %||% 1
  delta <- control$cgcv_max_log10_step %||% Inf
  param <- .cgcv_parameterization(control)
  outer_maxit <- as.integer(control$outer_maxit %||% 20L)
  if (isTRUE(rho >= 1 - 1e-15) && !is.finite(delta)) {
    warning(
      "cgcv_update='outer_simultaneous' with cgcv_damping=1 and ",
      "cgcv_max_log10_step=Inf is undamped Jacobi; consider ",
      "cgcv_damping=0.25 and cgcv_max_log10_step=1.",
      call. = FALSE
    )
  }

  # Warm start
  if (is.null(init_cores)) {
    p <- ncol(basis[[1]])
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
    for (k in seq_len(d)) cores[[k]] <- cores[[k]] * 0.05
  } else {
    cores <- init_cores
  }

  n_eval <- 0L
  history <- list()
  proposal_hist <- list()
  t0 <- proc.time()[["elapsed"]]
  prev_lam <- lambda
  n_outer <- 0L
  fit_last <- NULL

  for (outer in seq_len(outer_maxit)) {
    fixed_spec <- list(method = "fixed", values = lambda, automatic = FALSE,
                       lambda0 = lambda)
    fit_last <- tt_pirls_fit_sequential(
      y, basis, fam, ranks, fixed_spec, control, penalty_order,
      init_cores = cores, offset = offset, weights = w_obs
    )
    cores <- fit_last$cores
    intercept <- fit_last$intercept
    eta <- fit_last$eta
    work <- glm_working(fam, y, eta, control = control)
    w <- work$weight * w_obs
    zc <- work$z - offset - intercept

    step <- .cgcv_simultaneous_step(
      cores, lambda, basis, zc, ranks, control,
      weight = w, penalty_order = penalty_order
    )
    n_eval <- n_eval + step$n_eval + fit_last$n_criterion_evals
    lambda <- step$lambda
    n_outer <- outer

    history[[outer]] <- list(
      outer = outer,
      deviance = fit_last$deviance,
      objective = fit_last$objective,
      lambda = lambda
    )
    prop_df <- step$proposals
    prop_df$outer <- outer
    proposal_hist[[outer]] <- prop_df

    if (isTRUE(control$trace)) {
      cat(sprintf(
        "  cGCV-outer %2d | dev=%.6g | λ=%s\n",
        outer, fit_last$deviance,
        paste(sprintf("%.3g", lambda), collapse = ",")
      ))
    }

    dlog <- max(abs(log(lambda) - log(pmax(prev_lam, 1e-12))))
    if (dlog < control$tol_lambda && outer > 1L) break
    prev_lam <- lambda
  }

  fixed_spec <- list(method = "fixed", values = lambda, automatic = FALSE,
                     lambda0 = lambda)
  fit_final <- tt_pirls_fit_sequential(
    y, basis, fam, ranks, fixed_spec, control, penalty_order,
    init_cores = cores, offset = offset, weights = w_obs
  )

  lambda0_table <- NULL
  if (identical(param, "scale_anisotropy") &&
      identical(control$cgcv_lambda0_method %||% "fixed_start", "log_grid")) {
    grid <- control$cgcv_lambda0_grid
    if (is.null(grid)) {
      b <- control$lambda_bounds
      grid <- exp(seq(log(b[1]), log(b[2]), length.out = 9L))
    }
    sa <- .cgcv_scale_anisotropy_from_lambda(fit_final$lambda)
    sel <- .cgcv_select_lambda0_grid(sa$omega, grid, fit_fixed = function(lam) {
      fs <- list(method = "fixed", values = lam, automatic = FALSE, lambda0 = lam)
      ff <- tt_pirls_fit_sequential(
        y, basis, fam, ranks, fs, control, penalty_order,
        init_cores = fit_final$cores, offset = offset, weights = w_obs
      )
      list(criterion = ff$deviance, fit = ff)
    })
    lambda0_table <- sel$table
    fs <- list(method = "fixed", values = sel$lambda, automatic = FALSE,
               lambda0 = sel$lambda)
    fit_final <- tt_pirls_fit_sequential(
      y, basis, fam, ranks, fs, control, penalty_order,
      init_cores = fit_final$cores, offset = offset, weights = w_obs
    )
  }

  prop_all <- if (length(proposal_hist)) do.call(rbind, proposal_hist) else NULL
  fit_final$n_outer <- n_outer
  fit_final$n_criterion_evals <- n_eval + fit_final$n_criterion_evals
  fit_final$method_lambda <- "cGCV"
  fit_final$history_outer <- history
  fit_final$cgcv <- list(
    update = "outer_simultaneous",
    parameterization = param,
    damping = rho,
    max_log10_step = delta,
    proposals = prop_all,
    lambda0_table = lambda0_table,
    trace = prop_all
  )
  fit_final$elapsed <- proc.time()[["elapsed"]] - t0
  fit_final
}

.pirls_hist_row <- function(it, obj, y, eta, family, work, als_sweeps,
                            proposed_step, accepted_step, n_halve_it, line_ok,
                            weights = NULL) {
  mu <- invlink_eta(family, eta)
  data.frame(
    pirls = it,
    objective = obj$value,
    nll = obj$nll,
    penalty = obj$penalty,
    deviance = glm_deviance(family, y, mu, weights = weights),
    max_abs_eta = max(abs(eta)),
    min_probability = min(mu),
    max_probability = max(mu),
    min_weight = min(work$weight),
    median_weight = stats::median(work$weight),
    als_sweeps = als_sweeps,
    proposed_step = proposed_step,
    accepted_step = accepted_step,
    n_step_halvings = n_halve_it,
    line_ok = line_ok,
    stringsAsFactors = FALSE
  )
}

#' GLM PIRLS via Rcpp when available.
#'
#' Bernoulli with true-objective step-halving uses the R path (FIX 1–2).
#'
#' @keywords internal
tt_pirls_fit_rcpp <- function(y, basis, family, ranks, lambda_spec, control,
                              penalty_order = 2, init_cores = NULL, offset = NULL,
                              weights = NULL) {
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  key <- family_key(family)
  do_halving <- identical(key, "bernoulli") &&
    isTRUE(control$pirls_step_halving %||% control$damping %||% TRUE)
  if (do_halving || !exists("tt_glm_pirls_cgcv_cpp", mode = "function")) {
    return(tt_pirls_fit(y, basis, family, ranks, lambda_spec, control,
                        penalty_order, init_cores = init_cores,
                        offset = offset, weights = w))
  }
  d <- length(basis)
  p <- ncol(basis[[1]])
  lambda0 <- lambda_spec$values %||% lambda_spec$lambda0
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- tt_core_penalties_from_basis(ranks, basis, penalty_order)
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
    select_lambda = identical(lambda_spec$method, "cGCV"),
    weights = w,
    offset = offset
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
    convergence = list(
      overall = is.finite(fit$deviance),
      pirls = is.finite(fit$deviance),
      als = TRUE,
      reason = "rcpp_legacy",
      n_pirls = fit$n_pirls,
      n_als_sweeps = NA_integer_,
      n_step_halvings = 0L
    ),
    method_lambda = lambda_spec$method,
    optimizer = "ALS",
    backend = "Rcpp"
  )
}
