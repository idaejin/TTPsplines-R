#' Gaussian TT-ALS with fixed or cGCV λ (R backend).
#'
#' Always uses the classical global discrete P-spline penalty on Θ via
#' \(P_k^{\mathrm{full}}=A_k^\top S_{\boldsymbol\lambda} A_k\).
#'
#' cGCV modes ([tt_control()] `cgcv_update`):
#' - `"sequential"` (default): Gauss–Seidel core/λ updates (legacy dynamics).
#' - `"outer_simultaneous"`: fit all cores at fixed λ → freeze → Jacobi
#'   proposals → damped / trust-region λ update.
#' @keywords internal
tt_als_fit <- function(y, basis, ranks, lambda_spec, control, penalty_order = 2,
                       init_cores = NULL, offset = NULL, weights = NULL,
                       linear = NULL, smooth = NULL) {
  method <- lambda_spec$method
  if (identical(method, "cGCV") &&
      identical(.cgcv_update_mode(control), "outer_simultaneous")) {
    return(tt_als_fit_cgcv_outer(
      y, basis, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores, offset = offset, weights = weights,
      linear = linear, smooth = smooth
    ))
  }
  tt_als_fit_sequential(
    y, basis, ranks, lambda_spec, control, penalty_order,
    init_cores = init_cores, offset = offset, weights = weights,
    linear = linear, smooth = smooth
  )
}

#' Sequential (Gauss–Seidel) ALS — fixed λ or immediate cGCV.
#' @keywords internal
#' @noRd
tt_als_fit_sequential <- function(y, basis, ranks, lambda_spec, control,
                                  penalty_order = 2, init_cores = NULL,
                                  offset = NULL, weights = NULL,
                                  linear = NULL, smooth = NULL) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  method <- lambda_spec$method
  lambda <- lambda_spec$values %||% lambda_spec$lambda0
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  linear <- normalize_linear(linear, length(y))
  smooth <- normalize_smooth(smooth, length(y))
  ab0 <- tt_update_intercept_beta(y, offset, f = 0, linear = linear, weights = w)
  intercept <- ab0$intercept
  beta <- ab0$beta
  if (!is.null(smooth)) {
    add <- tt_refresh_additive(y, offset, f_tt = 0, linear = linear,
                               smooth = smooth, weights = w, control = control)
    intercept <- add$intercept
    beta <- add$beta
    smooth <- add$smooth
  }
  yc <- y - offset - intercept - tt_linear_contrib(linear, beta) -
    tt_smooth_contrib(smooth)
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- tt_core_penalties_from_basis(ranks, basis, penalty_order)
  cyclic <- attr(basis, "cyclic")
  check_q_descent <- identical(method, "fixed")
  q_tol <- as.numeric(control$objective_tol %||% 1e-10)
  q_max_increase <- 0
  q_violations <- 0L
  do_trace <- identical(method, "cGCV") && isTRUE(control$cgcv_trace %||% TRUE)
  margin_order <- .cgcv_margin_order(control$cgcv_margin_order, d)
  rho <- control$cgcv_damping %||% 1
  delta <- control$cgcv_max_log10_step %||% Inf
  bounds <- control$lambda_bounds

  n_eval <- 0L
  n_sweeps <- 0L
  history <- list()
  cgcv_trace <- list()
  prev_lam <- lambda
  prev_eta <- NULL
  t0 <- proc.time()[["elapsed"]]
  use_spec <- isTRUE(control$use_spectral_gcv)

  for (sw in seq_len(control$max_sweeps)) {
    use_cache <- isTRUE(control$design_interface_cache %||% TRUE)
    use_ltr <- use_cache && .tt_is_ltr_order(margin_order, d)
    use_rtl <- use_cache && .tt_is_rtl_order(margin_order, d)
    if (use_ltr) {
      R_all <- .tt_design_prepare_right(cores, basis)
      L_cur <- matrix(1, nrow(basis[[1]]), 1)
    } else if (use_rtl) {
      L_all <- .tt_design_prepare_left(cores, basis)
      R_cur <- matrix(1, nrow(basis[[1]]), 1)
    }
    for (k in margin_order) {
      if (check_q_descent) {
        q_old <- tt_gaussian_Q(
          y, cores, intercept, basis, lambda,
          offset = offset, weights = w,
          penalty_order = penalty_order, cyclic = cyclic,
          linear = linear, beta = beta, smooth = smooth
        )$value
      }
      if (use_ltr) {
        Left <- L_cur
        Right <- R_all[[k]]
      } else if (use_rtl) {
        Left <- L_all[[k]]
        Right <- R_cur
      } else {
        Left <- NULL
        Right <- NULL
      }
      built <- .cgcv_core_workspace(
        cores, k, lambda, basis, yc, ranks, control,
        weight = w, penalty_order = penalty_order,
        use_spectral = identical(method, "cGCV") && isTRUE(control$use_spectral_gcv),
        compute_op_norms = do_trace,
        Left = Left,
        Right = Right
      )
      Pk <- built$P_own
      penalties[[k]] <- Pk
      ws <- built$workspace
      lambda_old_k <- lambda[k]
      gcv_old <- if (do_trace && identical(method, "cGCV")) {
        .cgcv_eval_at(ws, lambda_old_k)
      } else {
        NULL
      }
      upd <- update_lambda(method, ws)
      n_eval <- n_eval + upd$n_eval

      # Optional damping / trust even in sequential mode (rho=1, Inf = no-op)
      if (identical(method, "cGCV") && (rho < 1 - 1e-15 || is.finite(delta))) {
        step <- .cgcv_damped_trust_update(
          lambda_old = lambda_old_k,
          lambda_tilde = upd$lambda,
          rho = rho,
          max_log10_step = delta,
          bounds = bounds
        )
        lam_new <- step$lambda_new[[1L]]
        # Refit core at clipped λ
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
        rss_use <- NA_real_
        tilde <- upd$lambda
        if (identical(method, "cGCV")) {
          fit_tmp <- .cgcv_eval_at(ws, lam_new)
          rss_use <- fit_tmp$rss
        }
      }

      cores[[k]] <- array(g_use, c(ranks[k], p, ranks[k + 1L]))
      lambda[k] <- lam_new

      if (use_ltr && k < d) {
        L_cur <- .tt_design_left_absorb(L_cur, cores[[k]], basis[[k]])
      } else if (use_rtl && k > 1L) {
        R_cur <- .tt_design_right_absorb(R_cur, cores[[k]], basis[[k]])
      }

      if (do_trace) {
        cgcv_trace[[length(cgcv_trace) + 1L]] <- data.frame(
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
          boundary = .lambda_boundary_status(lam_new, bounds),
          mode = "sequential",
          stringsAsFactors = FALSE
        )
      }

      if (check_q_descent) {
        q_new <- tt_gaussian_Q(
          y, cores, intercept, basis, lambda,
          offset = offset, weights = w,
          penalty_order = penalty_order, cyclic = cyclic,
          linear = linear, beta = beta, smooth = smooth
        )$value
        dq <- q_new - q_old
        if (is.finite(dq) && dq > q_max_increase) q_max_increase <- dq
        if (is.finite(dq) && dq > q_tol * max(1, abs(q_old))) {
          q_violations <- q_violations + 1L
        }
      }
    }
    # Refresh additive (smooth + α, β) given current TT surface; update residual
    f <- tt_contraction(cores, basis)
    add <- tt_refresh_additive(y, offset, f, linear = linear, smooth = smooth,
                               weights = w, control = control)
    intercept <- add$intercept
    beta <- add$beta
    smooth <- add$smooth
    yc <- y - offset - intercept - tt_linear_contrib(linear, beta) -
      tt_smooth_contrib(smooth)

    n_sweeps <- sw
    eta <- tt_eta(offset, intercept, cores, basis,
                  linear = linear, beta = beta, smooth = smooth)
    rss <- sum(w * (y - eta)^2)
    pen_val <- tt_global_penalty_value(
      cores, lambda, penalty_order = penalty_order, cyclic = cyclic
    ) + tt_smooth_penalty_value(smooth)
    obj <- 0.5 * rss + pen_val
    d_eta <- if (sw == 1L) NA_real_ else sqrt(mean((eta - prev_eta)^2))
    history[[sw]] <- list(
      sweep = sw, rss = rss, objective = obj, penalty = pen_val,
      lambda = lambda, d_eta = d_eta
    )
    prev_eta <- eta
    if (control$trace) {
      cat(sprintf("  ALS sweep %2d | obj=%.6g | RSS=%.6g | lambda=%s\n",
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

  eta <- tt_eta(offset, intercept, cores, basis,
                linear = linear, beta = beta, smooth = smooth)
  cgcv_df <- if (length(cgcv_trace)) do.call(rbind, cgcv_trace) else NULL
  list(
    cores = cores,
    intercept = intercept,
    beta = beta,
    linear = linear,
    smooth = smooth,
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
    cgcv = list(
      update = "sequential",
      parameterization = .cgcv_parameterization(control),
      damping = rho,
      max_log10_step = delta,
      margin_order = margin_order,
      trace = cgcv_df,
      proposals = NULL
    ),
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

#' Outer simultaneous cGCV for Gaussian ALS.
#'
#' Fit all cores at fixed λ → Jacobi proposals → damped / trust update.
#' @keywords internal
#' @noRd
tt_als_fit_cgcv_outer <- function(y, basis, ranks, lambda_spec, control,
                                  penalty_order = 2, init_cores = NULL,
                                  offset = NULL, weights = NULL,
                                  linear = NULL, smooth = NULL) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  lambda <- as.numeric(lambda_spec$values %||% lambda_spec$lambda0)
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  linear <- normalize_linear(linear, length(y))
  smooth <- normalize_smooth(smooth, length(y))
  ab0 <- tt_update_intercept_beta(y, offset, f = 0, linear = linear, weights = w)
  intercept <- ab0$intercept
  beta <- ab0$beta
  if (!is.null(smooth)) {
    add <- tt_refresh_additive(y, offset, f_tt = 0, linear = linear,
                               smooth = smooth, weights = w, control = control)
    intercept <- add$intercept
    beta <- add$beta
    smooth <- add$smooth
  }
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  cyclic <- attr(basis, "cyclic")
  outer_maxit <- as.integer(control$outer_maxit %||% 20L)
  fit_sweeps <- as.integer(control$cgcv_fit_sweeps %||% control$max_sweeps)
  rho <- control$cgcv_damping %||% 1
  delta <- control$cgcv_max_log10_step %||% Inf
  if (isTRUE(rho >= 1 - 1e-15) && !is.finite(delta)) {
    warning(
      "cgcv_update='outer_simultaneous' with cgcv_damping=1 and ",
      "cgcv_max_log10_step=Inf is undamped Jacobi; consider ",
      "cgcv_damping=0.25 and cgcv_max_log10_step=1.",
      call. = FALSE
    )
  }

  # Optional overall-scale grid for scale–anisotropy
  param <- .cgcv_parameterization(control)
  n_eval <- 0L
  history <- list()
  proposal_hist <- list()
  t0 <- proc.time()[["elapsed"]]
  prev_lam <- lambda
  n_outer <- 0L

  for (outer in seq_len(outer_maxit)) {
    # ---- A. Fit all cores at fixed λ ----
    fixed_spec <- list(method = "fixed", values = lambda, automatic = FALSE,
                       lambda0 = lambda)
    ctrl_fit <- control
    ctrl_fit$max_sweeps <- fit_sweeps
    fit_fixed <- tt_als_fit_sequential(
      y, basis, ranks, fixed_spec, ctrl_fit, penalty_order,
      init_cores = cores, offset = offset, weights = w,
      linear = linear, smooth = smooth
    )
    cores <- fit_fixed$cores
    intercept <- fit_fixed$intercept
    beta <- fit_fixed$beta
    smooth <- fit_fixed$smooth
    penalties <- fit_fixed$penalties

    # ---- B/C. Frozen Jacobi proposals + damp/trust ----
    yc <- y - offset - intercept - tt_linear_contrib(linear, beta) -
      tt_smooth_contrib(smooth)
    step <- .cgcv_simultaneous_step(
      cores, lambda, basis, yc, ranks, control,
      weight = w, penalty_order = penalty_order
    )
    n_eval <- n_eval + step$n_eval + fit_fixed$n_criterion_evals
    lambda <- step$lambda
    n_outer <- outer

    eta <- tt_eta(offset, intercept, cores, basis,
                  linear = linear, beta = beta, smooth = smooth)
    rss <- sum(w * (y - eta)^2)
    pen_val <- tt_global_penalty_value(
      cores, lambda, penalty_order = penalty_order, cyclic = cyclic
    ) + tt_smooth_penalty_value(smooth)
    obj <- 0.5 * rss + pen_val
    history[[outer]] <- list(
      outer = outer, rss = rss, objective = obj, penalty = pen_val,
      lambda = lambda, deviance = rss
    )
    prop_df <- step$proposals
    prop_df$outer <- outer
    proposal_hist[[outer]] <- prop_df

    if (isTRUE(control$trace)) {
      cat(sprintf(
        "  cGCV-outer %2d | obj=%.6g | lambda=%s\n",
        outer, obj, paste(sprintf("%.3g", lambda), collapse = ",")
      ))
    }

    dlog <- max(abs(log(lambda) - log(pmax(prev_lam, 1e-12))))
    if (dlog < control$tol_lambda && outer > 1L) break
    prev_lam <- lambda
  }

  # Final polish at converged λ
  fixed_spec <- list(method = "fixed", values = lambda, automatic = FALSE,
                     lambda0 = lambda)
  fit_final <- tt_als_fit_sequential(
    y, basis, ranks, fixed_spec, control, penalty_order,
    init_cores = cores, offset = offset, weights = w,
    linear = linear, smooth = smooth
  )

  # Optional λ0 grid after anisotropy settled
  lambda0_table <- NULL
  if (identical(param, "scale_anisotropy") &&
      identical(control$cgcv_lambda0_method %||% "fixed_start", "log_grid")) {
    grid <- control$cgcv_lambda0_grid
    if (is.null(grid)) {
      b <- control$lambda_bounds
      grid <- exp(seq(log(b[1]), log(b[2]), length.out = 11L))
    }
    sa <- .cgcv_scale_anisotropy_from_lambda(fit_final$lambda)
    sel <- .cgcv_select_lambda0_grid(sa$omega, grid, fit_fixed = function(lam) {
      fs <- list(method = "fixed", values = lam, automatic = FALSE, lambda0 = lam)
      ff <- tt_als_fit_sequential(
        y, basis, ranks, fs, control, penalty_order,
        init_cores = fit_final$cores, offset = offset, weights = w,
        linear = linear, smooth = smooth
      )
      list(criterion = ff$deviance, fit = ff)
    })
    lambda0_table <- sel$table
    fs <- list(method = "fixed", values = sel$lambda, automatic = FALSE,
               lambda0 = sel$lambda)
    fit_final <- tt_als_fit_sequential(
      y, basis, ranks, fs, control, penalty_order,
      init_cores = fit_final$cores, offset = offset, weights = w,
      linear = linear, smooth = smooth
    )
    lambda <- fit_final$lambda
  }

  prop_all <- if (length(proposal_hist)) do.call(rbind, proposal_hist) else NULL
  list(
    cores = fit_final$cores,
    intercept = fit_final$intercept,
    beta = fit_final$beta,
    linear = linear,
    smooth = fit_final$smooth,
    lambda = fit_final$lambda,
    ranks = ranks,
    eta = fit_final$eta,
    mu = fit_final$mu,
    deviance = fit_final$deviance,
    n_sweeps = fit_final$n_sweeps,
    n_outer = n_outer,
    n_pirls = NA_integer_,
    n_criterion_evals = n_eval + fit_final$n_criterion_evals,
    history = history,
    penalties = fit_final$penalties,
    penalty_mode = "global",
    cgcv = list(
      update = "outer_simultaneous",
      parameterization = param,
      damping = rho,
      max_log10_step = delta,
      proposals = prop_all,
      lambda0_table = lambda0_table,
      trace = prop_all
    ),
    q_descent = list(checked = FALSE),
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = TRUE,
    method_lambda = "cGCV",
    optimizer = "ALS"
  )
}

#' Rcpp entry for Gaussian ALS — redirects to the global-penalty R path.
#'
#' Legacy full-sweep Rcpp ALS used the own-margin surrogate and is no longer
#' a production fitter. Penalty helpers remain available via Rcpp.
#' @keywords internal
tt_als_fit_rcpp <- function(y, basis, ranks, lambda_spec, control, penalty_order = 2,
                            init_cores = NULL, offset = NULL, weights = NULL,
                            linear = NULL, smooth = NULL) {
  out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                    init_cores = init_cores, offset = offset, weights = weights,
                    linear = linear, smooth = smooth)
  out$backend <- "R"
  out
}
