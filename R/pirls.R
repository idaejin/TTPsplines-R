#' GLM PIRLS + weighted TT-ALS (fixed or cGCV λ), R backend.
#'
#' Bernoulli uses consistent W/z working variance and optional true-objective
#' outer step-halving (see [tt_control()]).
#'
#' @keywords internal
tt_pirls_fit <- function(y, basis, family, ranks, lambda_spec, control,
                         penalty_order = 2, init_cores = NULL, offset = NULL) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  method <- lambda_spec$method
  lambda <- lambda_spec$values %||% lambda_spec$lambda0
  fam <- normalize_family(family)
  key <- family_key(fam)
  offset <- normalize_offset(offset, length(y))
  intercept <- init_intercept(fam, y, offset = offset)
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
    for (k in seq_len(d)) cores[[k]] <- cores[[k]] * 0.05
  } else {
    cores <- init_cores
  }
  penalties <- lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p, ranks[k + 1L], penalty_order)
  })

  eta <- tt_eta(offset, intercept, cores, basis)
  mu <- invlink_eta(fam, eta)
  dev <- glm_deviance(fam, y, mu)
  obj <- tt_glm_penalized_objective(y, cores, intercept, basis, penalties, lambda, fam, offset = offset)

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
    w <- work$weight
    z <- work$z

    # --- inner weighted ALS ---
    for (sw in seq_len(control$als_sweeps_per_pirls)) {
      zc <- z - offset - intercept
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
      intercept <- sum(w * (z - offset - f)) / max(sum(w), 1e-12)
      n_als_sweeps_total <- n_als_sweeps_total + 1L
    }

    cores_cand <- cores
    intercept_cand <- intercept
    eta_cand <- tt_eta(offset, intercept_cand, cores_cand, basis)
    obj_cand <- tt_glm_penalized_objective(
      y, cores_cand, intercept_cand, basis, penalties, lambda, fam, offset = offset
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
          offset = offset
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
        dev <- glm_deviance(fam, y, mu)
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
          line_ok = FALSE
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
    dev <- glm_deviance(fam, y, mu)
    n_pirls <- it
    hist_rows[[length(hist_rows) + 1L]] <- .pirls_hist_row(
      it, obj, y, eta, fam, work, control$als_sweeps_per_pirls,
      proposed_step = 1, accepted_step = accepted_step, n_halve_it = n_halve_it,
      line_ok = TRUE
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
      rel_dev <- abs(glm_deviance(fam, y, invlink_eta(fam, eta_old)) - dev) /
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
    mu = mu,
    deviance = dev,
    n_sweeps = NA_integer_,
    n_pirls = n_pirls,
    n_criterion_evals = n_eval,
    history = history,
    penalties = penalties,
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = convergence$overall,
    convergence = convergence,
    method_lambda = method,
    optimizer = "ALS",
    backend = "R",
    objective = obj$value
  )
}

.pirls_hist_row <- function(it, obj, y, eta, family, work, als_sweeps,
                            proposed_step, accepted_step, n_halve_it, line_ok) {
  mu <- invlink_eta(family, eta)
  data.frame(
    pirls = it,
    objective = obj$value,
    nll = obj$nll,
    penalty = obj$penalty,
    deviance = glm_deviance(family, y, mu),
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
                              penalty_order = 2, init_cores = NULL, offset = NULL) {
  offset <- normalize_offset(offset, length(y))
  if (any(offset != 0)) {
    out <- tt_pirls_fit(y, basis, family, ranks, lambda_spec, control,
                        penalty_order, init_cores = init_cores, offset = offset)
    out$backend <- "R"
    return(out)
  }
  key <- family_key(family)
  do_halving <- identical(key, "bernoulli") &&
    isTRUE(control$pirls_step_halving %||% control$damping %||% TRUE)
  if (do_halving || !exists("tt_glm_pirls_cgcv_cpp", mode = "function")) {
    return(tt_pirls_fit(y, basis, family, ranks, lambda_spec, control,
                        penalty_order, init_cores = init_cores, offset = offset))
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
