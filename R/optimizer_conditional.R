# Conditional GLM core solvers: Damped-Newton-ALS and Block-LBFGS-ALS.
# Same penalized likelihood as global LBFGS; alternating over TT cores.

# ---------------------------------------------------------------------------
# Shared conditional objective / gradient / Hessian
# ---------------------------------------------------------------------------

#' Negative log-likelihood + optional scalar intercept (no TT penalty).
#' @keywords internal
#' @noRd
.tt_glm_nll <- function(y, eta, family) {
  key <- family_key(family)
  if (identical(key, "gaussian")) {
    return(0.5 * sum((y - eta)^2))
  }
  if (identical(key, "bernoulli")) {
    mu <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12)
    return(-sum(y * log(mu) + (1 - y) * log(1 - mu)))
  }
  if (identical(key, "poisson")) {
    mu <- pmax(exp(pmin(pmax(eta, -20), 20)), 1e-12)
    return(sum(mu - y * log(mu)))
  }
  stop("Unsupported family in .tt_glm_nll", call. = FALSE)
}

#' Mean and diagonal GLM weight for Newton Hessian (canonical links).
#'
#' Uses the *true* variance function (no PIRLS weight floor). Tiny clipping
#' of mu only for numerical safety of Bernoulli/Poisson.
#' @keywords internal
#' @noRd
.tt_glm_mu_weight <- function(y, eta, family) {
  key <- family_key(family)
  if (identical(key, "gaussian")) {
    return(list(mu = eta, weight = rep(1, length(y)), score = y - eta))
  }
  if (identical(key, "bernoulli")) {
    mu <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12)
    w <- mu * (1 - mu)
    return(list(mu = mu, weight = w, score = y - mu))
  }
  if (identical(key, "poisson")) {
    mu <- pmax(exp(pmin(pmax(eta, -20), 20)), 1e-12)
    return(list(mu = mu, weight = mu, score = y - mu))
  }
  stop("Unsupported family", call. = FALSE)
}

#' Conditional core objective Q_k(g) = nll(intercept + Xk g) + (λ/2) g'P g.
#' @keywords internal
#' @noRd
.tt_conditional_qk <- function(g, y, intercept, Xk, Pk, lambda_k, family) {
  eta <- intercept + as.numeric(Xk %*% g)
  nll <- .tt_glm_nll(y, eta, family)
  pen <- 0.5 * lambda_k * as.numeric(crossprod(g, Pk %*% g))
  list(value = nll + pen, nll = nll, penalty = pen, eta = eta)
}

#' Conditional gradient ∇Q_k = -X'(y-μ) + λ P g  (equiv. X'(μ-y) + λPg).
#' @keywords internal
#' @noRd
.tt_conditional_qk_grad <- function(g, y, intercept, Xk, Pk, lambda_k, family) {
  eta <- intercept + as.numeric(Xk %*% g)
  mw <- .tt_glm_mu_weight(y, eta, family)
  # score = y - mu; NLL grad_eta = -(y-mu); grad_g = X'(-(y-mu)) = X'(mu-y)
  grad_nll <- as.numeric(crossprod(Xk, mw$mu - y))
  grad_pen <- as.numeric(lambda_k * (Pk %*% g))
  list(
    grad = grad_nll + grad_pen,
    eta = eta,
    mu = mw$mu,
    weight = mw$weight,
    min_weight = min(mw$weight)
  )
}

#' Conditional Hessian H_k = X' W X + λ P.
#' @keywords internal
#' @noRd
.tt_conditional_qk_hess <- function(Xk, weight, Pk, lambda_k, ridge = 0) {
  # crossprod(X, w * X) = X' diag(w) X
  XtWX <- crossprod(Xk, Xk * weight)
  H <- XtWX + lambda_k * Pk
  if (ridge > 0) {
    H <- H + ridge * diag(nrow(H))
  }
  H
}

#' One damped Newton update of a single TT core (true-objective Armijo).
#' @keywords internal
#' @noRd
.tt_damped_newton_core <- function(g, y, intercept, Xk, Pk, lambda_k, family,
                                   control) {
  c_arm <- control$dn_armijo_c %||% 1e-4
  rho <- control$dn_step_factor %||% control$step_factor %||% 0.5
  step_min <- control$dn_step_min %||% control$step_min %||% 1e-12
  ridge <- control$dn_ridge %||% 0
  # report ridge if used (non-silent via return)
  q0 <- .tt_conditional_qk(g, y, intercept, Xk, Pk, lambda_k, family)
  gi <- .tt_conditional_qk_grad(g, y, intercept, Xk, Pk, lambda_k, family)
  if (!is.finite(q0$value) || any(!is.finite(gi$grad))) {
    return(list(
      g = g, accepted = FALSE, alpha = 0, n_backtrack = 0L,
      q_before = q0$value, q_after = q0$value,
      grad_norm = NA_real_, min_weight = gi$min_weight,
      max_abs_eta = max(abs(gi$eta)),
      reason = "non-finite objective/gradient", ridge = ridge
    ))
  }
  H <- .tt_conditional_qk_hess(Xk, gi$weight, Pk, lambda_k, ridge = ridge)
  delta <- tryCatch(
    as.numeric(solve_spd_ridge(H, -gi$grad)),
    error = function(e) NULL
  )
  if (is.null(delta) || any(!is.finite(delta))) {
    return(list(
      g = g, accepted = FALSE, alpha = 0, n_backtrack = 0L,
      q_before = q0$value, q_after = q0$value,
      grad_norm = max(abs(gi$grad)), min_weight = gi$min_weight,
      max_abs_eta = max(abs(gi$eta)),
      reason = "Hessian / Newton solve failed", ridge = ridge
    ))
  }
  dir_deriv <- sum(gi$grad * delta)
  if (!(dir_deriv < 0)) {
    # not a descent direction (numerical); reject
    return(list(
      g = g, accepted = FALSE, alpha = 0, n_backtrack = 0L,
      q_before = q0$value, q_after = q0$value,
      grad_norm = max(abs(gi$grad)), min_weight = gi$min_weight,
      max_abs_eta = max(abs(gi$eta)),
      reason = "not a descent direction", ridge = ridge
    ))
  }

  alpha <- 1
  n_bt <- 0L
  accepted <- FALSE
  g_new <- g
  q_new <- q0$value
  while (alpha + 1e-18 >= step_min) {
    g_try <- g + alpha * delta
    q_try <- .tt_conditional_qk(g_try, y, intercept, Xk, Pk, lambda_k, family)
    if (is.finite(q_try$value) &&
        q_try$value <= q0$value + c_arm * alpha * dir_deriv) {
      g_new <- g_try
      q_new <- q_try$value
      accepted <- TRUE
      break
    }
    alpha <- alpha * rho
    n_bt <- n_bt + 1L
  }
  list(
    g = if (accepted) g_new else g,
    accepted = accepted,
    alpha = if (accepted) alpha else 0,
    n_backtrack = n_bt,
    q_before = q0$value,
    q_after = if (accepted) q_new else q0$value,
    grad_norm = max(abs(gi$grad)),
    min_weight = gi$min_weight,
    max_abs_eta = max(abs(gi$eta)),
    reason = if (accepted) "ok" else "line search failed",
    ridge = ridge
  )
}

#' L-BFGS minimization of Q_k (analytical grad).
#' @keywords internal
#' @noRd
.tt_block_lbfgs_core <- function(g, y, intercept, Xk, Pk, lambda_k, family,
                                control) {
  maxit <- as.integer(control$block_lbfgs_maxit %||% 50L)
  q0 <- .tt_conditional_qk(g, y, intercept, Xk, Pk, lambda_k, family)
  fn <- function(th) {
    .tt_conditional_qk(th, y, intercept, Xk, Pk, lambda_k, family)$value
  }
  gr <- function(th) {
    .tt_conditional_qk_grad(th, y, intercept, Xk, Pk, lambda_k, family)$grad
  }
  opt <- stats::optim(
    g, fn, gr,
    method = "L-BFGS-B",
    control = list(
      maxit = maxit,
      factr = 1e7,
      pgtol = 1e-8,
      trace = 0L
    )
  )
  q1 <- .tt_conditional_qk(opt$par, y, intercept, Xk, Pk, lambda_k, family)
  list(
    g = as.numeric(opt$par),
    q_before = q0$value,
    q_after = q1$value,
    n_iter = opt$counts[["function"]],
    convergence = opt$convergence,
    accepted = is.finite(q1$value) && q1$value <= q0$value + 1e-8,
    max_abs_eta = max(abs(q1$eta)),
    reason = if (opt$convergence %in% c(0L, 1L)) "ok" else "L-BFGS nonconvergence"
  )
}

# ---------------------------------------------------------------------------
# Fit wrappers
# ---------------------------------------------------------------------------

.tt_conditional_init <- function(y, basis, ranks, lambda_spec, control,
                                 penalty_order, init_cores, family) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  method <- lambda_spec$method
  if (!identical(method, "fixed")) {
    stop("Damped-Newton-ALS / LBFGS-ALS currently require fixed lambda ",
         "(no cGCV in this gate).", call. = FALSE)
  }
  lambda <- as.numeric(lambda_spec$values %||% lambda_spec$lambda0)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p, ranks[k + 1L], penalty_order)
  })
  fam <- normalize_family(family)
  intercept <- init_intercept(fam, y)
  list(
    d = d, p = p, ranks = ranks, cores = cores, penalties = penalties,
    lambda = lambda, intercept = intercept, fam = fam, family = fam
  )
}

#' Damped-Newton alternating core updates (true-objective Armijo).
#'
#' Distinct from PIRLS-ALS: each Newton step is accepted only if the true
#' conditional penalized likelihood decreases (Armijo). Intercept fixed at
#' GLM init (same convention as global LBFGS).
#'
#' @keywords internal
#' @noRd
tt_damped_newton_als_fit <- function(y, basis, family, ranks, lambda_spec,
                                     control, penalty_order = 2,
                                     init_cores = NULL) {
  t0 <- proc.time()[["elapsed"]]
  st <- .tt_conditional_init(
    y, basis, ranks, lambda_spec, control, penalty_order, init_cores, family
  )
  cores <- st$cores
  intercept <- st$intercept
  penalties <- st$penalties
  lambda <- st$lambda
  fam <- st$fam
  d <- st$d
  ranks <- st$ranks
  p <- st$p

  max_sweeps <- as.integer(control$dn_max_sweeps %||% control$max_sweeps %||% 40L)
  tol <- as.numeric(control$tol %||% 1e-8)
  n_backtrack_total <- 0L
  n_core_updates <- 0L
  alphas <- numeric()
  hist <- list()
  reason <- "maxit"
  converged <- FALSE

  obj <- tt_glm_penalized_objective(y, cores, intercept, basis, penalties, lambda, fam)

  for (sw in seq_len(max_sweeps)) {
    obj_old <- obj
    for (k in seq_len(d)) {
      L <- left_interfaces(cores, basis)
      R <- right_interfaces(cores, basis)
      Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
      g <- as.numeric(cores[[k]])
      upd <- .tt_damped_newton_core(
        g, y, intercept, Xk, penalties[[k]], lambda[k], fam, control
      )
      n_backtrack_total <- n_backtrack_total + upd$n_backtrack
      if (isTRUE(upd$accepted)) {
        cores[[k]] <- array(upd$g, c(ranks[k], p, ranks[k + 1L]))
        n_core_updates <- n_core_updates + 1L
        alphas <- c(alphas, upd$alpha)
      }
    }
    obj <- tt_glm_penalized_objective(y, cores, intercept, basis, penalties, lambda, fam)
    if (isTRUE(control$trace)) {
      cat(sprintf(
        "  DN-ALS sweep %2d | L=%.6g | max|eta|=%.3g\n",
        sw, obj$value, max(abs(obj$eta))
      ))
    }
    hist[[sw]] <- list(
      sweep = sw, objective = obj$value, max_abs_eta = max(abs(obj$eta))
    )
    if (!is.finite(obj$value)) {
      reason <- "non-finite objective"
      converged <- FALSE
      break
    }
    rel <- abs(obj_old$value - obj$value) / max(1, abs(obj_old$value))
    if (sw > 1L && rel < tol) {
      reason <- "relative objective change"
      converged <- TRUE
      break
    }
    if (sw == max_sweeps) reason <- "maxit"
  }

  eta <- obj$eta
  mu <- invlink_eta(fam, eta)
  list(
    cores = cores,
    intercept = intercept,
    lambda = lambda,
    ranks = ranks,
    eta = eta,
    mu = mu,
    deviance = glm_deviance(fam, y, mu),
    n_sweeps = length(hist),
    n_pirls = NA_integer_,
    n_criterion_evals = NA_integer_,
    n_opt_iter = as.integer(n_core_updates),
    n_outer = NA_integer_,
    history = hist,
    penalties = penalties,
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = isTRUE(converged) && is.finite(obj$value),
    convergence = list(
      overall = isTRUE(converged) && is.finite(obj$value),
      pirls = NA,
      als = TRUE,
      reason = reason,
      n_pirls = NA_integer_,
      n_als_sweeps = length(hist),
      n_step_halvings = as.integer(n_backtrack_total),
      n_opt_iter = as.integer(n_core_updates),
      median_alpha = if (length(alphas)) stats::median(alphas) else NA_real_,
      n_backtrack = as.integer(n_backtrack_total)
    ),
    method_lambda = "fixed",
    optimizer = "Damped-Newton-ALS",
    backend = "R",
    objective = obj$value
  )
}

#' Block-wise L-BFGS alternating over TT cores (direct conditional likelihood).
#'
#' Distinct from global `LBFGS` (joint over all cores). Intercept fixed at
#' GLM init (same as global LBFGS).
#'
#' @keywords internal
#' @noRd
tt_lbfgs_als_fit <- function(y, basis, family, ranks, lambda_spec, control,
                             penalty_order = 2, init_cores = NULL) {
  t0 <- proc.time()[["elapsed"]]
  st <- .tt_conditional_init(
    y, basis, ranks, lambda_spec, control, penalty_order, init_cores, family
  )
  cores <- st$cores
  intercept <- st$intercept
  penalties <- st$penalties
  lambda <- st$lambda
  fam <- st$fam
  d <- st$d
  ranks <- st$ranks
  p <- st$p

  max_sweeps <- as.integer(control$block_lbfgs_sweeps %||% control$max_sweeps %||% 40L)
  tol <- as.numeric(control$tol %||% 1e-8)
  n_opt_iter <- 0L
  n_accepted <- 0L
  hist <- list()
  reason <- "maxit"
  converged <- FALSE

  obj <- tt_glm_penalized_objective(y, cores, intercept, basis, penalties, lambda, fam)

  for (sw in seq_len(max_sweeps)) {
    obj_old <- obj
    for (k in seq_len(d)) {
      L <- left_interfaces(cores, basis)
      R <- right_interfaces(cores, basis)
      Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
      g <- as.numeric(cores[[k]])
      upd <- .tt_block_lbfgs_core(
        g, y, intercept, Xk, penalties[[k]], lambda[k], fam, control
      )
      n_opt_iter <- n_opt_iter + (upd$n_iter %||% 0L)
      if (isTRUE(upd$accepted) ||
          (is.finite(upd$q_after) && upd$q_after <= upd$q_before + 1e-8)) {
        cores[[k]] <- array(upd$g, c(ranks[k], p, ranks[k + 1L]))
        n_accepted <- n_accepted + 1L
      }
    }
    obj <- tt_glm_penalized_objective(y, cores, intercept, basis, penalties, lambda, fam)
    if (isTRUE(control$trace)) {
      cat(sprintf(
        "  LBFGS-ALS sweep %2d | L=%.6g | max|eta|=%.3g\n",
        sw, obj$value, max(abs(obj$eta))
      ))
    }
    hist[[sw]] <- list(
      sweep = sw, objective = obj$value, max_abs_eta = max(abs(obj$eta))
    )
    if (!is.finite(obj$value)) {
      reason <- "non-finite objective"
      converged <- FALSE
      break
    }
    rel <- abs(obj_old$value - obj$value) / max(1, abs(obj_old$value))
    if (sw > 1L && rel < tol) {
      reason <- "relative objective change"
      converged <- TRUE
      break
    }
    if (sw == max_sweeps) reason <- "maxit"
  }

  eta <- obj$eta
  mu <- invlink_eta(fam, eta)
  list(
    cores = cores,
    intercept = intercept,
    lambda = lambda,
    ranks = ranks,
    eta = eta,
    mu = mu,
    deviance = glm_deviance(fam, y, mu),
    n_sweeps = length(hist),
    n_pirls = NA_integer_,
    n_criterion_evals = NA_integer_,
    n_opt_iter = as.integer(n_opt_iter),
    n_outer = NA_integer_,
    history = hist,
    penalties = penalties,
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = isTRUE(converged) && is.finite(obj$value),
    convergence = list(
      overall = isTRUE(converged) && is.finite(obj$value),
      pirls = NA,
      als = TRUE,
      reason = reason,
      n_pirls = NA_integer_,
      n_als_sweeps = length(hist),
      n_step_halvings = 0L,
      n_opt_iter = as.integer(n_opt_iter),
      n_accepted_cores = as.integer(n_accepted)
    ),
    method_lambda = "fixed",
    optimizer = "LBFGS-ALS",
    backend = "R",
    objective = obj$value
  )
}
