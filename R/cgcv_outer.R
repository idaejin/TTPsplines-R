# Outer / simultaneous / damped cGCV under the classical global penalty.
#
# Statistical model unchanged: Q(Θ; λ) = loss + J_λ^global(Θ).
# This file only changes how λ is updated (Jacobi / damping / trust /
# scale–anisotropy), not the form of J_λ.

#' Spectral (operator) norm of a symmetric matrix (dense, small m).
#' @keywords internal
#' @noRd
.matrix_op_norm <- function(A) {
  A <- as.matrix(A)
  if (!nrow(A)) return(0)
  A <- (A + t(A)) / 2
  ev <- tryCatch(
    eigen(A, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  if (anyNA(ev)) {
    return(as.numeric(norm(A, type = "F")))
  }
  max(abs(as.numeric(ev)))
}

#' Normalize margin update order (permutation of 1:d).
#' @keywords internal
#' @noRd
.cgcv_margin_order <- function(order, d) {
  d <- as.integer(d)
  if (is.null(order)) return(seq_len(d))
  order <- as.integer(order)
  if (length(order) != d || anyDuplicated(order) ||
      !setequal(order, seq_len(d))) {
    stop("`cgcv_margin_order` must be a permutation of 1:d.", call. = FALSE)
  }
  order
}

#' Resolve cGCV update mode.
#' @keywords internal
#' @noRd
.cgcv_update_mode <- function(control) {
  mode <- control$cgcv_update %||% "sequential"
  mode <- as.character(mode)[[1L]]
  match.arg(mode, c("sequential", "outer_simultaneous"))
}

#' Resolve cGCV parameterization.
#' @keywords internal
#' @noRd
.cgcv_parameterization <- function(control) {
  p <- control$cgcv_parameterization %||% "free"
  p <- as.character(p)[[1L]]
  match.arg(p, c("free", "scale_anisotropy"))
}

#' Damping + trust-region clip on log λ.
#'
#' \deqn{\log\lambda^{damp}=(1-\rho)\log\lambda^{\mathrm{old}}+\rho\log\tilde\lambda}
#' then clip so \eqn{|\Delta\log_{10}\lambda|\le\Delta_{\max}}, then project
#' onto \code{bounds}.
#'
#' @keywords internal
#' @noRd
.cgcv_damped_trust_update <- function(lambda_old,
                                      lambda_tilde,
                                      rho = 1,
                                      max_log10_step = Inf,
                                      bounds = c(1e-4, 1e4)) {
  lambda_old <- as.numeric(lambda_old)
  lambda_tilde <- as.numeric(lambda_tilde)
  stopifnot(length(lambda_old) == length(lambda_tilde))
  rho <- as.numeric(rho)[[1L]]
  if (!is.finite(rho) || rho <= 0 || rho > 1) {
    stop("`cgcv_damping` (rho) must lie in (0, 1].", call. = FALSE)
  }
  max_log10_step <- as.numeric(max_log10_step)[[1L]]
  if (!is.finite(max_log10_step) || max_log10_step < 0) {
    max_log10_step <- Inf
  }
  bounds <- as.numeric(bounds)
  lo <- log(bounds[1])
  hi <- log(bounds[2])

  log_old <- log(pmax(lambda_old, .Machine$double.xmin))
  log_til <- log(pmax(lambda_tilde, .Machine$double.xmin))
  log_damp <- (1 - rho) * log_old + rho * log_til

  # Trust region on log10
  if (is.finite(max_log10_step)) {
    max_ln <- max_log10_step * log(10)
    delta <- log_damp - log_old
    delta <- pmin(pmax(delta, -max_ln), max_ln)
    log_clip <- log_old + delta
  } else {
    log_clip <- log_damp
  }

  log_final <- pmin(pmax(log_clip, lo), hi)
  lambda_new <- exp(log_final)
  boundary <- .lambda_boundary_status(lambda_new, bounds)

  data.frame(
    margin = seq_along(lambda_old),
    lambda_old = lambda_old,
    lambda_tilde = lambda_tilde,
    lambda_damped = exp(pmin(pmax(log_damp, lo), hi)),
    lambda_new = lambda_new,
    log10_old = log_old / log(10),
    log10_tilde = log_til / log(10),
    log10_new = log_final / log(10),
    rho = rho,
    max_log10_step = ifelse(is.finite(max_log10_step), max_log10_step, NA_real_),
    boundary = boundary,
    stringsAsFactors = FALSE
  )
}

#' Geometric-mean / anisotropy factorization of a λ vector.
#' @keywords internal
#' @noRd
.cgcv_scale_anisotropy_from_lambda <- function(lambda) {
  lambda <- as.numeric(lambda)
  stopifnot(all(lambda > 0), all(is.finite(lambda)))
  log_lam <- log(lambda)
  log_g <- mean(log_lam)
  list(
    lambda0 = exp(log_g),
    omega = exp(log_lam - log_g),
    log_omega = log_lam - log_g
  )
}

#' Recentering log-ω so mean is zero (prod ω = 1).
#' @keywords internal
#' @noRd
.cgcv_recenter_log_omega <- function(log_omega) {
  log_omega - mean(log_omega)
}

#' Apply free→scale-anisotropy map for one outer proposal.
#'
#' Uses free tilde λ only for relative anisotropy; overall scale λ0 is held
#' unless `lambda0_new` is supplied.
#' @keywords internal
#' @noRd
.cgcv_update_scale_anisotropy <- function(lambda_old,
                                          lambda_tilde,
                                          rho = 0.25,
                                          max_log10_step = 1,
                                          bounds = c(1e-4, 1e4),
                                          lambda0_new = NULL) {
  old <- .cgcv_scale_anisotropy_from_lambda(lambda_old)
  til <- .cgcv_scale_anisotropy_from_lambda(lambda_tilde)
  log_om_new <- (1 - rho) * old$log_omega + rho * til$log_omega
  log_om_new <- .cgcv_recenter_log_omega(log_om_new)
  omega_new <- exp(log_om_new)
  lam0 <- if (is.null(lambda0_new)) old$lambda0 else as.numeric(lambda0_new)[[1L]]
  # Trust-region / bounds applied to the reconstructed free λ
  prop <- .cgcv_damped_trust_update(
    lambda_old = lambda_old,
    lambda_tilde = lam0 * omega_new,
    rho = 1, # already damped in ω-space
    max_log10_step = max_log10_step,
    bounds = bounds
  )
  # Rebuild ω from clipped λ so prod ω = 1 remains exact
  final <- .cgcv_scale_anisotropy_from_lambda(prop$lambda_new)
  prop$lambda0_old <- old$lambda0
  prop$lambda0_new <- final$lambda0
  prop$omega_old <- old$omega
  prop$omega_new <- final$omega
  prop$parameterization <- "scale_anisotropy"
  prop
}

#' Build conditional cGCV / core-update workspace for margin k (global P_full).
#'
#' @param use_spectral If `NULL`, uses `control$use_spectral_gcv`. Pass
#'   `FALSE` for fixed-λ solves (spectral factorization is only useful for
#'   repeated cGCV evaluations).
#' @param compute_op_norms Operator-norm diagnostics for cGCV traces (expensive
#'   eigen); leave `FALSE` on the hot ALS/PIRLS path.
#' @keywords internal
#' @noRd
.cgcv_core_workspace <- function(cores, k, lambda, basis, target, ranks,
                                 control, weight = NULL, penalty_order = 2L,
                                 use_spectral = NULL,
                                 compute_op_norms = FALSE) {
  L <- left_interfaces(cores, basis)
  R <- right_interfaces(cores, basis)
  Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
  pen_k <- tt_conditional_penalty_full(
    cores, k, lambda,
    penalty_order = penalty_order,
    cyclic = attr(basis, "cyclic")
  )
  use_spec <- if (is.null(use_spectral)) {
    isTRUE(control$use_spectral_gcv)
  } else {
    isTRUE(use_spectral)
  }
  ws <- make_core_workspace(
    target, Xk, pen_k$P_own, lambda[k],
    control$lambda_bounds, control$tol_lambda,
    weight = weight,
    use_spectral = use_spec,
    P0 = pen_k$P_other
  )
  list(
    workspace = ws,
    P_own = pen_k$P_own,
    P_other = pen_k$P_other,
    P_own_op = if (isTRUE(compute_op_norms)) .matrix_op_norm(pen_k$P_own) else NA_real_,
    P_other_op = if (isTRUE(compute_op_norms)) .matrix_op_norm(pen_k$P_other) else NA_real_
  )
}

#' Evaluate cGCV / ed / RSS at a given λ_k without writing cores.
#' Uses `ws$spectral` built once in [make_core_workspace()].
#' @keywords internal
#' @noRd
.cgcv_eval_at <- function(ws, lambda_k) {
  cache <- .cgcv_spectral_from_workspace(ws)
  if (is.null(cache)) {
    .conditional_gcv(ws$yw, ws$Xw, ws$S, ws$P, ws$b, lambda_k, P0 = ws$P0)
  } else {
    .conditional_gcv_spectral(
      ws$yw, ws$Xw, ws$S, ws$P, ws$b, lambda_k, cache, P0 = ws$P0
    )
  }
}

#' Jacobi proposals: for every margin, tilde λ_k from the same frozen fit.
#'
#' Does **not** write λ or cores. Optionally returns dense grid curves.
#' @keywords internal
#' @noRd
.cgcv_propose_all <- function(cores, lambda, basis, target, ranks, control,
                              weight = NULL, penalty_order = 2L,
                              grid = NULL, eval_old = TRUE) {
  d <- length(cores)
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  rows <- vector("list", d)
  curves <- if (is.null(grid)) NULL else vector("list", d)
  n_eval <- 0L

  for (k in seq_len(d)) {
    built <- .cgcv_core_workspace(
      cores, k, lambda, basis, target, ranks, control,
      weight = weight, penalty_order = penalty_order,
      use_spectral = isTRUE(control$use_spectral_gcv),
      compute_op_norms = TRUE
    )
    ws <- built$workspace
    upd <- update_lambda_cgcv(ws)
    n_eval <- n_eval + upd$n_eval
    old_fit <- if (isTRUE(eval_old)) .cgcv_eval_at(ws, lambda[k]) else NULL
    status <- .lambda_boundary_status(upd$lambda, bounds)
    rows[[k]] <- data.frame(
      margin = k,
      lambda_old = lambda[k],
      lambda_tilde = upd$lambda,
      log10_old = log10(lambda[k]),
      log10_tilde = log10(upd$lambda),
      ed_tilde = upd$ed,
      gcv_tilde = upd$value,
      gcv_old = if (is.null(old_fit)) NA_real_ else old_fit$value,
      ed_old = if (is.null(old_fit)) NA_real_ else old_fit$ed,
      rss_tilde = NA_real_,
      P_other_op = built$P_other_op,
      P_own_op = built$P_own_op,
      lambda_P_own_op = upd$lambda * built$P_own_op,
      boundary_tilde = status,
      stringsAsFactors = FALSE
    )
    # attach rss if available from last fit object
    fit_til <- .cgcv_eval_at(ws, upd$lambda)
    rows[[k]]$rss_tilde <- fit_til$rss

    if (!is.null(grid)) {
      gvals <- vapply(grid, function(lam) .cgcv_eval_at(ws, lam)$value, numeric(1))
      curves[[k]] <- data.frame(
        margin = k, lambda = grid, log10_lambda = log10(grid),
        gcv = gvals, stringsAsFactors = FALSE
      )
    }
  }

  proposals <- do.call(rbind, rows)
  list(
    proposals = proposals,
    curves = if (is.null(curves)) NULL else do.call(rbind, curves),
    n_eval = n_eval
  )
}

#' Frozen-fit conditional cGCV curves (diagnostic; no λ updates).
#'
#' Fit (or reuse) a TT model at fixed \(\boldsymbol\lambda\), freeze cores,
#' then evaluate each margin's conditional GCV on a log-grid using the same
#' frozen \(P_{k,-k}\) / design.
#'
#' @param fit A `"ttpspline"` fit (preferably fixed-λ) or the list returned by
#'   an internal ALS/PIRLS fitter with `cores`, `lambda`, `intercept`.
#' @param y,X,basis Optional; taken from `fit` when available.
#' @param grid Numeric λ grid (default log-spaced over `fit$lambda_bounds`).
#' @param weight Optional PIRLS/Gaussian weights for the working response.
#' @param z Optional working response (GLM); default uses `y` / Gaussian.
#' @return List with `proposals`, `curves`, `lambda_frozen`, `penalty_mode`.
#' @export
tt_cgcv_frozen_curves <- function(fit,
                                  y = NULL,
                                  X = NULL,
                                  basis = NULL,
                                  grid = NULL,
                                  weight = NULL,
                                  z = NULL,
                                  offset = NULL,
                                  penalty_order = NULL) {
  stopifnot(!is.null(fit$cores), !is.null(fit$lambda))
  cores <- fit$cores
  lambda <- as.numeric(fit$lambda)
  d <- length(cores)
  ranks <- tt_ranks_from_cores(cores)
  if (is.null(basis)) {
    if (is.null(fit$knots) || is.null(fit$X)) {
      stop("Pass basis= or a full ttpspline fit with knots/X.", call. = FALSE)
    }
    basis <- eval_marginal_bases(fit$X, fit$knots, fit$degree,
                                 cyclic = fit$cyclic)
  }
  if (is.null(y)) y <- fit$y
  if (is.null(y)) stop("`y` required.", call. = FALSE)
  intercept <- fit$intercept %||% 0
  offset <- normalize_offset(offset %||% fit$offset, length(y))
  po <- as.integer(penalty_order %||% fit$penalty_order %||% 2L)
  bounds <- fit$lambda_bounds %||% fit$control$lambda_bounds %||% c(1e-4, 1e4)
  if (is.null(grid)) {
    grid <- exp(seq(log(bounds[1]), log(bounds[2]), length.out = 61L))
  }
  control <- fit$control %||% tt_control(lambda_bounds = bounds)
  control$lambda_bounds <- bounds

  if (is.null(z)) {
    target <- y - offset - intercept
  } else {
    target <- as.numeric(z) - offset - intercept
  }

  prop <- .cgcv_propose_all(
    cores, lambda, basis, target, ranks, control,
    weight = weight, penalty_order = po, grid = grid, eval_old = TRUE
  )
  list(
    proposals = prop$proposals,
    curves = prop$curves,
    lambda_frozen = lambda,
    penalty_mode = fit$penalty_mode %||% "global",
    grid = grid,
    n_eval = prop$n_eval
  )
}

#' Ranks chain from a list of TT cores.
#' @keywords internal
#' @noRd
tt_ranks_from_cores <- function(cores) {
  d <- length(cores)
  ranks <- integer(d + 1L)
  ranks[1] <- dim(cores[[1]])[1]
  for (k in seq_len(d)) ranks[k + 1L] <- dim(cores[[k]])[3]
  ranks
}

#' One Jacobi proposal + damped/trust (or scale–anisotropy) λ update.
#' @keywords internal
#' @noRd
.cgcv_simultaneous_step <- function(cores, lambda, basis, target, ranks,
                                    control, weight = NULL,
                                    penalty_order = 2L) {
  rho <- control$cgcv_damping %||% 1
  delta <- control$cgcv_max_log10_step %||% Inf
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  param <- .cgcv_parameterization(control)

  prop <- .cgcv_propose_all(
    cores, lambda, basis, target, ranks, control,
    weight = weight, penalty_order = penalty_order, grid = NULL
  )
  tilde <- prop$proposals$lambda_tilde

  if (identical(param, "scale_anisotropy")) {
    upd <- .cgcv_update_scale_anisotropy(
      lambda_old = lambda,
      lambda_tilde = tilde,
      rho = rho,
      max_log10_step = delta,
      bounds = bounds,
      lambda0_new = NULL
    )
  } else {
    upd <- .cgcv_damped_trust_update(
      lambda_old = lambda,
      lambda_tilde = tilde,
      rho = rho,
      max_log10_step = delta,
      bounds = bounds
    )
    upd$parameterization <- "free"
  }

  # Merge proposal diagnostics
  upd$gcv_old <- prop$proposals$gcv_old
  upd$gcv_tilde <- prop$proposals$gcv_tilde
  upd$ed_old <- prop$proposals$ed_old
  upd$ed_tilde <- prop$proposals$ed_tilde
  upd$P_other_op <- prop$proposals$P_other_op
  upd$P_own_op <- prop$proposals$P_own_op
  upd$lambda_P_own_op <- prop$proposals$lambda_P_own_op

  list(
    lambda = as.numeric(upd$lambda_new),
    proposals = upd,
    n_eval = prop$n_eval
  )
}

#' Select overall scale λ0 on a log-grid with fixed anisotropy ω.
#'
#' For each candidate λ0, set λ = λ0·ω and call `fit_fixed(lambda)`.
#' `fit_fixed` must return a list with numeric `criterion` (lower better)
#' and optional `fit` payload.
#' @keywords internal
#' @noRd
.cgcv_select_lambda0_grid <- function(omega, grid, fit_fixed) {
  omega <- as.numeric(omega)
  omega <- omega / exp(mean(log(omega)))
  rows <- lapply(as.numeric(grid), function(l0) {
    lam <- l0 * omega
    out <- fit_fixed(lam)
    data.frame(
      lambda0 = l0,
      criterion = as.numeric(out$criterion),
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  best <- which.min(tab$criterion)
  list(
    lambda0 = tab$lambda0[best],
    omega = omega,
    lambda = tab$lambda0[best] * omega,
    table = tab
  )
}
