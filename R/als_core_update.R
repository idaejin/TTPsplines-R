#' Single-core Gaussian ALS update under global \(P_k^{\mathrm{full}}\) (fixed λ).
#'
#' **P1 building block** for a future fixed-λ C++ fitter used by global-GCV.
#' Updates **one** margin `k` only: no cGCV, no sweep loop, Gaussian,
#' `null_space = "joint"` (no null projection). Numeric λ only.
#'
#' @param y Response (length `n`).
#' @param cores Length-`d` list of TT cores (not mutated; a clone is updated).
#' @param intercept Scalar intercept used to form the residual.
#' @param basis Length-`d` list of B-spline bases (`n × p`).
#' @param k Margin index in `1:d`.
#' @param lambda Numeric length-`d` (or scalar recycled) smoothing parameters.
#' @param offset,weights Optional offset / observation weights.
#' @param penalty_order Difference penalty order (default 2).
#' @param backend `"R"` (reference via [.cgcv_core_workspace] +
#'   [update_lambda_fixed]) or `"Rcpp"` ([tt_als_core_update_global_cpp]).
#' @param y_centered Optional precomputed residual
#'   `y - offset - intercept` (skips recomputation when supplied).
#' @return List with `S`, `b`, `P_own`, `P_other`, `P_full`, `g`, updated
#'   `cores`, global `rss` / `penalty` / `objective` (= \(L_{\mathrm{pen}}\)),
#'   and conditional diagnostics.
#' @keywords internal
#' @noRd
tt_als_core_update_global <- function(y, cores, intercept, basis, k, lambda,
                                      offset = NULL, weights = NULL,
                                      penalty_order = 2L,
                                      backend = c("R", "Rcpp"),
                                      y_centered = NULL) {
  backend <- match.arg(backend)
  d <- length(cores)
  k <- as.integer(k)
  if (k < 1L || k > d) stop("`k` out of range.", call. = FALSE)
  y <- as.numeric(y)
  n <- length(y)
  offset <- normalize_offset(offset, n)
  w <- normalize_weights(weights, n)
  intercept <- as.numeric(intercept)[1L]
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (length(lambda) != d) stop("`lambda` length must be 1 or d.", call. = FALSE)

  yc <- if (is.null(y_centered)) {
    as.numeric(y - offset - intercept)
  } else {
    as.numeric(y_centered)
  }
  if (length(yc) != n) stop("`y_centered` length mismatch.", call. = FALSE)

  # Deep-copy cores so callers' lists are never mutated
  cores0 <- lapply(cores, function(C) {
    array(as.numeric(C), dim = dim(C))
  })

  if (identical(backend, "Rcpp")) {
    if (!exists("tt_als_core_update_global_cpp", mode = "function")) {
      stop("tt_als_core_update_global_cpp not available.", call. = FALSE)
    }
    DtD_list <- tt_DtD_list(cores0, penalty_order, cyclic = attr(basis, "cyclic"))
    raw <- tt_als_core_update_global_cpp(
      yc = yc,
      cores_list = cores0,
      basis_list = basis,
      k = k,
      lambda = lambda,
      DtD_list = DtD_list,
      weight = w
    )
    g <- as.numeric(raw$g)
    S <- raw$S
    b <- as.numeric(raw$b)
    P_own <- raw$P_own
    P_other <- raw$P_other
    P_full <- raw$P_full
    rss_cond <- as.numeric(raw$rss_cond)
    penalty_cond <- as.numeric(raw$penalty_cond)
    q_cond <- as.numeric(raw$q_cond)
  } else {
    # Reference: same path as fixed-λ ALS core visit (dense Gram = legacy)
    ranks <- integer(d + 1L)
    ranks[1L] <- dim(cores0[[1L]])[1L]
    for (j in seq_len(d)) ranks[j + 1L] <- dim(cores0[[j]])[3L]
    ctrl <- tt_control(
      gram_method = "legacy",
      use_spectral_gcv = FALSE,
      compute_edf = FALSE
    )
    built <- .cgcv_core_workspace(
      cores0, k, lambda, basis, yc, ranks, ctrl,
      weight = w, penalty_order = penalty_order,
      use_spectral = FALSE, compute_op_norms = FALSE
    )
    # Expose full P for the gate (workspace stores P_own + P0 separately)
    pen_k <- list(
      P_own = built$P_own,
      P_other = built$P_other,
      P_full = as.numeric(lambda[k]) * built$P_own + built$P_other
    )
    pen_k$P_full <- 0.5 * (pen_k$P_full + t(pen_k$P_full))
    upd <- update_lambda_fixed(built$workspace)
    g <- as.numeric(upd$g)
    S <- built$workspace$S
    b <- as.numeric(built$workspace$b)
    P_own <- pen_k$P_own
    P_other <- pen_k$P_other
    P_full <- pen_k$P_full
    yw <- built$workspace$yw
    rss_cond <- as.numeric(
      sum(yw^2) - 2 * crossprod(b, g) + crossprod(g, S %*% g)
    )
    rss_cond <- max(rss_cond, 0)
    penalty_cond <- as.numeric(0.5 * crossprod(g, P_full %*% g))
    q_cond <- 0.5 * rss_cond + penalty_cond
  }

  cores_new <- cores0
  dm <- dim(cores_new[[k]])
  cores_new[[k]] <- array(g, dm)

  q <- tt_gaussian_Q(
    y, cores_new, intercept, basis, lambda,
    offset = offset, weights = w,
    penalty_order = penalty_order,
    cyclic = attr(basis, "cyclic")
  )

  list(
    S = S,
    b = b,
    P_own = P_own,
    P_other = P_other,
    P_full = P_full,
    g = g,
    cores = cores_new,
    rss = q$rss,
    penalty = q$penalty,
    objective = q$value,
    rss_cond = rss_cond,
    penalty_cond = penalty_cond,
    q_cond = q_cond,
    k = k,
    lambda = lambda,
    backend = backend
  )
}

#' One Gauss–Seidel ALS sweep under global \(P_k^{\mathrm{full}}\) (fixed λ).
#'
#' **P2 building block.** Visits every margin in `margin_order` (default `1:d`).
#' No cGCV, no multi-sweep loop, no intercept refresh inside the sweep
#' (matches one pass of fixed-λ ALS with frozen `yc`). Gaussian /
#' `null_space = "joint"` only.
#'
#' @param margin_order Integer permutation of `1:d` (LTR default; use `d:1`
#'   for RTL). Arbitrary permutations rebuild interfaces each core.
#' @inheritParams tt_als_core_update_global
#' @return List with updated `cores`, `f` (TT fit), global `rss` / `penalty` /
#'   `objective`, and `margin_order`.
#' @keywords internal
#' @noRd
tt_als_sweep_global <- function(y, cores, intercept, basis, lambda,
                                offset = NULL, weights = NULL,
                                penalty_order = 2L,
                                margin_order = NULL,
                                backend = c("R", "Rcpp"),
                                y_centered = NULL) {
  backend <- match.arg(backend)
  d <- length(cores)
  y <- as.numeric(y)
  n <- length(y)
  offset <- normalize_offset(offset, n)
  w <- normalize_weights(weights, n)
  intercept <- as.numeric(intercept)[1L]
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (length(lambda) != d) stop("`lambda` length must be 1 or d.", call. = FALSE)
  margin_order <- .cgcv_margin_order(margin_order, d)

  yc <- if (is.null(y_centered)) {
    as.numeric(y - offset - intercept)
  } else {
    as.numeric(y_centered)
  }
  if (length(yc) != n) stop("`y_centered` length mismatch.", call. = FALSE)

  cores0 <- lapply(cores, function(C) array(as.numeric(C), dim = dim(C)))

  if (identical(backend, "Rcpp")) {
    if (!exists("tt_als_sweep_global_cpp", mode = "function")) {
      stop("tt_als_sweep_global_cpp not available.", call. = FALSE)
    }
    DtD_list <- tt_DtD_list(cores0, penalty_order, cyclic = attr(basis, "cyclic"))
    raw <- tt_als_sweep_global_cpp(
      yc = yc,
      cores_list = cores0,
      basis_list = basis,
      lambda = lambda,
      DtD_list = DtD_list,
      weight = w,
      margin_order = as.integer(margin_order)
    )
    cores_new <- lapply(seq_len(d), function(j) {
      array(as.numeric(raw$cores[[j]]), dim = dim(cores0[[j]]))
    })
    f <- as.numeric(raw$f)
    order_out <- as.integer(raw$margin_order)
  } else {
    cores_new <- cores0
    for (k in margin_order) {
      step <- tt_als_core_update_global(
        y, cores_new, intercept, basis, k, lambda,
        offset = offset, weights = w, penalty_order = penalty_order,
        backend = "R", y_centered = yc
      )
      cores_new <- step$cores
    }
    f <- as.numeric(tt_contraction(cores_new, basis))
    order_out <- as.integer(margin_order)
  }

  q <- tt_gaussian_Q(
    y, cores_new, intercept, basis, lambda,
    offset = offset, weights = w,
    penalty_order = penalty_order,
    cyclic = attr(basis, "cyclic")
  )

  list(
    cores = cores_new,
    f = f,
    rss = q$rss,
    penalty = q$penalty,
    objective = q$value,
    eta = offset + intercept + f,
    lambda = lambda,
    margin_order = order_out,
    backend = backend
  )
}

#' Fixed-λ Gaussian ALS fitter under global \(P_k^{\mathrm{full}}\) (P3).
#'
#' Multi-sweep Gauss–Seidel with intercept refresh after each sweep.
#' No cGCV, no `linear=` / `smooth=`, `null_space = "joint"` only.
#' Stopping: relative RSS change `< tol` after sweep `> 2` (same as R ALS).
#'
#' @param max_sweeps Maximum sweeps.
#' @param tol Relative RSS tolerance (default `1e-8`).
#' @inheritParams tt_als_sweep_global
#' @return List with `cores`, `intercept`, `f`/`eta`, `rss`/`penalty`/`objective`,
#'   `n_sweeps`, `converged`, `history`, `backend`.
#' @keywords internal
#' @noRd
tt_als_fit_fixed_global <- function(y, cores, basis, lambda,
                                    offset = NULL, weights = NULL,
                                    penalty_order = 2L,
                                    margin_order = NULL,
                                    max_sweeps = 50L,
                                    tol = 1e-8,
                                    backend = c("R", "Rcpp")) {
  backend <- match.arg(backend)
  d <- length(cores)
  y <- as.numeric(y)
  n <- length(y)
  offset <- normalize_offset(offset, n)
  w <- normalize_weights(weights, n)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (length(lambda) != d) stop("`lambda` length must be 1 or d.", call. = FALSE)
  margin_order <- .cgcv_margin_order(margin_order, d)
  max_sweeps <- as.integer(max_sweeps)[1L]
  tol <- as.numeric(tol)[1L]
  if (max_sweeps < 1L) stop("`max_sweeps` must be >= 1.", call. = FALSE)

  cores0 <- lapply(cores, function(C) array(as.numeric(C), dim = dim(C)))

  if (identical(backend, "Rcpp")) {
    if (!exists("tt_als_fit_fixed_global_cpp", mode = "function")) {
      stop("tt_als_fit_fixed_global_cpp not available.", call. = FALSE)
    }
    DtD_list <- tt_DtD_list(cores0, penalty_order, cyclic = attr(basis, "cyclic"))
    raw <- tt_als_fit_fixed_global_cpp(
      y = y,
      cores_list = cores0,
      basis_list = basis,
      lambda = lambda,
      DtD_list = DtD_list,
      weight = w,
      offset = offset,
      max_sweeps = max_sweeps,
      tol = tol,
      margin_order = as.integer(margin_order)
    )
    cores_new <- lapply(seq_len(d), function(j) {
      array(as.numeric(raw$cores[[j]]), dim = dim(cores0[[j]]))
    })
    hist <- as.data.frame(raw$history)
    return(list(
      cores = cores_new,
      intercept = as.numeric(raw$intercept)[1L],
      f = as.numeric(raw$f),
      eta = as.numeric(raw$eta),
      mu = as.numeric(raw$mu),
      rss = as.numeric(raw$rss),
      penalty = as.numeric(raw$penalty),
      objective = as.numeric(raw$objective),
      lambda = as.numeric(raw$lambda),
      n_sweeps = as.integer(raw$n_sweeps),
      converged = isTRUE(raw$converged),
      convergence_reason = as.character(raw$convergence_reason),
      margin_order = as.integer(raw$margin_order),
      history = hist,
      backend = "Rcpp",
      method = "global_fixed_fit"
    ))
  }

  # R reference: same outer loop as C++ (P2 sweeps + intercept refresh)
  intercept <- tt_update_intercept_beta(y, offset, f = 0, weights = w)$intercept
  history <- vector("list", max_sweeps)
  prev_rss <- Inf
  prev_eta <- NULL
  n_sweeps <- 0L
  converged <- FALSE
  reason <- "max_sweeps"
  cores_cur <- cores0
  f <- eta <- NULL
  rss <- pen <- obj <- NA_real_

  for (sw in seq_len(max_sweeps)) {
    swp <- tt_als_sweep_global(
      y, cores_cur, intercept, basis, lambda,
      offset = offset, weights = w, penalty_order = penalty_order,
      margin_order = margin_order, backend = "R"
    )
    cores_cur <- swp$cores
    f <- swp$f
    intercept <- tt_update_intercept_beta(
      y, offset, f = f, weights = w
    )$intercept
    eta <- offset + intercept + f
    rss <- sum(w * (y - eta)^2)
    pen <- tt_global_penalty_value(
      cores_cur, lambda, penalty_order = penalty_order,
      cyclic = attr(basis, "cyclic")
    )
    obj <- 0.5 * rss + pen
    d_eta <- if (is.null(prev_eta)) NA_real_ else sqrt(mean((eta - prev_eta)^2))
    history[[sw]] <- data.frame(
      sweep = sw, rss = rss, objective = obj, penalty = pen, d_eta = d_eta,
      stringsAsFactors = FALSE
    )
    n_sweeps <- sw
    prev_eta <- eta
    if (sw > 2L) {
      if (abs(prev_rss - rss) / max(1, abs(prev_rss)) < tol) {
        converged <- TRUE
        reason <- "tol_rss"
        break
      }
    }
    prev_rss <- rss
  }

  list(
    cores = cores_cur,
    intercept = intercept,
    f = as.numeric(f),
    eta = as.numeric(eta),
    mu = as.numeric(eta),
    rss = rss,
    penalty = pen,
    objective = obj,
    lambda = lambda,
    n_sweeps = n_sweeps,
    converged = converged,
    convergence_reason = reason,
    margin_order = as.integer(margin_order),
    history = do.call(rbind, history[seq_len(n_sweeps)]),
    backend = "R",
    method = "global_fixed_fit"
  )
}
