# Conditional Marra–Wood style uncertainty for TT P-splines.
# Level 1 only: Var{f̂(x) | r̂, λ̂}. Gauge handled by left-orthogonal cores.
# See docs/TT_INFERENCE_MARRA_WOOD.md.

#' Left-orthogonal TT gauge fix (sequential QR on left unfoldings).
#'
#' Absorbs triangular factors into the next core so the contracted map is
#' unchanged. Used as the local identifiable parameterization for inference.
#'
#' @keywords internal
#' @noRd
tt_left_orthogonalize <- function(cores) {
  d <- length(cores)
  out <- lapply(cores, function(g) {
    a <- array(as.numeric(g), dim(g))
    storage.mode(a) <- "double"
    a
  })
  if (d < 2L) return(out)
  for (k in seq_len(d - 1L)) {
    dm <- dim(out[[k]])
    rl <- dm[1L]
    p <- dm[2L]
    rr <- dm[3L]
    M <- matrix(out[[k]], rl * p, rr)
    if (rr < 1L || nrow(M) < rr) next
    # Use non-pivoted QR so R stays upper-triangular and M = Q %*% R.
    # (LAPACK pivoted QR would require an explicit unpivot absorb.)
    qrM <- qr(M, LAPACK = FALSE)
    Q <- qr.Q(qrM, complete = FALSE)
    R <- qr.R(qrM)
    # Sign convention for uniqueness of the QR gauge
    s <- sign(diag(R))
    s[!is.finite(s) | s == 0] <- 1
    S <- diag(s, nrow = length(s))
    Q <- Q %*% S
    R <- S %*% R
    out[[k]] <- array(Q, c(rl, p, rr))
    gnext <- out[[k + 1L]]
    dn <- dim(gnext)
    # Matricize next core as (left-rank) × (p * right-rank)
    Mn <- matrix(gnext, dn[1L], dn[2L] * dn[3L])
    if (nrow(R) != nrow(Mn)) {
      stop("Gauge absorb size mismatch at core ", k, ".", call. = FALSE)
    }
    out[[k + 1L]] <- array(R %*% Mn, dn)
  }
  out
}

#' Apply an invertible gauge transform at interface `iface` (1..d-1).
#'
#' G_k <- G_k · A (right index), G_{k+1} <- A^{-1} · G_{k+1} (left index).
#' Contracted coefficient tensor is unchanged.
#'
#' @keywords internal
#' @noRd
tt_apply_gauge <- function(cores, iface, A) {
  d <- length(cores)
  stopifnot(iface >= 1L, iface <= d - 1L)
  A <- as.matrix(A)
  Ai <- solve(A)
  out <- lapply(cores, function(g) array(as.numeric(g), dim(g)))
  gk <- out[[iface]]
  dm <- dim(gk)
  # Right multiply: for each (a,j), row vec of length rr *= A
  M <- matrix(aperm(gk, c(1L, 2L, 3L)), dm[1L] * dm[2L], dm[3L])
  M <- M %*% A
  out[[iface]] <- array(M, dm)
  gkp <- out[[iface + 1L]]
  dn <- dim(gkp)
  Mn <- matrix(gkp, dn[1L], dn[2L] * dn[3L])
  Mn <- Ai %*% Mn
  out[[iface + 1L]] <- array(Mn, dn)
  out
}

#' Rebuild per-core P-spline penalties from a fit object.
#' @keywords internal
#' @noRd
.tt_fit_penalties <- function(object) {
  if (!is.null(object$penalties)) return(object$penalties)
  ranks <- as.integer(object$rank)
  p <- as.integer(object$k)
  d <- as.integer(object$d)
  po <- as.integer(object$penalty_order %||% 2L)
  tt_core_penalties(ranks, p, po, cyclic = object$cyclic)
}

#' Training bases for a fitted object.
#' @keywords internal
#' @noRd
.tt_fit_basis <- function(object) {
  if (!is.null(object$basis)) return(object$basis)
  if (is.null(object$X)) {
    stop(
      "Training covariates fit$X are required for inference. ",
      "Refit with the current package version.",
      call. = FALSE
    )
  }
  eval_marginal_bases(object$X, object$knots, object$degree,
                      cyclic = object$cyclic)
}

#' Row Jacobian of η = α + f_TT at newdata (columns: intercept + packed cores).
#'
#' Uses the same stacked conditional-design map as joint EDF, evaluated at
#' gauge-fixed cores.
#'
#' @keywords internal
#' @noRd
.tt_predict_jacobian <- function(cores, intercept_col = TRUE, basis) {
  Jf <- tt_stacked_jacobian(cores, basis, weight = NULL)
  if (!isTRUE(intercept_col)) return(Jf)
  cbind(1, Jf)
}

#' Finite-difference check of prediction Jacobian (one newdata row).
#' @keywords internal
#' @noRd
.tt_predict_jacobian_fd <- function(cores, intercept, basis_row, eps = 1e-6) {
  theta <- .tt_pack_cores(cores)
  template <- cores
  f0 <- as.numeric(intercept + tt_contraction(cores, basis_row))
  npar <- length(theta)
  J <- numeric(npar)
  for (j in seq_len(npar)) {
    tp <- theta
    tm <- theta
    tp[j] <- tp[j] + eps
    tm[j] <- tm[j] - eps
    fp <- as.numeric(intercept + tt_contraction(.tt_unpack_cores(tp, template), basis_row))
    fm <- as.numeric(intercept + tt_contraction(.tt_unpack_cores(tm, template), basis_row))
    J[j] <- (fp - fm) / (2 * eps)
  }
  c(1, J) # intercept column
}

#' Build / refresh conditional inference factorization on a fit.
#'
#' Gaussian Gate 1/2: H = X'X with X = [1 | J], S = blkdiag(0, S_λ),
#' H_p = H + S + ridge, V_B = σ² H_p^{-1}, V_F = σ² H_p^{-1} H H_p^{-1}.
#' Cores are left-orthogonalized before forming J (gauge-fixed local coords).
#'
#' @keywords internal
#' @noRd
tt_prepare_inference <- function(object, force = FALSE) {
  stopifnot(inherits(object, "ttpspline"))
  cache <- object$._inf
  if (is.environment(cache) && !isTRUE(force) && isTRUE(cache$ready) &&
      !is.null(cache$data)) {
    object$inference <- cache$data
    return(object)
  }
  if (!isTRUE(force) && !is.null(object$inference) &&
      isTRUE(object$inference$ready)) {
    return(object)
  }
  key <- object$family_key %||% family_key(object$family)
  if (!identical(key, "gaussian")) {
    stop(
      "Pointwise uncertainty is implemented for Gaussian fits first ",
      "(Gate 1/2). GLM extension is not enabled yet.",
      call. = FALSE
    )
  }
  t0 <- proc.time()[[3L]]
  basis <- .tt_fit_basis(object)
  penalties <- .tt_fit_penalties(object)
  cores0 <- object$cores
  cores <- tt_left_orthogonalize(cores0)
  # Sanity: gauge fix must preserve fitted surface
  f0 <- tt_contraction(cores0, basis)
  f1 <- tt_contraction(cores, basis)
  if (max(abs(f0 - f1)) > 1e-8 * (1 + max(abs(f0)))) {
    stop("Left-orthogonal gauge fix changed the TT contraction.", call. = FALSE)
  }

  npar <- sum(vapply(cores, length, integer(1)))
  n <- length(object$y)
  max_npar <- as.integer(object$control$edf_max_npar %||% 2500L)
  if (npar > max_npar) {
    stop(
      sprintf(
        "Inference skipped: npar_TT=%d exceeds edf_max_npar=%d.",
        npar, max_npar
      ),
      call. = FALSE
    )
  }
  if (as.numeric(n) * as.numeric(npar) > 5e7) {
    stop("Inference skipped: dense Jacobian would be too large.", call. = FALSE)
  }

  J <- tt_stacked_jacobian(cores, basis, weight = NULL)
  X <- cbind(1, J)
  S_cores <- tt_block_penalty(penalties, as.numeric(object$lambda))
  m <- ncol(X)
  S <- matrix(0, m, m)
  S[2:m, 2:m] <- S_cores

  H <- crossprod(X)
  ridge <- ridge_scale(H, multiplier = 1e-9)
  Hp <- H + S + ridge * diag(m)

  chol_Hp <- tryCatch(chol(Hp), error = function(e) NULL)
  if (is.null(chol_Hp)) {
    # Escalate ridge (gauge near-null directions)
    for (fac in c(1e2, 1e3, 1e4, 1e5, 1e6)) {
      Hp_try <- H + S + (fac * ridge) * diag(m)
      chol_Hp <- tryCatch(chol(Hp_try), error = function(e) NULL)
      if (!is.null(chol_Hp)) {
        Hp <- Hp_try
        ridge <- fac * ridge
        break
      }
    }
  }
  if (is.null(chol_Hp)) {
    stop(
      "Penalized Hessian H_p is not numerically positive definite after ",
      "gauge fixing and ridge. Inference aborted (STOP condition).",
      call. = FALSE
    )
  }

  # EDF of augmented smoother (includes intercept)
  infl <- solve_spd(Hp, H)
  edf_inf <- sum(diag(infl))
  if (!is.finite(edf_inf) || edf_inf < 1) {
    edf_inf <- if (is.finite(object$edf)) object$edf + 1 else 2
  }

  off <- object$offset %||% rep(0, n)
  eta <- as.numeric(off + object$intercept + f1)
  rss <- sum((object$y - eta)^2)
  df_res <- max(n - edf_inf, 1)
  sigma2 <- rss / df_res
  if (!is.finite(sigma2) || sigma2 <= 0) {
    sigma2 <- max(rss / max(n - 1, 1), .Machine$double.eps)
  }

  setup_time <- proc.time()[[3L]] - t0
  object$penalties <- penalties
  inf <- list(
    ready = TRUE,
    method = "bayesian_penalized_spline",
    parameterization = "left_orthogonal_TT_cores_plus_intercept",
    gauge = "left-orthogonal QR (sequential)",
    family = "gaussian",
    conditional_on_rank = TRUE,
    conditional_on_lambda = TRUE,
    smoothing_uncertainty = FALSE,
    rank_uncertainty = FALSE,
    scale = sigma2,
    scale_estimator = "RSS / (n - edf_aug)",
    edf_aug = edf_inf,
    ridge = ridge,
    npar_packed = npar,
    npar_aug = m,
    npar_intrinsic = object$npar_tt_intrinsic %||% (npar - tt_gauge_dim(object$rank)),
    dim_note = paste0(
      "Factorization lives in packed TT coordinates after left-orthogonal ",
      "gauge fix (size npar_TT+1 including intercept). Identifiable ",
      "directions align with dim(M_r)=npar_TT-sum r_k^2; gauge null ",
      "directions are controlled by the penalty + ridge."
    ),
    cores_gauge = cores,
    chol_Hp = chol_Hp,
    H_data = H,
    H_pen = S,
    setup_time = setup_time,
    lambda_method = object$lambda_method %||% "fixed"
  )
  object$inference <- inf
  if (is.environment(cache)) {
    cache$ready <- TRUE
    cache$data <- inf
  }
  object
}

#' Solve H_p^{-1} B using stored Cholesky factor.
#' @keywords internal
#' @noRd
.tt_Hp_solve <- function(inference, B) {
  R <- inference$chol_Hp
  if (is.null(dim(B))) B <- matrix(B, ncol = 1L)
  backsolve(R, forwardsolve(t(R), B))
}

#' Pointwise SE of η at rows of Jacobian J_new (m × npar_aug).
#' @keywords internal
#' @noRd
.tt_se_from_jacobian <- function(inference, J_new,
                                 type = c("bayesian", "frequentist")) {
  type <- match.arg(type)
  sigma2 <- inference$scale
  # W = H_p^{-1} J'  (npar × m)
  W <- .tt_Hp_solve(inference, t(J_new))
  if (identical(type, "bayesian")) {
    # diag(J H_p^{-1} J') = colSums(J' * W) after rearrange: rowSums(J_new * t(W))
    q <- rowSums(J_new * t(W))
  } else {
    # V_F ∝ H_p^{-1} H H_p^{-1}; q_i = w_i' H w_i
    H <- inference$H_data
    q <- vapply(seq_len(ncol(W)), function(i) {
      wi <- W[, i]
      as.numeric(crossprod(wi, H %*% wi))
    }, numeric(1))
  }
  q <- pmax(as.numeric(q), 0)
  sqrt(sigma2 * q)
}

#' Approximate Bayesian / frequentist covariance of packed parameters.
#'
#' Default is the Bayesian penalized-spline covariance
#' \eqn{V_B=\hat\sigma^2 H_p^{-1}} with
#' \eqn{H_p=H+S_\lambda} in left-orthogonal TT coordinates (plus intercept).
#' Frequentist sandwich \eqn{V_F=\hat\sigma^2 H_p^{-1} H H_p^{-1}} is optional.
#'
#' Uncertainty is **conditional** on the fitted TT rank and smoothing
#' parameters. `unconditional = TRUE` is not implemented.
#'
#' @param object A `"ttpspline"` fit (Gaussian for v1).
#' @param type `"bayesian"` (default) or `"frequentist"`.
#' @param unconditional Must be `FALSE` in v1.
#' @param ... Unused.
#' @return Dense covariance matrix for `(intercept, packed TT cores)` in the
#'   gauge-fixed parameterization stored on `object$inference`.
#' @export
vcov.ttpspline <- function(object,
                           type = c("bayesian", "frequentist"),
                           unconditional = FALSE,
                           ...) {
  type <- match.arg(type)
  if (isTRUE(unconditional)) {
    stop(
      "smoothing-parameter uncertainty not implemented ",
      "(unconditional = TRUE). v1 inference is conditional on fitted lambda.",
      call. = FALSE
    )
  }
  object <- tt_prepare_inference(object)
  inf <- object$inference
  m <- inf$npar_aug
  # Columns of I: solve H_p^{-1}
  Hp_inv <- .tt_Hp_solve(inf, diag(m))
  if (identical(type, "bayesian")) {
    return(inf$scale * Hp_inv)
  }
  H <- inf$H_data
  # H_p^{-1} H H_p^{-1}
  mid <- H %*% Hp_inv
  inf$scale * (Hp_inv %*% mid)
}

#' @rdname predict.ttpspline
#' @param se.fit If `TRUE`, return list with `fit` and `se.fit` (link-scale SE).
#' @param interval `"none"` or `"confidence"` (pointwise, conditional).
#' @param level Confidence level for intervals.
#' @param vcov_type `"bayesian"` (default for `se.fit`) or `"frequentist"`.
#' @param full_cov Reserved; must be `FALSE` (no m×m prediction covariance).
#' @export
predict.ttpspline <- function(object,
                              newdata = NULL,
                              type = c("link", "response"),
                              offset = NULL,
                              se.fit = FALSE,
                              interval = c("none", "confidence"),
                              level = 0.95,
                              vcov_type = c("bayesian", "frequentist"),
                              full_cov = FALSE,
                              ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  vcov_type <- match.arg(vcov_type)
  if (isTRUE(full_cov)) {
    stop("full_cov = TRUE is not implemented; use diagonal SE only.", call. = FALSE)
  }
  want_se <- isTRUE(se.fit) || identical(interval, "confidence")

  if (is.null(newdata)) {
    eta <- object$linear.predictors
    Xnew <- object$X
    if (want_se && is.null(Xnew)) {
      stop("newdata or fit$X required when se.fit / interval is requested.",
           call. = FALSE)
    }
  } else {
    Xnew <- as.matrix(newdata)
    if (ncol(Xnew) != object$d) {
      stop("newdata must have ", object$d, " columns.", call. = FALSE)
    }
    if (anyNA(Xnew)) stop("NA in newdata not supported.", call. = FALSE)
    off <- normalize_offset(offset, nrow(Xnew))
    basis <- eval_marginal_bases(Xnew, object$knots, object$degree,
                                 cyclic = object$cyclic)
    eta <- tt_eta(off, object$intercept, object$cores, basis)
  }

  if (!want_se) {
    if (identical(type, "link") || identical(object$family_key, "gaussian")) {
      return(eta)
    }
    return(invlink_eta(object$family, eta))
  }

  object <- tt_prepare_inference(object)
  cores_g <- object$inference$cores_gauge
  if (is.null(newdata)) {
    basis_new <- .tt_fit_basis(object)
  } else {
    basis_new <- eval_marginal_bases(Xnew, object$knots, object$degree,
                                     cyclic = object$cyclic)
  }
  J_new <- .tt_predict_jacobian(cores_g, intercept_col = TRUE, basis = basis_new)
  se_eta <- .tt_se_from_jacobian(object$inference, J_new, type = vcov_type)

  z <- stats::qnorm(1 - (1 - level) / 2)
  lower_eta <- eta - z * se_eta
  upper_eta <- eta + z * se_eta

  if (identical(type, "link") || identical(object$family_key, "gaussian")) {
    fit_out <- eta
    lower <- lower_eta
    upper <- upper_eta
    se_out <- se_eta
  } else {
    # Transform intervals on the link scale (not SE on the response)
    fit_out <- invlink_eta(object$family, eta)
    lower <- invlink_eta(object$family, lower_eta)
    upper <- invlink_eta(object$family, upper_eta)
    se_out <- se_eta # still link-scale; documented
  }

  out <- list(fit = fit_out, se.fit = se_out)
  if (identical(interval, "confidence")) {
    out$lower <- lower
    out$upper <- upper
    out$level <- level
    out$interval <- "confidence"
    out$vcov_type <- vcov_type
    out$scale <- "link"
  }
  out
}

#' Truncate TT cores to a lower (or equal) rank by sequential SVD.
#'
#' Diagnostic / warm-start helper: left-orthogonalize, then truncate bond
#' dimensions right-to-left to the target [tt_rank()] chain. The contracted
#' coefficient map is approximated (not exact unless singular values beyond
#' the cut are zero). Use to warm-start a lower-rank [ttps()] fit from a
#' higher-rank solution, e.g. Ishigami \(r=3\to r=2\).
#'
#' @param cores List of TT core arrays (length `d`), as in a `"ttpspline"` fit
#'   (`fit$cores`), or from [tt_initialize()].
#' @param rank Target rank (scalar, length `d-1`, or full chain); see [tt_rank()].
#' @return List of truncated cores with attribute `ranks`.
#' @examples
#' init3 <- tt_initialize(d = 3, rank = 3, k = 6, seed = 1)
#' init2 <- tt_truncate_rank(init3, rank = 2)
#' attr(init2, "ranks")
#' @export
tt_truncate_rank <- function(cores, rank) {
  if (!is.list(cores) || length(cores) < 2L) {
    stop("`cores` must be a length-d list of TT cores.", call. = FALSE)
  }
  d <- length(cores)
  cores <- lapply(cores, function(g) {
    a <- array(as.numeric(g), dim(g))
    storage.mode(a) <- "double"
    a
  })
  target <- tt_rank(rank, d = d)
  cores <- tt_left_orthogonalize(cores)

  # Right-to-left bond truncation to target left-ranks of cores 2..d
  for (k in d:2) {
    r_keep <- as.integer(target[k])
    gk <- cores[[k]]
    dm <- dim(gk)
    rl <- dm[1L]
    p <- dm[2L]
    rr <- dm[3L]
    if (r_keep >= rl) next
    if (r_keep < 1L) {
      stop("Target rank chain has non-positive bond at position ", k, ".",
           call. = FALSE)
    }
    # Matricize as left-bond × (p * right-bond)
    M <- matrix(gk, rl, p * rr)
    r_keep <- min(r_keep, nrow(M), ncol(M))
    sv <- svd(M, nu = r_keep, nv = r_keep)
    US <- sv$u %*% diag(sv$d[seq_len(r_keep)], nrow = r_keep)
    # New core k: r_keep × p × rr from V
    cores[[k]] <- array(t(sv$v), c(r_keep, p, rr))
    # Absorb U S into previous core's right bond
    gprev <- cores[[k - 1L]]
    dp <- dim(gprev)
    if (dp[3L] != rl) {
      stop("Rank mismatch absorbing into core ", k - 1L, ".", call. = FALSE)
    }
    Mp <- matrix(gprev, dp[1L] * dp[2L], dp[3L])
    Mp <- Mp %*% US
    cores[[k - 1L]] <- array(Mp, c(dp[1L], dp[2L], r_keep))
  }

  ranks <- integer(d + 1L)
  ranks[1L] <- dim(cores[[1L]])[1L]
  for (k in seq_len(d)) ranks[k + 1L] <- dim(cores[[k]])[3L]
  if (ranks[1L] != 1L || ranks[d + 1L] != 1L) {
    warning("Truncated cores do not have boundary ranks 1; check SVD path.",
            call. = FALSE)
  }
  attr(cores, "ranks") <- as.integer(ranks)
  cores
}
