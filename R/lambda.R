# Modular smoothing / lambda engines (fixed, cGCV; reserved hooks: cFS, cREML).

#' Parse public `lambda` into an internal specification.
#'
#' @param lambda Scalar, length-`d` vector, or `"cGCV"` (also reserved `"cFS"`/`"cREML"`).
#' @param d Number of margins.
#' @param control Optional [tt_control()] (uses `lambda_start` for cGCV init).
#' @return List with `method`, `values`, `automatic`.
#' @keywords internal
parse_lambda_spec <- function(lambda, d, control = NULL) {
  d <- as.integer(d)
  start <- 1
  if (!is.null(control) && !is.null(control$lambda_start)) {
    start <- control$lambda_start
  }

  if (is.character(lambda)) {
    method <- match.arg(lambda, c("cGCV", "cFS", "cREML", "fixed"))
    if (method %in% c("cFS", "cREML")) {
      stop(
        "lambda = '", method,
        "' is not implemented yet (planned: TT-cFS / cREML).",
        call. = FALSE
      )
    }
    if (identical(method, "fixed")) {
      stop("Use a numeric lambda for fixed smoothing, or lambda = \"cGCV\".",
           call. = FALSE)
    }
    values <- rep(as.numeric(start), length.out = d)
    .validate_lambda_values(values, d)
    return(list(method = "cGCV", values = values, automatic = TRUE,
                lambda0 = values)) # lambda0 alias for older callers
  }

  if (is.numeric(lambda)) {
    if (length(lambda) == 1L) {
      values <- rep(as.numeric(lambda), d)
    } else if (length(lambda) == d) {
      values <- as.numeric(lambda)
    } else {
      stop(
        "Fixed anisotropic lambda must have length 1 (isotropic) or length d = ",
        d, "; got length ", length(lambda), ".",
        call. = FALSE
      )
    }
    .validate_lambda_values(values, d)
    return(list(method = "fixed", values = values, automatic = FALSE,
                lambda0 = values))
  }

  stop("lambda must be numeric or \"cGCV\".", call. = FALSE)
}

.validate_lambda_values <- function(values, d) {
  if (length(values) != d) stop("Internal lambda length mismatch.", call. = FALSE)
  if (any(!is.finite(values)) || any(values <= 0)) {
    stop("All lambda values must be finite and strictly positive.", call. = FALSE)
  }
  invisible(values)
}

#' Update one core's λ and coefficients under a smoothing method.
#'
#' Workspace must cache `S`, `b`, `P` (and weighted `Xw`,`yw` for GCV RSS).
#' Future releases may add `cFS` / `cREML` here without changing [ttps()].
#'
#' @keywords internal
update_lambda <- function(method, workspace, ...) {
  method <- match.arg(method, c("fixed", "cGCV", "cFS", "cREML"))
  switch(
    method,
    fixed = update_lambda_fixed(workspace, ...),
    cGCV = update_lambda_cgcv(workspace, ...),
    cFS = stop("cFS not implemented yet.", call. = FALSE),
    cREML = stop("cREML not implemented yet.", call. = FALSE)
  )
}

update_lambda_fixed <- function(workspace, ...) {
  lam <- workspace$lambda0
  P0 <- workspace$P0
  M <- workspace$S + lam * workspace$P
  if (!is.null(P0)) M <- M + P0
  # Global mode (P0 present): prefer exact SPD solve so fixed-λ ALS is an
  # exact conditional minimizer of Q (P4–P6). Own-margin / legacy keeps the
  # ridge path for R↔Rcpp parity with gaussian_core_update_cpp.
  if (!is.null(P0)) {
    g <- tryCatch(
      as.numeric(solve_spd(M, workspace$b)),
      error = function(e) solve_spd_ridge(M, workspace$b)
    )
    if (!all(is.finite(g))) g <- solve_spd_ridge(M, workspace$b)
  } else {
    g <- solve_spd_ridge(M, workspace$b)
  }
  # n_eval counts smoothing-*criterion* evaluations (cGCV/cFS/…), not core solves
  list(lambda = lam, g = g, value = NA_real_, ed = NA_real_, n_eval = 0L)
}

.ed_S <- function(S, P, lambda, P0 = NULL) {
  m <- nrow(S)
  ridge <- ridge_scale(S, multiplier = 1e-6)
  M <- S + lambda * P + ridge * diag(m)
  if (!is.null(P0)) M <- M + P0
  Minv_S <- tryCatch(
    solve_spd(M, S),
    error = function(e) {
      Mr <- S + lambda * P
      if (!is.null(P0)) Mr <- Mr + P0
      solve_spd_ridge(Mr, S, base_ridge = ridge)
    }
  )
  sum(diag(Minv_S))
}

.conditional_gcv <- function(yw, Xw, S, P, b, lambda, P0 = NULL) {
  n <- length(yw)
  M <- S + lambda * P
  if (!is.null(P0)) M <- M + P0
  g <- solve_spd_ridge(M, b)
  # Prefer quadratic form (same as spectral path); Xw kept for back-compat.
  rss <- as.numeric(sum(yw^2) - 2 * crossprod(b, g) + crossprod(g, S %*% g))
  if (!is.finite(rss)) {
    rss <- sum((yw - as.numeric(Xw %*% g))^2)
  }
  rss <- max(rss, 0)
  ed <- .ed_S(S, P, lambda, P0 = P0)
  denom <- (n - ed)^2
  value <- if (!is.finite(denom) || denom < 1e-12) Inf else n * rss / denom
  list(value = value, g = g, ed = ed, rss = rss)
}

#' Spectral factorization workspace for frozen conditional cGCV.
#'
#' For fixed \eqn{A = S + P_{-k} + \varepsilon I} and \eqn{B = P_{kk}},
#' factor once
#' \deqn{A = R^\top R,\quad
#'   R^{-\top} B R^{-1} = U\,\mathrm{diag}(d)\,U^\top}
#' and cache
#' \eqn{\tilde b = U^\top R^{-\top} b} and
#' \eqn{s_j = (U^\top R^{-\top} S R^{-1} U)_{jj}} so that each \eqn{\lambda}
#' evaluation is
#' \deqn{g(\lambda)=R^{-1}U(I+\lambda D)^{-1}\tilde b,\quad
#'   \mathrm{ed}(\lambda)=\sum_j s_j/(1+\lambda d_j),\quad
#'   \mathrm{RSS}=c-2b^\top g+g^\top S g.}
#'
#' @param S Gram matrix \eqn{X^\top W X}.
#' @param P Own-margin penalty \eqn{P_{kk}}.
#' @param b Cross-product \eqn{X^\top W z}.
#' @param yw Optional weighted response (for \eqn{c=z^\top W z} and \eqn{n}).
#' @param P0 Optional fixed offset \eqn{P_{k,-k}}.
#' @param max_m Skip factorization when \eqn{m>\texttt{max_m}}.
#' @return List workspace, or `NULL` if factorization fails / skipped.
#' @keywords internal
.make_gcv_spectral_cache <- function(S, P, b = NULL, yw = NULL, P0 = NULL,
                                     max_m = 1500L) {
  m <- nrow(S)
  if (m < 1L || m > as.integer(max_m)) return(NULL)
  if (is.null(b)) {
    # Back-compat: older callers passed only (S, P, P0).
    b <- rep(0, m)
  }
  b <- as.numeric(b)
  # Match .ed_S ridge scale so spectral edf ≡ tr((A+λB)^{-1} S).
  ridge <- ridge_scale(S, multiplier = 1e-6)
  A <- S + ridge * diag(m)
  if (!is.null(P0)) A <- A + P0
  R <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  # C = R^{-T} P R^{-1}: first Y = P R^{-1}, then C = R^{-T} Y.
  # (Previously used R^{-1} R^{-T} P = A^{-1} P, which is the wrong congruence.)
  Yt <- forwardsolve(t(R), t(P)) # Yt = R^{-T} P^T ⇒ Y = P R^{-1}
  Y <- t(Yt)
  C <- forwardsolve(t(R), Y)
  C <- (C + t(C)) / 2
  eg <- tryCatch(eigen(C, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eg)) return(NULL)
  U <- Re(eg$vectors)
  dvals <- pmax(Re(eg$values), 0)
  # transformed_b = U^T R^{-T} b
  Rb <- forwardsolve(t(R), b)
  transformed_b <- as.numeric(crossprod(U, Rb))
  # ed_weights = diag(U^T R^{-T} S R^{-1} U)
  St_t <- forwardsolve(t(R), t(S)) # St = S R^{-1}
  St <- t(St_t)
  St <- forwardsolve(t(R), St)     # St = R^{-T} S R^{-1}
  St <- (St + t(St)) / 2
  SU <- St %*% U
  ed_weights <- as.numeric(colSums(U * SU))
  zWz <- if (is.null(yw)) NA_real_ else sum(as.numeric(yw)^2)
  n <- if (is.null(yw)) NA_integer_ else length(yw)
  list(
    R = R,
    values = dvals,
    vectors = U,
    transformed_b = transformed_b,
    ed_weights = ed_weights,
    zWz = zWz,
    S = S,
    b = b,
    ridge = ridge,
    n = n
  )
}

#' Evaluate conditional GCV from a spectral workspace (no n×m matvec).
#' @keywords internal
.conditional_gcv_spectral <- function(yw, Xw, S, P, b, lambda, cache,
                                      P0 = NULL) {
  if (is.null(cache)) {
    return(.conditional_gcv(yw, Xw, S, P, b, lambda, P0 = P0))
  }
  n <- cache$n
  if (is.null(n) || !is.finite(n)) n <- length(yw)
  inv <- 1 / (1 + lambda * cache$values)
  coef_u <- cache$transformed_b * inv
  g <- as.numeric(backsolve(cache$R, cache$vectors %*% coef_u))
  ed <- sum(cache$ed_weights * inv)
  zWz <- cache$zWz
  if (!is.finite(zWz)) zWz <- sum(as.numeric(yw)^2)
  bb <- cache$b
  SS <- cache$S
  rss <- as.numeric(zWz - 2 * crossprod(bb, g) + crossprod(g, SS %*% g))
  if (!is.finite(rss) && !is.null(Xw)) {
    rss <- sum((yw - as.numeric(Xw %*% g))^2)
  }
  rss <- max(rss, 0)
  denom <- (n - ed)^2
  value <- if (!is.finite(denom) || denom < 1e-12) Inf else n * rss / denom
  list(value = value, g = drop(g), ed = ed, rss = rss)
}

#' Resolve / attach spectral cache on a core workspace (build at most once).
#' @keywords internal
.cgcv_spectral_from_workspace <- function(workspace) {
  if (!isTRUE(workspace$use_spectral)) return(NULL)
  cache <- workspace$spectral
  if (!is.null(cache)) return(cache)
  .make_gcv_spectral_cache(
    workspace$S, workspace$P, workspace$b,
    yw = workspace$yw, P0 = workspace$P0
  )
}

update_lambda_cgcv <- function(workspace, ...) {
  bounds <- workspace$bounds
  tol <- workspace$tol
  yw <- workspace$yw
  Xw <- workspace$Xw
  S <- workspace$S
  P <- workspace$P
  P0 <- workspace$P0
  b <- workspace$b
  cache <- .cgcv_spectral_from_workspace(workspace)
  n_eval <- 0L
  obj <- function(ll) {
    n_eval <<- n_eval + 1L
    lam <- exp(ll)
    if (is.null(cache)) {
      .conditional_gcv(yw, Xw, S, P, b, lam, P0 = P0)$value
    } else {
      .conditional_gcv_spectral(yw, Xw, S, P, b, lam, cache, P0 = P0)$value
    }
  }
  opt <- stats::optimize(obj, interval = log(bounds), tol = tol)
  lam <- exp(opt$minimum)
  fit <- if (is.null(cache)) {
    .conditional_gcv(yw, Xw, S, P, b, lam, P0 = P0)
  } else {
    .conditional_gcv_spectral(yw, Xw, S, P, b, lam, cache, P0 = P0)
  }
  list(
    lambda = lam, g = fit$g, value = fit$value, ed = fit$ed,
    n_eval = n_eval + 1L
  )
}

#' Build cached conditional workspace for core k (Gaussian or weighted).
#' @param P0 Optional fixed penalty offset (for exact P_k^full cross-margin terms).
#' @keywords internal
make_core_workspace <- function(zc, X, P, lambda0, bounds, tol,
                                weight = NULL, use_spectral = FALSE,
                                P0 = NULL) {
  if (is.null(weight)) {
    S <- crossprod(X)
    b <- as.numeric(crossprod(X, zc))
    yw <- zc
    Xw <- X
  } else {
    sw <- sqrt(pmax(as.numeric(weight), 0))
    Xw <- X * sw
    yw <- zc * sw
    S <- crossprod(Xw)
    b <- as.numeric(crossprod(Xw, yw))
  }
  ws <- list(
    S = S, b = b, P = P, P0 = P0, X = X, Xw = Xw, yw = yw,
    lambda0 = lambda0, bounds = bounds, tol = tol,
    use_spectral = isTRUE(use_spectral),
    spectral = NULL
  )
  if (isTRUE(use_spectral)) {
    ws$spectral <- .make_gcv_spectral_cache(S, P, b, yw = yw, P0 = P0)
  }
  ws
}

#' Classify each λ relative to cGCV search bounds.
#'
#' Labels: `"lower"`, `"upper"`, `"interior"` (or `"NA"` if non-finite).
#' Nearness is judged on a multiplicative / log scale (default 2%).
#'
#' @param lambda Numeric vector of smoothing parameters.
#' @param bounds Length-2 `c(lambda_min, lambda_max)`.
#' @param rel Relative tolerance (e.g. `0.02` ⇒ within 2% of a bound).
#' @return Character vector of labels, same length as `lambda`.
#' @keywords internal
#' @noRd
.lambda_boundary_status <- function(lambda, bounds, rel = 0.02) {
  lambda <- as.numeric(lambda)
  bounds <- as.numeric(bounds)
  if (length(bounds) != 2L || any(!is.finite(bounds)) || bounds[1] <= 0 ||
      bounds[2] <= bounds[1]) {
    return(rep(NA_character_, length(lambda)))
  }
  lo <- bounds[1]
  hi <- bounds[2]
  log_tol <- abs(log1p(rel))
  vapply(lambda, function(l) {
    if (!is.finite(l) || l <= 0) return(NA_character_)
    if (abs(log(l) - log(lo)) <= log_tol || l <= lo * (1 + rel)) return("lower")
    if (abs(log(l) - log(hi)) <= log_tol || l >= hi / (1 + rel)) return("upper")
    "interior"
  }, character(1))
}

#' Attach λ-boundary diagnostics and optionally warn.
#' @keywords internal
#' @noRd
.tt_lambda_boundary_info <- function(lambda, method, control) {
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  status <- .lambda_boundary_status(lambda, bounds)
  at_bound <- isTRUE(identical(method, "cGCV")) &&
    any(status %in% c("lower", "upper"), na.rm = TRUE)
  if (at_bound && isTRUE(control$warn_lambda_boundary %||% TRUE)) {
    idx <- which(status %in% c("lower", "upper"))
    parts <- sprintf("lambda[%d]->%s (%.6g)", idx, status[idx], lambda[idx])
    warning(
      "Note: ", paste(parts, collapse = "; "),
      " near the cGCV search boundaries [",
      format(bounds[1], digits = 4), ", ", format(bounds[2], digits = 4), "]. ",
      "Interpret directional λ cautiously; consider widening lambda_bounds ",
      "or inspecting the surface.",
      call. = FALSE
    )
  }
  list(
    lambda_bounds = as.numeric(bounds),
    lambda_boundary = status,
    lambda_at_boundary = isTRUE(at_bound)
  )
}
