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
  rss <- sum((yw - as.numeric(Xw %*% g))^2)
  ed <- .ed_S(S, P, lambda, P0 = P0)
  denom <- (n - ed)^2
  value <- if (!is.finite(denom) || denom < 1e-12) Inf else n * rss / denom
  list(value = value, g = g, ed = ed, rss = rss)
}

#' Optional spectral cache for repeated GCV evals on fixed (S, P[, P0]).
#' Uses chol(S_ridge + P0) then eigen of transformed P (dense, small m).
#' @keywords internal
.make_gcv_spectral_cache <- function(S, P, P0 = NULL) {
  m <- nrow(S)
  if (m > 400L) return(NULL) # not worth for large cores
  ridge <- ridge_scale(S, multiplier = 1e-8)
  Sr <- S + ridge * diag(m)
  if (!is.null(P0)) Sr <- Sr + P0
  R <- tryCatch(chol(Sr), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  # Q = R^{-T} P R^{-1}
  tmp <- forwardsolve(t(R), P)
  Q <- backsolve(R, tmp)
  Q <- (Q + t(Q)) / 2
  eg <- tryCatch(eigen(Q, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eg)) return(NULL)
  list(R = R, values = pmax(Re(eg$values), 0), vectors = Re(eg$vectors),
       ridge = ridge)
}

.conditional_gcv_spectral <- function(yw, Xw, S, P, b, lambda, cache, P0 = NULL) {
  # Fall back if cache missing
  if (is.null(cache)) return(.conditional_gcv(yw, Xw, S, P, b, lambda, P0 = P0))
  n <- length(yw)
  R <- cache$R
  # Solve (S_ridge + P0 + lam P) g = b via spectral of Q
  rhs <- forwardsolve(t(R), b)
  U <- cache$vectors
  d <- cache$values
  coef_u <- as.numeric(crossprod(U, rhs)) / (1 + lambda * d)
  Rg <- U %*% coef_u
  g <- backsolve(R, Rg)
  rss <- sum((yw - as.numeric(Xw %*% g))^2)
  # ed ≈ sum 1/(1+lam d_i) for the whitened problem (approx with ridge)
  ed <- sum(1 / (1 + lambda * d))
  denom <- (n - ed)^2
  value <- if (!is.finite(denom) || denom < 1e-12) Inf else n * rss / denom
  list(value = value, g = drop(g), ed = ed, rss = rss)
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
  use_spectral <- isTRUE(workspace$use_spectral)
  cache <- if (use_spectral) .make_gcv_spectral_cache(S, P, P0 = P0) else NULL
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
  list(
    S = S, b = b, P = P, P0 = P0, X = X, Xw = Xw, yw = yw,
    lambda0 = lambda0, bounds = bounds, tol = tol,
    use_spectral = isTRUE(use_spectral)
  )
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
