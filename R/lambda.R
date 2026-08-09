# Modular smoothing / λ engines (Paper 1: fixed, cGCV; Paper 2 hooks: cFS, cREML).

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
        "' is reserved for Paper 2 (TT-cFS / cREML) and not implemented yet.",
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
#' Paper 2 may add `cFS` / `cREML` here without changing [ttpspline()].
#'
#' @keywords internal
update_lambda <- function(method, workspace, ...) {
  method <- match.arg(method, c("fixed", "cGCV", "cFS", "cREML"))
  switch(
    method,
    fixed = update_lambda_fixed(workspace, ...),
    cGCV = update_lambda_cgcv(workspace, ...),
    cFS = stop("cFS not implemented (Paper 2).", call. = FALSE),
    cREML = stop("cREML not implemented (Paper 2).", call. = FALSE)
  )
}

update_lambda_fixed <- function(workspace, ...) {
  lam <- workspace$lambda0
  g <- solve_spd_ridge(workspace$S + lam * workspace$P, workspace$b)
  list(lambda = lam, g = g, value = NA_real_, ed = NA_real_, n_eval = 1L)
}

.ed_S <- function(S, P, lambda) {
  m <- nrow(S)
  ridge <- ridge_scale(S, multiplier = 1e-6)
  M <- S + lambda * P + ridge * diag(m)
  Minv_S <- tryCatch(
    solve_spd(M, S),
    error = function(e) solve_spd_ridge(S + lambda * P, S, base_ridge = ridge)
  )
  sum(diag(Minv_S))
}

.conditional_gcv <- function(yw, Xw, S, P, b, lambda) {
  n <- length(yw)
  g <- solve_spd_ridge(S + lambda * P, b)
  rss <- sum((yw - as.numeric(Xw %*% g))^2)
  ed <- .ed_S(S, P, lambda)
  denom <- (n - ed)^2
  value <- if (!is.finite(denom) || denom < 1e-12) Inf else n * rss / denom
  list(value = value, g = g, ed = ed, rss = rss)
}

#' Optional spectral cache for repeated GCV evals on fixed (S, P).
#' Uses chol(S_ridge) then eigen of transformed P (dense, small m).
#' @keywords internal
.make_gcv_spectral_cache <- function(S, P) {
  m <- nrow(S)
  if (m > 400L) return(NULL) # not worth for large cores
  ridge <- ridge_scale(S, multiplier = 1e-8)
  Sr <- S + ridge * diag(m)
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

.conditional_gcv_spectral <- function(yw, Xw, S, P, b, lambda, cache) {
  # Fall back if cache missing
  if (is.null(cache)) return(.conditional_gcv(yw, Xw, S, P, b, lambda))
  n <- length(yw)
  R <- cache$R
  # Solve (S_ridge + lam P) g = b via spectral of Q
  # (R'R + lam P) g = b  <=>  (I + lam Q) R g = R^{-T} b
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
  b <- workspace$b
  use_spectral <- isTRUE(workspace$use_spectral)
  cache <- if (use_spectral) .make_gcv_spectral_cache(S, P) else NULL
  n_eval <- 0L
  obj <- function(ll) {
    n_eval <<- n_eval + 1L
    lam <- exp(ll)
    if (is.null(cache)) {
      .conditional_gcv(yw, Xw, S, P, b, lam)$value
    } else {
      .conditional_gcv_spectral(yw, Xw, S, P, b, lam, cache)$value
    }
  }
  opt <- stats::optimize(obj, interval = log(bounds), tol = tol)
  lam <- exp(opt$minimum)
  fit <- if (is.null(cache)) {
    .conditional_gcv(yw, Xw, S, P, b, lam)
  } else {
    .conditional_gcv_spectral(yw, Xw, S, P, b, lam, cache)
  }
  list(
    lambda = lam, g = fit$g, value = fit$value, ed = fit$ed,
    n_eval = n_eval + 1L
  )
}

#' Build cached conditional workspace for core k (Gaussian or weighted).
#' @keywords internal
make_core_workspace <- function(zc, X, P, lambda0, bounds, tol,
                                weight = NULL, use_spectral = FALSE) {
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
    S = S, b = b, P = P, X = X, Xw = Xw, yw = yw,
    lambda0 = lambda0, bounds = bounds, tol = tol,
    use_spectral = isTRUE(use_spectral)
  )
}
