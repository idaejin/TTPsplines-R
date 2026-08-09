# Modular smoothing / λ engines (Paper 1: fixed, cGCV; Paper 2 hooks: cFS, cREML).

parse_lambda_spec <- function(lambda, d) {
  if (is.character(lambda)) {
    method <- match.arg(lambda, c("cGCV", "cFS", "cREML", "fixed"))
    if (method %in% c("cFS", "cREML")) {
      stop(
        "lambda = '", method, "' is reserved for Paper 2 and not implemented in v0.",
        call. = FALSE
      )
    }
    if (identical(method, "fixed")) {
      stop("Use a numeric lambda for fixed smoothing, or lambda = \"cGCV\".", call. = FALSE)
    }
    return(list(method = "cGCV", lambda0 = rep(1, d)))
  }
  if (is.numeric(lambda)) {
    return(list(method = "fixed", lambda0 = rep(as.numeric(lambda), length.out = d)))
  }
  stop("lambda must be numeric or \"cGCV\".", call. = FALSE)
}

#' Update one core's λ and coefficients under a smoothing method.
#'
#' Workspace fields used: yc_or_zc, X, S, b, P, lambda0, bounds, tol, weight (optional).
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
  g <- solve_spd_ridge(
    workspace$S + lam * workspace$P,
    workspace$b
  )
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

update_lambda_cgcv <- function(workspace, ...) {
  bounds <- workspace$bounds
  tol <- workspace$tol
  yw <- workspace$yw
  Xw <- workspace$Xw
  S <- workspace$S
  P <- workspace$P
  b <- workspace$b
  n_eval <- 0L
  obj <- function(ll) {
    n_eval <<- n_eval + 1L
    .conditional_gcv(yw, Xw, S, P, b, exp(ll))$value
  }
  opt <- stats::optimize(obj, interval = log(bounds), tol = tol)
  lam <- exp(opt$minimum)
  fit <- .conditional_gcv(yw, Xw, S, P, b, lam)
  list(
    lambda = lam, g = fit$g, value = fit$value, ed = fit$ed,
    n_eval = n_eval + 1L
  )
}

#' Build cached conditional workspace for core k (Gaussian or weighted).
#' @keywords internal
make_core_workspace <- function(zc, X, P, lambda0, bounds, tol, weight = NULL) {
  if (is.null(weight)) {
    S <- crossprod(X)
    b <- crossprod(X, zc)
    yw <- zc
    Xw <- X
  } else {
    sw <- sqrt(pmax(as.numeric(weight), 0))
    Xw <- X * sw
    yw <- zc * sw
    S <- crossprod(Xw)
    b <- crossprod(Xw, yw)
  }
  list(
    S = S, b = b, P = P, X = X, Xw = Xw, yw = yw,
    lambda0 = lambda0, bounds = bounds, tol = tol
  )
}
