# Minimal GLAM helpers for fixed-λ grid benchmarks (Currie–Durbán–Eilers).
# Scattered-data users should prefer [ttpspline()].

#' @keywords internal
#' @noRd
glam_H <- function(M, A) {
  d <- dim(A)
  Amat <- matrix(A, nrow = d[1])
  array(M %*% Amat, c(nrow(M), d[-1]))
}

#' @keywords internal
#' @noRd
glam_Rotate <- function(A) {
  d <- dim(A)
  if (is.null(d) || length(d) == 1L) return(A)
  aperm(A, c(seq_along(d)[-1L], 1L))
}

#' @keywords internal
#' @noRd
glam_RH <- function(M, A) glam_Rotate(glam_H(M, A))

#' @keywords internal
#' @noRd
glam_linear_predictor <- function(Theta, B_list) {
  out <- Theta
  for (k in seq_along(B_list)) out <- glam_RH(B_list[[k]], out)
  out
}

#' @keywords internal
#' @noRd
glam_Bt_y <- function(Y, B_list) {
  out <- Y
  for (k in seq_along(B_list)) out <- glam_RH(t(B_list[[k]]), out)
  out
}

#' @keywords internal
#' @noRd
glam_penalty <- function(p_vec, lambda, penalty_order = 2L) {
  d <- length(p_vec)
  lambda <- rep(as.numeric(lambda), length.out = d)
  npar <- prod(p_vec)
  P <- matrix(0, npar, npar)
  for (k in seq_len(d)) {
    Dd <- difference_penalty(p_vec[k], penalty_order)
    blocks <- vector("list", d)
    for (j in seq_len(d)) {
      blocks[[j]] <- if (j == k) Dd else diag(p_vec[j])
    }
    K <- blocks[[1L]]
    for (j in 2:d) K <- kronecker(blocks[[j]], K)
    P <- P + lambda[k] * K
  }
  P
}

#' @keywords internal
#' @noRd
glam_gram_unweighted <- function(B_list) {
  G <- crossprod(B_list[[1L]])
  if (length(B_list) == 1L) return(G)
  for (k in 2:length(B_list)) G <- kronecker(crossprod(B_list[[k]]), G)
  G
}

#' Weighted Gram X'WX for d=2 without forming X (Currie).
#' @keywords internal
#' @noRd
glam_xtwx_d2 <- function(B1, B2, W) {
  p1 <- ncol(B1); p2 <- ncol(B2); n2 <- nrow(B2)
  stopifnot(nrow(B1) == nrow(W), ncol(W) == n2)
  XtWX <- matrix(0, p1 * p2, p1 * p2)
  for (v in seq_len(n2)) {
    G1 <- crossprod(B1, B1 * W[, v])
    bv <- B2[v, ]
    XtWX <- XtWX + kronecker(tcrossprod(bv), G1)
  }
  XtWX
}

#' Weighted Gram for d=3 without forming X.
#' @keywords internal
#' @noRd
glam_xtwx_d3 <- function(B1, B2, B3, W) {
  p1 <- ncol(B1); p2 <- ncol(B2); p3 <- ncol(B3)
  n3 <- nrow(B3)
  stopifnot(dim(W)[1] == nrow(B1), dim(W)[2] == nrow(B2), dim(W)[3] == n3)
  XtWX <- matrix(0, p1 * p2 * p3, p1 * p2 * p3)
  for (u in seq_len(n3)) {
    G12 <- glam_xtwx_d2(B1, B2, W[, , u])
    bu <- B3[u, ]
    XtWX <- XtWX + kronecker(tcrossprod(bu), G12)
  }
  XtWX
}

#' Weighted Gram dispatcher (d = 1,2,3).
#' @keywords internal
#' @noRd
glam_xtwx <- function(B_list, W) {
  d <- length(B_list)
  if (d == 1L) {
    return(crossprod(B_list[[1L]], B_list[[1L]] * as.numeric(W)))
  }
  if (d == 2L) return(glam_xtwx_d2(B_list[[1]], B_list[[2]], W))
  if (d == 3L) return(glam_xtwx_d3(B_list[[1]], B_list[[2]], B_list[[3]], W))
  stop("glam_xtwx: only d = 1,2,3 implemented.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.glam_solve <- function(A, b) {
  ridge <- ridge_scale(A, multiplier = 1e-8)
  m <- nrow(A)
  drop(solve_spd(A + ridge * diag(m), b))
}

#' Build marginal B-spline bases on grid axes (Currie / GLAM setup).
#'
#' @param axes List of numeric vectors (one per margin), or a single vector for d=1.
#' @param k Basis size per margin (scalar or length-`d`).
#' @param degree B-spline degree.
#' @return List with `B` (list of bases), `knots`, `axes`, `p`.
#' @export
#' @examples
#' axes <- list(age = 20:50, year = 1990:2010)
#' bb <- glam_grid_bases(axes, k = 8)
#' length(bb$B)
glam_grid_bases <- function(axes, k = 10, degree = 3L) {
  if (is.numeric(axes) && is.null(dim(axes))) axes <- list(axes)
  if (!is.list(axes)) stop("`axes` must be a list of grid coordinates.", call. = FALSE)
  d <- length(axes)
  k <- rep(as.integer(k), length.out = d)
  degree <- as.integer(degree)
  knots <- vector("list", d)
  B <- vector("list", d)
  for (j in seq_len(d)) {
    gj <- as.numeric(axes[[j]])
    knots[[j]] <- make_knots(gj, k = k[j], degree = degree)
    B[[j]] <- bspline_basis(gj, knots[[j]], degree = degree)
  }
  names(B) <- names(axes)
  names(knots) <- names(axes)
  list(B = B, knots = knots, axes = axes, p = k, degree = degree)
}

#' Fixed-λ Gaussian GLAM on a d-way grid (Currie–Durbán–Eilers).
#'
#' Exposed for compression benchmarks; scattered-data users should
#' prefer [ttpspline()].
#'
#' @param Y d-way array.
#' @param B_list Marginal bases (`n_k × p_k`).
#' @param lambda Fixed smoothing (scalar or length-`d`).
#' @param penalty_order Difference penalty order.
#' @export
glam_fit_gaussian <- function(Y, B_list, lambda = 1, penalty_order = 2L) {
  d <- length(B_list)
  stopifnot(length(dim(Y)) == d)
  for (k in seq_len(d)) stopifnot(nrow(B_list[[k]]) == dim(Y)[k])
  p_vec <- vapply(B_list, ncol, integer(1))
  lambda <- rep(as.numeric(lambda), length.out = d)
  intercept <- mean(Y)
  Yc <- Y - intercept
  XtX <- glam_gram_unweighted(B_list)
  P <- glam_penalty(p_vec, lambda, penalty_order)
  coef <- .glam_solve(XtX + P, as.numeric(glam_Bt_y(Yc, B_list)))
  Theta <- array(coef, p_vec)
  mu <- intercept + glam_linear_predictor(Theta, B_list)
  list(
    Theta = Theta, intercept = intercept, mu = mu, eta = mu,
    lambda = lambda, npar = prod(p_vec), p = p_vec,
    family = "gaussian", method = "GLAM-fixed",
    penalty_order = as.integer(penalty_order)
  )
}

#' Fixed-λ Poisson GLAM on a grid (Currie–Durbán–Eilers PIRLS).
#'
#' Fits a Poisson log-linear P-spline on a regular multiway array using
#' GLAM / rotated-H algebra (no explicit Kronecker design). Supports an
#' optional exposure/`offset` array (e.g. mortality person-years).
#' Method of Currie, Durbán and Eilers (2006).
#'
#' @param Y d-way count array.
#' @param B_list Marginal bases (`n_k × p_k`), e.g. from [glam_grid_bases()].
#' @param lambda Fixed smoothing (scalar or length-`d`).
#' @param offset Optional d-way array added on the linear predictor
#'   (typically `log(exposure)`). Default `0`.
#' @param penalty_order Difference penalty order.
#' @param pirls_maxit Maximum PIRLS iterations.
#' @param tol Relative deviance change for stopping.
#' @param trace Print PIRLS progress.
#' @return A list with `Theta`, `intercept`, `eta`, `mu`, `deviance`, etc.
#' @references
#' Currie, I. D., Durban, M. and Eilers, P. H. C. (2006).
#' Generalized linear array models with applications to multidimensional
#' smoothing. *Journal of the Royal Statistical Society: Series B*.
#' @export
#' @examples
#' data(glam_poisson)
#' bb <- glam_grid_bases(list(age = glam_poisson$age, year = glam_poisson$year), k = 8)
#' fit <- glam_fit_poisson(glam_poisson$Y, bb$B, lambda = c(10, 1),
#'                         offset = log(glam_poisson$exposure))
#' fit$deviance
glam_fit_poisson <- function(Y, B_list, lambda = 1, offset = NULL,
                             penalty_order = 2L, pirls_maxit = 25L,
                             tol = 1e-8, trace = FALSE) {
  d <- length(B_list)
  if (is.null(dim(Y)) || length(dim(Y)) != d) {
    stop("Y must be a d-way array matching length(B_list).", call. = FALSE)
  }
  for (k in seq_len(d)) {
    if (nrow(B_list[[k]]) != dim(Y)[k]) {
      stop("nrow(B_list[[", k, "]]) must equal dim(Y)[", k, "].", call. = FALSE)
    }
  }
  if (d > 3L) {
    stop("glam_fit_poisson currently supports d = 1,2,3 (Currie grid).", call. = FALSE)
  }
  p_vec <- vapply(B_list, ncol, integer(1))
  lambda <- rep(as.numeric(lambda), length.out = d)
  P <- glam_penalty(p_vec, lambda, penalty_order)
  fam <- stats::poisson()
  y_vec <- as.numeric(Y)
  if (is.null(offset)) {
    off <- array(0, dim(Y))
  } else {
    off <- as.array(offset)
    if (!identical(as.integer(dim(off)), as.integer(dim(Y)))) {
      stop("`offset` must have the same dimensions as Y.", call. = FALSE)
    }
  }
  off_vec <- as.numeric(off)

  # initialise at log(mean rate) on the offset-adjusted scale
  rate0 <- mean(y_vec / pmax(exp(off_vec), 1e-12))
  intercept <- log(max(rate0, 1e-8))
  Theta <- array(0, p_vec)
  eta <- off + intercept
  mu <- array(exp(pmin(pmax(as.numeric(eta), -20), 20)), dim(Y))
  dev <- glm_deviance(fam, y_vec, as.numeric(mu))
  hist <- list()

  for (it in seq_len(as.integer(pirls_maxit))) {
    work <- glm_working(fam, y_vec, as.numeric(eta))
    W <- array(work$weight, dim(Y))
    z <- array(work$z, dim(Y))
    # η = offset + α + Xθ  ⇒  X'W X θ = X'W (z - offset - α)
    XtWX <- glam_xtwx(B_list, W)
    rhs <- as.numeric(glam_Bt_y(W * (z - off - intercept), B_list))
    coef <- .glam_solve(XtWX + P, rhs)
    Theta <- array(coef, p_vec)
    f <- as.numeric(glam_linear_predictor(Theta, B_list))
    intercept <- sum(work$weight * (work$z - off_vec - f)) /
      max(sum(work$weight), 1e-12)
    eta <- off + intercept + glam_linear_predictor(Theta, B_list)
    mu <- array(exp(pmin(pmax(as.numeric(eta), -20), 20)), dim(Y))
    dev_new <- glm_deviance(fam, y_vec, as.numeric(mu))
    hist[[it]] <- list(pirls = it, deviance = dev_new)
    if (isTRUE(trace)) {
      cat(sprintf("  GLAM-Poisson PIRLS %2d | dev=%.6g\n", it, dev_new))
    }
    rel <- abs(dev - dev_new) / max(1, abs(dev))
    dev <- dev_new
    if (it > 2L && rel < tol) break
  }

  list(
    Theta = Theta,
    intercept = intercept,
    eta = eta,
    mu = mu,
    offset = off,
    lambda = lambda,
    p = p_vec,
    npar = prod(p_vec),
    deviance = dev,
    n_pirls = length(hist),
    history = hist,
    family = "poisson",
    method = "GLAM-PIRLS-fixed",
    penalty_order = as.integer(penalty_order)
  )
}

#' Simulate a Currie–Durbán–Eilers-style Poisson age × year array.
#'
#' Builds a smooth log-rate surface on an age–period grid, applies exposures,
#' and draws Poisson counts — the classical GLAM mortality-style demo setting.
#'
#' @param n_age,n_year Grid sizes.
#' @param seed RNG seed.
#' @return A list with `Y`, `exposure`, `age`, `year`, `eta` (true log mean),
#'   `mu` (true mean), suitable for [glam_fit_poisson()].
#' @export
simulate_glam_poisson <- function(n_age = 41L, n_year = 31L, seed = 44L) {
  set.seed(as.integer(seed))
  age <- seq(20, 60, length.out = as.integer(n_age))
  year <- seq(1980, 2010, length.out = as.integer(n_year))
  A <- (age - mean(age)) / sd(age)
  Tm <- (year - mean(year)) / sd(year)
  eta <- outer(
    0.15 * A + 0.35 * A^2,
    -0.25 * Tm
  ) + outer(
    sin(2 * pi * (age - min(age)) / diff(range(age))),
    0.4 * cos(2 * pi * (year - min(year)) / diff(range(year)))
  )
  eta <- eta - mean(eta) + log(2)
  # exposure grows mildly with year
  exposure <- outer(rep(200, length(age)), 1 + 0.02 * (year - min(year)))
  mu <- exposure * exp(eta)
  Y <- array(stats::rpois(length(mu), mu), dim(mu))
  dimnames(Y) <- list(age = as.character(round(age, 1)),
                      year = as.character(round(year, 1)))
  list(
    Y = Y,
    exposure = exposure,
    age = age,
    year = year,
    eta = eta,
    mu = mu,
    name = "glam_poisson",
    seed = as.integer(seed)
  )
}
