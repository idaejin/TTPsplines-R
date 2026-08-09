#' Open knot sequence for a univariate B-spline basis.
#' @keywords internal
make_knots <- function(x, k, degree = 3) {
  stopifnot(is.numeric(x), length(x) > 0, k > degree)
  xr <- range(x, na.rm = TRUE)
  inner <- seq(xr[1], xr[2], length.out = k - degree + 1)
  inner <- inner[-c(1, length(inner))]
  c(rep(xr[1], degree + 1), inner, rep(xr[2], degree + 1))
}

#' Evaluate a B-spline basis matrix.
#' @keywords internal
bspline_basis <- function(x, knots, degree = 3) {
  splines::splineDesign(knots, x, ord = degree + 1, outer.ok = TRUE)
}

#' Build marginal bases for columns of X (scattered observations).
#'
#' The design is observation-wise: each row of X gets its own basis
#' evaluations. No Cartesian grid is required.
#'
#' @keywords internal
build_marginal_bases <- function(X, k = 10, degree = 3, knots = NULL) {
  X <- as.matrix(X)
  d <- ncol(X)
  k <- rep(as.integer(k), length.out = d)
  if (is.null(knots)) {
    knots <- lapply(seq_len(d), function(j) {
      make_knots(X[, j], k = k[j], degree = degree)
    })
  }
  stopifnot(length(knots) == d)
  basis <- lapply(seq_len(d), function(j) {
    bspline_basis(X[, j], knots[[j]], degree = degree)
  })
  # Equalize column counts if k was scalar (already); if unequal p, ALS needs equal p
  p <- vapply(basis, ncol, integer(1))
  if (length(unique(p)) != 1L) {
    stop("All margins must have the same number of basis functions (k) in v0.")
  }
  list(
    basis = basis,
    knots = knots,
    degree = as.integer(degree),
    k = as.integer(p[1]),
    d = d,
    n = nrow(X)
  )
}

#' Evaluate bases on newdata using training knots.
#' @keywords internal
eval_marginal_bases <- function(Xnew, knots, degree) {
  Xnew <- as.matrix(Xnew)
  stopifnot(ncol(Xnew) == length(knots))
  lapply(seq_len(ncol(Xnew)), function(j) {
    bspline_basis(Xnew[, j], knots[[j]], degree = degree)
  })
}
