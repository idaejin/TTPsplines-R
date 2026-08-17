#' Default spline domain for a margin (open or cyclic).
#'
#' If all finite values lie in the unit interval 0-1 (within a tiny tolerance), use the
#' unit interval so prediction on a nice 0-1 grid does not fall outside
#' the knot span (which would zero the B-spline basis under `outer.ok`).
#' Otherwise use the empirical range.
#' @keywords internal
#' @noRd
default_basis_domain <- function(x) {
  xr <- range(x, na.rm = TRUE)
  if (xr[1] >= -1e-9 && xr[2] <= 1 + 1e-9) {
    return(c(0, 1))
  }
  xr
}

#' Open knot sequence for a univariate B-spline basis.
#' @keywords internal
make_knots <- function(x, k, degree = 3, xl = NULL, xr = NULL) {
  stopifnot(is.numeric(x), length(x) > 0, k > degree)
  if (is.null(xl) || is.null(xr)) {
    dom <- default_basis_domain(x)
    if (is.null(xl)) xl <- dom[1]
    if (is.null(xr)) xr <- dom[2]
  }
  if (!(is.finite(xl) && is.finite(xr) && xr > xl)) {
    stop("make_knots needs finite xl < xr.", call. = FALSE)
  }
  inner <- seq(xl, xr, length.out = k - degree + 1)
  inner <- inner[-c(1, length(inner))]
  c(rep(xl, degree + 1), inner, rep(xr, degree + 1))
}

#' Cyclic B-spline basis on period (xl, xr) with `k` periodic functions (Eilers-style).
#'
#' Equally spaced extended knots; wrap the last `degree` columns into the first
#' so the basis is periodic with period `xr - xl`.
#' @keywords internal
#' @noRd
cyclic_bspline_basis <- function(x, knots = NULL, degree = 3,
                                 k = NULL, xl = NULL, xr = NULL) {
  degree <- as.integer(degree)
  if (!is.null(knots)) {
    if (is.null(xl)) {
      xl <- attr(knots, "xl")
      if (is.null(xl)) xl <- min(as.numeric(knots))
    }
    if (is.null(xr)) {
      xr <- attr(knots, "xr")
      if (is.null(xr)) xr <- max(as.numeric(knots))
    }
    if (is.null(k)) {
      k <- attr(knots, "k")
      if (is.null(k)) k <- length(as.numeric(knots)) - 2L * degree - 1L
    }
  }
  if (is.null(xl) || is.null(xr) || is.null(k)) {
    stop("cyclic_bspline_basis needs knots with xl/xr/k or explicit args.",
         call. = FALSE)
  }
  k <- as.integer(k)
  if (k < degree + 1L) stop("cyclic basis needs k >= degree + 1.", call. = FALSE)
  if (any(!is.finite(x))) stop("non-finite x in cyclic basis.", call. = FALSE)
  if (any(x < xl - 1e-10 | x > xr + 1e-10)) {
    stop("x outside cyclic period [", xl, ", ", xr, "].", call. = FALSE)
  }
  x <- pmin(pmax(as.numeric(x), xl), xr)
  dx <- (xr - xl) / k
  knots_ext <- seq(xl - degree * dx, xr + degree * dx, by = dx)
  B <- splines::splineDesign(knots_ext, x, ord = degree + 1L, outer.ok = TRUE)
  # ncol(B) = k + degree; wrap last `degree` into first `degree`
  if (degree > 0L) {
    B[, seq_len(degree)] <- B[, seq_len(degree), drop = FALSE] +
      B[, k + seq_len(degree), drop = FALSE]
  }
  B[, seq_len(k), drop = FALSE]
}

#' Evaluate a (non-cyclic) B-spline basis matrix.
#' @keywords internal
bspline_basis <- function(x, knots, degree = 3) {
  splines::splineDesign(knots, x, ord = degree + 1, outer.ok = TRUE)
}

#' Equally spaced knot metadata for a cyclic margin (period xl to xr).
#' @keywords internal
#' @noRd
make_cyclic_knots <- function(xl, xr, k, degree = 3) {
  k <- as.integer(k)
  degree <- as.integer(degree)
  if (!(is.finite(xl) && is.finite(xr) && xr > xl)) {
    stop("cyclic basis requires finite xl < xr.", call. = FALSE)
  }
  if (k < degree + 1L) {
    stop("cyclic basis needs k >= degree + 1.", call. = FALSE)
  }
  # Store period endpoints + k; evaluation rebuilds the extended knot sequence.
  structure(
    c(xl, xr),
    cyclic = TRUE,
    xl = xl,
    xr = xr,
    k = k,
    degree = degree
  )
}

#' Normalize cyclic margin flags to length d.
#' @keywords internal
#' @noRd
normalize_cyclic <- function(cyclic, d) {
  d <- as.integer(d)
  if (is.null(cyclic)) return(rep(FALSE, d))
  cyclic <- as.logical(cyclic)
  if (length(cyclic) == 1L) cyclic <- rep(cyclic, d)
  if (length(cyclic) != d) {
    stop("`cyclic` must be length 1 or d = ", d, ".", call. = FALSE)
  }
  if (anyNA(cyclic)) stop("`cyclic` contains NA.", call. = FALSE)
  cyclic
}

#' Period (xl, xr) for a cyclic margin.
#' @keywords internal
#' @noRd
cyclic_period_range <- function(x, period = NULL) {
  if (!is.null(period)) {
    period <- as.numeric(period)
    if (length(period) != 2L || !(period[2] > period[1])) {
      stop("`period` must be c(xl, xr) with xr > xl.", call. = FALSE)
    }
    return(period)
  }
  dom <- default_basis_domain(x)
  if (diff(dom) < .Machine$double.eps) {
    stop("cyclic margin has zero range; pass period = c(xl, xr).", call. = FALSE)
  }
  dom
}

#' Build marginal bases for columns of X (scattered observations).
#'
#' The design is observation-wise: each row of X gets its own basis
#' evaluations. No Cartesian grid is required.
#'
#' @param cyclic Logical scalar or length-`d` flags: use a **circular**
#'   B-spline basis (and, at penalty construction, a circular difference
#'   penalty) on that margin. For hour-of-day, map to the unit interval 0-1
#'   (fraction of day) and set the corresponding flag to `TRUE` (period
#'   defaults to 0-1).
#' @param period Optional length-`d` list of `c(xl,xr)` periods for cyclic
#'   margins (`NULL` entries use `cyclic_period_range()`).
#' @keywords internal
build_marginal_bases <- function(X, k = 10, degree = 3, knots = NULL,
                                 cyclic = NULL, period = NULL) {
  X <- as.matrix(X)
  d <- ncol(X)
  k <- rep(as.integer(k), length.out = d)
  cyclic <- normalize_cyclic(cyclic, d)
  if (!is.null(period) && !is.list(period)) {
    stop("`period` must be a list of length d (or NULL).", call. = FALSE)
  }
  if (is.null(period)) period <- vector("list", d)
  if (length(period) != d) {
    stop("`period` must have length d = ", d, ".", call. = FALSE)
  }

  if (is.null(knots)) {
    knots <- vector("list", d)
    for (j in seq_len(d)) {
      if (cyclic[j]) {
        pr <- cyclic_period_range(X[, j], period[[j]])
        knots[[j]] <- structure(
          make_cyclic_knots(pr[1], pr[2], k = k[j], degree = degree),
          cyclic = TRUE,
          xl = pr[1],
          xr = pr[2]
        )
      } else {
        knots[[j]] <- structure(
          make_knots(X[, j], k = k[j], degree = degree),
          cyclic = FALSE
        )
      }
    }
  } else {
    stopifnot(length(knots) == d)
    for (j in seq_len(d)) {
      if (is.null(attr(knots[[j]], "cyclic"))) {
        attr(knots[[j]], "cyclic") <- cyclic[j]
      }
      if (isTRUE(attr(knots[[j]], "cyclic"))) {
        if (is.null(attr(knots[[j]], "xl"))) {
          pr <- cyclic_period_range(X[, j], period[[j]])
          attr(knots[[j]], "xl") <- pr[1]
          attr(knots[[j]], "xr") <- pr[2]
        }
      }
      cyclic[j] <- isTRUE(attr(knots[[j]], "cyclic"))
    }
  }

  basis <- lapply(seq_len(d), function(j) {
    if (cyclic[j]) {
      cyclic_bspline_basis(X[, j], knots[[j]], degree = degree)
    } else {
      bspline_basis(X[, j], knots[[j]], degree = degree)
    }
  })
  p <- vapply(basis, ncol, integer(1))
  if (length(unique(p)) != 1L) {
    stop("All margins must have the same number of basis functions (k) in v0.")
  }
  structure(
    list(
      basis = `attr<-`(basis, "cyclic", cyclic),
      knots = knots,
      degree = as.integer(degree),
      k = as.integer(p[1]),
      d = d,
      n = nrow(X),
      cyclic = cyclic
    ),
    class = "tt_marginal_bases"
  )
}

#' Evaluate bases on newdata using training knots.
#' @keywords internal
eval_marginal_bases <- function(Xnew, knots, degree, cyclic = NULL) {
  Xnew <- as.matrix(Xnew)
  stopifnot(ncol(Xnew) == length(knots))
  d <- length(knots)
  if (is.null(cyclic)) {
    cyclic <- vapply(knots, function(kn) isTRUE(attr(kn, "cyclic")), logical(1))
  } else {
    cyclic <- normalize_cyclic(cyclic, d)
  }
  for (j in seq_len(d)) {
    if (cyclic[j]) next
    kn <- as.numeric(knots[[j]])
    xl <- min(kn)
    xr <- max(kn)
    xj <- Xnew[, j]
    if (any(xj < xl - 1e-10 | xj > xr + 1e-10, na.rm = TRUE)) {
      warning(
        sprintf(
          paste0(
            "newdata margin %d has values outside knot span [%.4g, %.4g]; ",
            "open B-splines evaluate to 0 there (prediction collapses toward ",
            "the intercept). Prefer knots covering the prediction domain ",
            "(unit-interval covariates now default to [0,1])."
          ),
          j, xl, xr
        ),
        call. = FALSE
      )
    }
  }
  lapply(seq_len(d), function(j) {
    if (cyclic[j]) {
      cyclic_bspline_basis(Xnew[, j], knots[[j]], degree = degree)
    } else {
      bspline_basis(Xnew[, j], knots[[j]], degree = degree)
    }
  })
}
