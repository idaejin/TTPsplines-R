#' Tensor-Train P-spline for multidimensional smooth / GLM regression.
#'
#' Fits a **non-additive** tensor-product P-spline surface whose coefficient
#' array is represented in Tensor-Train (TT) format. Observations may be
#' arbitrarily scattered in \eqn{[a,b]^d}; no data grid is required.
#'
#' @param y Numeric response (`0/1` for binomial; counts for Poisson).
#' @param X Numeric matrix / data frame of continuous covariates (`n × d`).
#' @param family A [stats::family()] object or one of
#'   `"gaussian"`, `"poisson"`, `"binomial"`.
#' @param rank TT rank (scalar or length `d-1`); see [tt_rank()].
#' @param k Number of B-spline basis functions per margin.
#' @param degree B-spline degree.
#' @param penalty_order Difference penalty order.
#' @param lambda Numeric (isotropic / anisotropic fixed) or `"cGCV"`.
#'   `"cFS"` / `"cREML"` are reserved for Paper 2.
#' @param control A [tt_control()] list.
#' @param knots Optional list of knot vectors (advanced).
#'
#' @return An object of class `"ttpspline"`.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' X <- matrix(runif(n * 3), n, 3)
#' f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + X[, 3]
#' y <- f + rnorm(n, 0, 0.3)
#' fit <- ttpspline(y, X, family = gaussian(), rank = 2, k = 6,
#'                  lambda = 1, control = tt_control(max_sweeps = 8))
#' tt_complexity(fit)
#'
#' @export
ttpspline <- function(y,
                      X,
                      family = stats::gaussian(),
                      rank = 3,
                      k = 10,
                      degree = 3,
                      penalty_order = 2,
                      lambda = "cGCV",
                      control = tt_control(),
                      knots = NULL) {
  cl <- match.call()
  fam <- normalize_family(family)
  y <- as.numeric(y)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (nrow(X) != length(y)) stop("nrow(X) must equal length(y).", call. = FALSE)
  if (anyNA(y) || anyNA(X)) stop("NA values are not supported in v0.", call. = FALSE)
  d <- ncol(X)
  if (d < 2L) stop("Need at least d = 2 covariates.", call. = FALSE)

  key <- family_key(fam)
  if (identical(key, "bernoulli")) {
    if (!all(y %in% c(0, 1))) stop("binomial/Bernoulli requires y in {0,1}.", call. = FALSE)
  }
  if (identical(key, "poisson") && any(y < 0)) {
    stop("poisson requires non-negative y.", call. = FALSE)
  }

  if (!inherits(control, "tt_control")) {
    control <- do.call(tt_control, as.list(control))
  }
  backend <- resolve_backend(control)
  ranks <- tt_rank(rank, d = d)
  lambda_spec <- parse_lambda_spec(lambda, d = d)

  bs <- build_marginal_bases(X, k = k, degree = degree, knots = knots)
  basis <- bs$basis
  p <- bs$k

  if (identical(key, "gaussian")) {
    raw <- if (identical(backend, "Rcpp")) {
      tt_als_fit_rcpp(y, basis, ranks, lambda_spec, control, penalty_order)
    } else {
      out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order)
      out$backend <- "R"
      out
    }
  } else {
    raw <- if (identical(backend, "Rcpp")) {
      tt_pirls_fit_rcpp(y, basis, fam, ranks, lambda_spec, control, penalty_order)
    } else {
      tt_pirls_fit(y, basis, fam, ranks, lambda_spec, control, penalty_order)
    }
  }

  npar_tt <- tt_npar(p, ranks)
  npar_full <- dense_npar(p, d)
  residuals_resp <- y - raw$mu

  structure(
    list(
      call = cl,
      family = fam,
      family_key = key,
      y = y,
      d = d,
      n = length(y),
      k = p,
      degree = bs$degree,
      knots = bs$knots,
      penalty_order = as.integer(penalty_order),
      cores = raw$cores,
      rank = ranks,
      lambda = as.numeric(raw$lambda),
      lambda_method = raw$method_lambda,
      intercept = raw$intercept,
      fitted.values = raw$mu,
      linear.predictors = raw$eta,
      residuals = residuals_resp,
      deviance = raw$deviance,
      edf = NA_real_, # joint EDF optional later
      npar_tt = npar_tt,
      npar_dense = npar_full,
      compression_ratio = npar_full / max(npar_tt, 1),
      converged = isTRUE(raw$converged),
      n_sweeps = raw$n_sweeps,
      n_pirls = raw$n_pirls,
      n_criterion_evals = raw$n_criterion_evals,
      history = raw$history,
      backend = raw$backend %||% backend,
      sparse_backend = control$sparse,
      timing = raw$elapsed,
      control = control,
      x_names = colnames(X),
      x_range = apply(X, 2, range)
    ),
    class = "ttpspline"
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
