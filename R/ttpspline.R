#' Tensor-Train P-spline for multidimensional smooth / GLM regression.
#'
#' Fits a **non-additive** tensor-product P-spline surface whose coefficient
#' array is represented in Tensor-Train (TT) format. Observations may be
#' arbitrarily scattered in \eqn{[a,b]^d}; no data grid is required.
#'
#' Three orthogonal choices:
#' \itemize{
#'   \item \code{optimizer}: how TT cores are estimated (`ALS`, `LBFGS`, `Adam`)
#'   \item \code{lambda}: fixed isotropic/anisotropic or automatic `"cGCV"`
#'   \item \code{backend}: computational engine (`auto` / `R` / `Rcpp` / `keras`)
#' }
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
#'   `"cFS"` / `"cREML"` are not implemented yet.
#' @param optimizer `"ALS"` (default), `"LBFGS"`, or `"Adam"` (optional Keras;
#'   not yet implemented — use ALS or LBFGS).
#' @param backend `"auto"`, `"R"`, `"Rcpp"`, or `"keras"`. Overridden by
#'   `control$backend` only when this argument is `"auto"` and control is not;
#'   prefer setting backend here or in [tt_control()].
#' @param init Optional TT cores from [tt_initialize()] for fair optimizer
#'   comparisons; `NULL` draws from `control$seed`.
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
                      optimizer = c("ALS", "LBFGS", "Adam"),
                      backend = c("auto", "R", "Rcpp", "keras"),
                      init = NULL,
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
  optimizer <- match.arg(optimizer)
  backend_arg <- match.arg(backend)
  # Prefer explicit ttpspline(backend=...) over control when not auto
  if (!identical(backend_arg, "auto")) {
    control$backend <- backend_arg
  }
  backend <- resolve_backend(control, optimizer = optimizer)
  ranks <- tt_rank(rank, d = d)
  lambda_spec <- parse_lambda_spec(lambda, d = d, control = control)

  if (!is.null(init)) {
    if (!is.list(init) || length(init) != d) {
      stop("`init` must be a length-d list of TT cores (see tt_initialize()).",
           call. = FALSE)
    }
    for (kk in seq_len(d)) {
      dm <- dim(init[[kk]])
      if (is.null(dm) || length(dm) != 3L ||
          dm[1] != ranks[kk] || dm[3] != ranks[kk + 1L]) {
        stop("init core ", kk, " has incompatible TT dimensions.", call. = FALSE)
      }
    }
  }

  bs <- build_marginal_bases(X, k = k, degree = degree, knots = knots)
  basis <- bs$basis
  p <- bs$k
  if (!is.null(init) && ncol(basis[[1]]) != dim(init[[1]])[2]) {
    stop("init cores k does not match basis size from argument k/knots.",
         call. = FALSE)
  }

  raw <- .ttpspline_dispatch(
    y = y,
    basis = basis,
    fam = fam,
    key = key,
    ranks = ranks,
    lambda_spec = lambda_spec,
    control = control,
    penalty_order = penalty_order,
    optimizer = optimizer,
    backend = backend,
    init_cores = init
  )

  npar_tt <- tt_npar(p, ranks)
  npar_full <- dense_npar(p, d)
  npar_tt_intrinsic <- npar_tt - tt_gauge_dim(ranks)
  residuals_resp <- y - raw$mu

  edf <- NA_real_
  edf_note <- "not computed"
  if (isTRUE(control$compute_edf) && !is.null(raw$penalties)) {
    w_edf <- tt_edf_weights(raw, key, y)
    edf <- tt_joint_edf(
      raw$cores, basis, raw$penalties, as.numeric(raw$lambda),
      weight = w_edf, max_npar = control$edf_max_npar
    )
    if (is.finite(edf)) {
      edf_note <- "joint linearized TT map at convergence"
    } else if (npar_tt > control$edf_max_npar) {
      edf_note <- sprintf("skipped (npar_TT=%d > edf_max_npar=%d)",
                          npar_tt, control$edf_max_npar)
    } else {
      edf_note <- "failed or skipped (see control$compute_edf / edf_max_npar)"
    }
  } else if (!isTRUE(control$compute_edf)) {
    edf_note <- "disabled (control$compute_edf = FALSE)"
  }

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
      rank_internal = ranks[-c(1L, length(ranks))],
      rank_max = max(ranks),
      lambda = as.numeric(raw$lambda),
      lambda_method = raw$method_lambda,
      intercept = raw$intercept,
      fitted.values = raw$mu,
      linear.predictors = raw$eta,
      residuals = residuals_resp,
      deviance = raw$deviance,
      edf = edf,
      edf_note = edf_note,
      npar_tt = npar_tt,
      npar_tt_intrinsic = npar_tt_intrinsic,
      npar_dense = npar_full,
      compression_ratio = npar_full / max(npar_tt, 1),
      converged = isTRUE(raw$converged),
      optimizer = raw$optimizer %||% optimizer,
      n_sweeps = raw$n_sweeps,
      n_pirls = raw$n_pirls,
      n_opt_iter = raw$n_opt_iter %||% NA_integer_,
      n_outer = raw$n_outer %||% NA_integer_,
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

#' @keywords internal
.ttpspline_dispatch <- function(y, basis, fam, key, ranks, lambda_spec,
                                control, penalty_order, optimizer, backend,
                                init_cores) {
  if (identical(optimizer, "Adam")) {
    return(tt_adam_fit(
      y, basis, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores, family = if (identical(key, "gaussian")) NULL else fam
    ))
  }

  if (identical(optimizer, "LBFGS")) {
    return(tt_lbfgs_fit(
      y, basis, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores,
      family = if (identical(key, "gaussian")) NULL else fam
    ))
  }

  # ALS (default)
  if (identical(key, "gaussian")) {
    if (identical(backend, "Rcpp")) {
      tt_als_fit_rcpp(y, basis, ranks, lambda_spec, control, penalty_order,
                      init_cores = init_cores)
    } else {
      out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                       init_cores = init_cores)
      out$backend <- "R"
      out
    }
  } else {
    if (identical(backend, "Rcpp")) {
      tt_pirls_fit_rcpp(y, basis, fam, ranks, lambda_spec, control, penalty_order,
                        init_cores = init_cores)
    } else {
      tt_pirls_fit(y, basis, fam, ranks, lambda_spec, control, penalty_order,
                  init_cores = init_cores)
    }
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
