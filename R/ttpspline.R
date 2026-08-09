#' Tensor-Train P-spline for multidimensional smooth / GLM regression.
#'
#' Fits a **non-additive** tensor-product P-spline surface whose coefficient
#' array is represented in Tensor-Train (TT) format. Observations may be
#' arbitrarily scattered in \eqn{[a,b]^d}; no data grid is required.
#'
#' Three orthogonal choices:
#' \itemize{
#'   \item \code{optimizer}: estimation philosophy —
#'     structure-aware \code{ALS} / \code{PIRLS-ALS}, or direct penalized
#'     likelihood \code{GD} / \code{LBFGS} / \code{Adam}.
#'     \code{auto} is a simple family-aware default (Gaussian \(\to\) ALS,
#'     Poisson \(\to\) PIRLS-ALS, binomial \(\to\) LBFGS); always overridable.
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
#' @param optimizer One of:
#'   \itemize{
#'     \item `"auto"` — documented family default:
#'       Gaussian \(\to\) `"ALS"`, Poisson \(\to\) `"PIRLS-ALS"`,
#'       binomial \(\to\) `"LBFGS"` (see `optimizer_requested` /
#'       `optimizer_used` on the fit);
#'     \item `"ALS"` / `"PIRLS-ALS"` — structure-aware ALS / PIRLS+ALS;
#'     \item `"Damped-Newton-ALS"` — conditional Newton + Armijo on true \(Q_k\);
#'     \item `"LBFGS-ALS"` — block/core-wise L-BFGS on each \(Q_k\) (≠ global LBFGS);
#'     \item `"GD"` — global first-order on \(\mathcal L\) (diagnostic);
#'     \item `"LBFGS"` — global quasi-Newton on \(\mathcal L\);
#'     \item `"hybrid"` — experimental ALS→LBFGS polish;
#'     \item `"Adam"` — optional Keras (not yet implemented).
#'   }
#' @param backend `"auto"`, `"R"`, `"Rcpp"`, or `"keras"`. Overridden by
#'   `control$backend` only when this argument is `"auto"` and control is not;
#'   prefer setting backend here or in [tt_control()].
#' @param init Optional TT cores from [tt_initialize()] for fair optimizer
#'   comparisons; `NULL` draws from `control$seed`.
#' @param control A [tt_control()] list.
#' @param monitor If `TRUE`, print iteration progress (sets `control$trace`).
#'   Convenient alias of `tt_control(monitor = TRUE)` / `tt_control(trace = TRUE)`.
#' @param knots Optional list of knot vectors (advanced).
#'
#' @return An object of class `"ttpspline"`.
#'
#' @examples
#' ## Packaged classic surfaces
#' data(ishigami)
#' X <- as.matrix(ishigami[, c("x1", "x2", "x3")])
#' fit <- ttpspline(ishigami$y, X, rank = 2, k = 6, lambda = 1,
#'                  control = tt_control(max_sweeps = 6, compute_edf = FALSE))
#' summary(fit)
#'
#' data(sobol_g)
#' data(friedman)
#'
#' ## On-the-fly simulation (Gaussian / Poisson / Bernoulli)
#' set.seed(1)
#' n <- 400
#' X <- matrix(runif(n * 3), n, 3)
#' f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + X[, 3]
#'
#' ## Gaussian (auto -> ALS)
#' y <- f + rnorm(n, 0, 0.3)
#' fit_g <- ttpspline(y, X, family = gaussian(), rank = 2, k = 6,
#'                    lambda = 1, control = tt_control(max_sweeps = 8))
#' tt_complexity(fit_g)
#'
#' ## Poisson (auto -> PIRLS-ALS)
#' yp <- rpois(n, exp(f - mean(f) + log(2)))
#' fit_p <- ttpspline(yp, X, family = poisson(), rank = 2, k = 6, lambda = 1,
#'                    control = tt_control(pirls_maxit = 15, compute_edf = FALSE))
#' fit_p$optimizer_used
#'
#' ## Bernoulli (auto -> LBFGS)
#' yb <- rbinom(n, 1, plogis(1.2 * (f - mean(f))))
#' fit_b <- ttpspline(yb, X, family = binomial(), rank = 2, k = 6, lambda = 5,
#'                    control = tt_control(lbfgs_maxit = 100, compute_edf = FALSE))
#' fit_b$optimizer_used
#' predict(fit_b, X[1:3, ], type = "response")
#'
#' ## Watch iteration progress
#' \dontrun{
#' fit_m <- ttpspline(ishigami$y, X, rank = 2, k = 6, lambda = 1,
#'                    monitor = TRUE, control = tt_control(max_sweeps = 8))
#' }
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
                      optimizer = c("auto", "ALS", "PIRLS-ALS",
                                    "Damped-Newton-ALS", "LBFGS-ALS",
                                    "GD", "LBFGS", "hybrid", "Adam"),
                      backend = c("auto", "R", "Rcpp", "keras"),
                      init = NULL,
                      control = tt_control(),
                      monitor = FALSE,
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
  if (isTRUE(monitor)) {
    control$trace <- TRUE
    control$monitor <- TRUE
  }
  opt_res <- .resolve_optimizer(match.arg(optimizer), key)
  optimizer_requested <- opt_res$requested
  optimizer_used <- opt_res$used
  optimizer_reason <- opt_res$reason
  optimizer_label <- optimizer_used
  # Internal dispatch token (PIRLS-ALS shares the ALS/PIRLS code path).
  optimizer <- opt_res$dispatch
  backend_arg <- match.arg(backend)
  # Prefer explicit ttpspline(backend=...) over control when not auto
  if (!identical(backend_arg, "auto")) {
    control$backend <- backend_arg
  }
  # Detailed ALS/PIRLS logs live on the R path; with monitor + backend=auto,
  # prefer R so iteration lines actually appear.
  if (isTRUE(control$trace) && identical(backend_arg, "auto") &&
      identical(control$backend, "auto") &&
      optimizer %in% c("ALS", "Damped-Newton-ALS", "LBFGS-ALS", "GD", "hybrid")) {
    control$backend <- "R"
  }
  backend <- resolve_backend(control, optimizer = optimizer)
  ranks <- tt_rank(rank, d = d)
  lambda_spec <- parse_lambda_spec(lambda, d = d, control = control)
  if (isTRUE(control$trace)) {
    lam_lab <- if (identical(lambda_spec$method, "cGCV")) {
      "cGCV"
    } else {
      paste(sprintf("%.3g", lambda_spec$values), collapse = ",")
    }
    cat(sprintf(
      "TTPsplines | family=%s | optimizer=%s | backend=%s | lambda=%s\n",
      fam$family, optimizer_used, backend, lam_lab
    ))
  }

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

  lam_info <- .tt_lambda_boundary_info(
    as.numeric(raw$lambda), raw$method_lambda, control
  )

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
      lambda_bounds = lam_info$lambda_bounds,
      lambda_boundary = lam_info$lambda_boundary,
      lambda_at_boundary = lam_info$lambda_at_boundary,
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
      convergence = raw$convergence %||% list(
        overall = isTRUE(raw$converged),
        pirls = NA,
        als = NA,
        reason = NA_character_
      ),
      optimizer = optimizer_used,
      optimizer_requested = optimizer_requested,
      optimizer_used = optimizer_used,
      optimizer_reason = optimizer_reason,
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

#' Resolve public optimizer choice to a dispatch token + transparency fields.
#'
#' Family-aware `auto` rules (v1; simple and documented):
#' Gaussian → ALS; Poisson → PIRLS-ALS; binomial → LBFGS.
#'
#' @keywords internal
#' @noRd
.resolve_optimizer <- function(optimizer, family_key) {
  requested <- optimizer
  if (!identical(requested, "auto")) {
    dispatch <- if (identical(requested, "PIRLS-ALS")) "ALS" else requested
    return(list(
      requested = requested,
      used = requested,
      reason = "user-specified",
      dispatch = dispatch
    ))
  }
  if (identical(family_key, "gaussian")) {
    used <- "ALS"
    reason <- "gaussian family default"
  } else if (identical(family_key, "poisson")) {
    used <- "PIRLS-ALS"
    reason <- "poisson family default"
  } else if (identical(family_key, "bernoulli")) {
    used <- "LBFGS"
    reason <- "binomial family default"
  } else {
    stop("Unsupported family for optimizer='auto'.", call. = FALSE)
  }
  dispatch <- if (identical(used, "PIRLS-ALS")) "ALS" else used
  list(
    requested = "auto",
    used = used,
    reason = reason,
    dispatch = dispatch
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

  if (identical(optimizer, "hybrid")) {
    return(tt_hybrid_fit(
      y, basis, fam, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores
    ))
  }

  if (identical(optimizer, "LBFGS")) {
    return(tt_lbfgs_fit(
      y, basis, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores,
      family = if (identical(key, "gaussian")) NULL else fam
    ))
  }

  if (identical(optimizer, "GD")) {
    return(tt_gd_fit(
      y, basis, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores,
      family = if (identical(key, "gaussian")) NULL else fam
    ))
  }

  if (identical(optimizer, "Damped-Newton-ALS")) {
    if (identical(key, "gaussian")) {
      # closed-form ALS is the natural Gaussian path; DN reduces to ridge ALS
      warning("Damped-Newton-ALS on Gaussian uses the same path as ALS.",
              call. = FALSE)
      out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                       init_cores = init_cores)
      out$optimizer <- "Damped-Newton-ALS"
      out$backend <- "R"
      return(out)
    }
    return(tt_damped_newton_als_fit(
      y, basis, fam, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores
    ))
  }

  if (identical(optimizer, "LBFGS-ALS")) {
    if (identical(key, "gaussian")) {
      warning("LBFGS-ALS on Gaussian falls back to ALS (closed form).",
              call. = FALSE)
      out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                       init_cores = init_cores)
      out$optimizer <- "LBFGS-ALS"
      out$backend <- "R"
      return(out)
    }
    return(tt_lbfgs_als_fit(
      y, basis, fam, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores
    ))
  }

  # ALS / PIRLS-ALS (default structure-aware path)
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
