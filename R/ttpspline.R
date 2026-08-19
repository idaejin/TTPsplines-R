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
#'   \item \code{backend}: `"R"` for ALS/PIRLS sweeps; `"Rcpp"` = kernel helpers
#'     only (not a full C++ ALS fitter); `"keras"` reserved for Adam
#' }
#'
#' @param y Numeric response (`0/1` for binomial; counts for Poisson).
#' @param X Numeric matrix / data frame of continuous covariates (`n × d`).
#' @param family A [stats::family()] object or one of
#'   `"gaussian"`, `"poisson"`, `"binomial"`.
#' @param rank TT rank (scalar or length `d-1`); see [tt_rank()].
#'   **Not** selected automatically — `rank = 2` means exactly rank 2.
#'   For data-driven choice use [tt_rank_select()] then [tt_rank_refit()].
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
#'   prefer setting backend here or in [tt_control()]. ALS / PIRLS always use
#'   the R sweep under the global penalty; `"Rcpp"` accelerates kernels only
#'   (not a full C++ ALS path). See [tt_control()].
#' @param init Optional TT cores from [tt_initialize()] for fair optimizer
#'   comparisons; `NULL` draws from `control$seed`.
#' @param control A [tt_control()] list.
#' @param monitor If `TRUE`, print iteration progress (sets `control$trace`).
#'   Convenient alias of `tt_control(monitor = TRUE)` / `tt_control(trace = TRUE)`.
#' @param knots Optional list of knot vectors (advanced).
#' @param offset Optional numeric vector of length `n` (or scalar) added to the
#'   linear predictor, e.g. `log(exposure)` for Poisson. Default `NULL` (zeros).
#' @param linear Optional unpenalized parametric design matrix (`n × p`), e.g.
#'   `model.matrix(~ 0 + dow, data)`. Do **not** include an intercept column
#'   (`ttps` already estimates one). Estimated jointly with the TT surface via
#'   block-coordinate updates (ALS / PIRLS-ALS).
#' @param smooth Optional additive 1D smooths (mgcv-like). A named list of
#'   specs, e.g.
#'   `list(time = list(x = d$time, bs = "ps", k = 20, m = 2))` or
#'   `list(hour = list(x = hour, bs = "cc", k = 12, m = 2, period = c(0, 24)))`.
#'   `bs`: `"ps"` (open P-spline) or `"cc"` (circular). `m` is the difference
#'   penalty order. Per-term smoothing: `lambda` (numeric or `"cGCV"`, default
#'   from `lambda_smooth`) **or** `target_edf` (choose \(\lambda\) so
#'   \(\mathrm{edf}(\lambda)\approx\) target; requires `m < target_edf <= k`).
#'   Unsupported for LBFGS / GD / hybrid / Adam / LBFGS-ALS / DN-ALS (error).
#' @param lambda_smooth Default smoother penalty for terms in `smooth` that
#'   omit `lambda` and `target_edf` (`"cGCV"` or a nonnegative scalar).
#' @param null_space How to treat the discrete P-spline penalty null space:
#'   \itemize{
#'     \item `"joint"` (default) — single TT on the full coefficient array
#'       (current behaviour; small uniform \(r\) may truncate \(\mathcal N(S)\)).
#'     \item `"profiled"` — **experimental** Gaussian null-space–preserving
#'       fit: profile \(\beta_0\) exactly and ALS-update cores on
#'       \(Q_0 y\) and \(Q_0 Z_k\) (empirical orthogonalization). Not yet
#'       available for GLM / cyclic / `linear` / `smooth`. Requires
#'       \(q^d \le\) `control$null_space_max_npar` (default 4096).
#'   }
#' @param weights Optional non-negative observation weights of length `n`
#'   (or scalar). `NULL` means all ones. Zero weights exclude observations from
#'   the likelihood (useful for masking empty array cells).
#' @param cyclic Logical scalar or length-`d` flags for **circular** P-spline
#'   margins (periodic B-spline basis + circular difference penalty). For
#'   hour-of-day, map to the unit interval 0-1 and set that margin to `TRUE`.
#' @param period Optional length-`d` list of `c(xl,xr)` for cyclic margins.
#'
#' @return An object of class `"ttpspline"`.
#'
#' @examples
#' ## Packaged classic surfaces
#' data(ishigami)
#' X <- as.matrix(ishigami[, c("x1", "x2", "x3")])
#' fit <- ttps(ishigami$y, X, rank = 2, k = 6, lambda = 1,
#'             control = tt_control(max_sweeps = 6, compute_edf = FALSE))
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
#' fit_g <- ttps(y, X, family = gaussian(), rank = 2, k = 6,
#'               lambda = 1, control = tt_control(max_sweeps = 8))
#' tt_complexity(fit_g)
#'
#' ## Poisson (auto -> PIRLS-ALS)
#' yp <- rpois(n, exp(f - mean(f) + log(2)))
#' fit_p <- ttps(yp, X, family = poisson(), rank = 2, k = 6, lambda = 1,
#'               control = tt_control(pirls_maxit = 15, compute_edf = FALSE))
#' fit_p$optimizer_used
#'
#' ## Bernoulli (auto -> LBFGS)
#' yb <- rbinom(n, 1, plogis(1.2 * (f - mean(f))))
#' fit_b <- ttps(yb, X, family = binomial(), rank = 2, k = 6, lambda = 5,
#'               control = tt_control(lbfgs_maxit = 100, compute_edf = FALSE))
#' fit_b$optimizer_used
#' predict(fit_b, X[1:3, ], type = "response")
#'
#' ## Watch iteration progress
#' \dontrun{
#' fit_m <- ttps(ishigami$y, X, rank = 2, k = 6, lambda = 1,
#'               monitor = TRUE, control = tt_control(max_sweeps = 8))
#' }
#'
#' @rdname ttps
#' @export
ttps <- function(y,
                 X = NULL,
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
                      knots = NULL,
                      offset = NULL,
                      linear = NULL,
                      smooth = NULL,
                      lambda_smooth = "cGCV",
                      weights = NULL,
                      cyclic = NULL,
                      period = NULL,
                      null_space = c("joint", "profiled"),
                      array = FALSE,
                      axes = NULL) {
  cl <- match.call()
  null_space <- match.arg(null_space)

  # ------------------------------------------------------------------
  # Array mode: y is a d-way array Y, X is replaced by axes list.
  # ------------------------------------------------------------------
  array_data_out <- NULL   # populated below when array = TRUE
  if (isTRUE(array)) {
    if (!is.array(y) || length(dim(y)) < 2L) {
      stop(
        "`array = TRUE` requires `y` to be a d-way array (d >= 2). ",
        "Provide Y as array(responses, dim = c(n1, ..., nd)).",
        call. = FALSE
      )
    }
    n_grid <- dim(y)
    d_arr  <- length(n_grid)
    if (!is.null(linear) || !is.null(smooth)) {
      stop("`array = TRUE` does not support `linear=` or `smooth=`.", call. = FALSE)
    }
    if (!identical(null_space, "joint")) {
      stop("`array = TRUE` only supports null_space = 'joint'.", call. = FALSE)
    }
    fam_tmp <- normalize_family(family)
    if (!identical(family_key(fam_tmp), "gaussian")) {
      stop("`array = TRUE` is currently Gaussian-only.", call. = FALSE)
    }
    if (is.null(axes)) {
      # Default: unit-interval grid for each margin
      axes <- lapply(n_grid, function(nk) seq(0, 1, length.out = nk))
    } else {
      if (!is.list(axes) || length(axes) != d_arr) {
        stop("`axes` must be a list of length d (one coordinate vector per margin).",
             call. = FALSE)
      }
      for (j in seq_len(d_arr)) {
        if (length(axes[[j]]) != n_grid[j]) {
          stop(sprintf("axes[[%d]] has length %d but dim(y)[%d] = %d.",
                       j, length(axes[[j]]), j, n_grid[j]),
               call. = FALSE)
        }
      }
    }
    if (is.null(names(axes))) names(axes) <- paste0("x", seq_len(d_arr))
    # Build scattered representation: expand.grid (dim1 fastest)
    idx <- expand.grid(lapply(n_grid, seq_len), KEEP.OUT.ATTRS = FALSE)
    X <- do.call(cbind, lapply(seq_len(d_arr), function(j) axes[[j]][idx[[j]]]))
    colnames(X) <- names(axes)
    y_sc <- as.numeric(y)  # dim1 fastest — matches expand.grid row order
    # Store array_data: will be attached after basis is built
    array_data_out <- list(
      Y        = y,       # original d-way array (centred later)
      n_grid   = n_grid,
      axes     = axes,
      k        = NA_integer_,   # updated per core visit in ALS
      L_all    = NULL,          # filled per sweep in ALS
      R_all    = NULL
    )
    y <- y_sc
  }

  fam <- normalize_family(family)
  if (is.null(X)) stop("`X` must be provided (or use `array = TRUE` with `y` as array).", call. = FALSE)
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
  offset <- normalize_offset(offset, length(y))
  weights <- normalize_weights(weights, length(y))
  linear <- normalize_linear(linear, length(y))
  smooth <- normalize_smooth(smooth, length(y), lambda_smooth = lambda_smooth)

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
  optimizer <- opt_res$dispatch
  # linear=/smooth= force structure-aware ALS/PIRLS (incl. binomial auto→LBFGS)
  if (!is.null(linear) || !is.null(smooth)) {
    unsupported <- c("LBFGS", "GD", "hybrid", "Adam",
                     "Damped-Newton-ALS", "LBFGS-ALS")
    if (identical(optimizer_requested, "auto") &&
        identical(optimizer_used, "LBFGS")) {
      optimizer_used <- "PIRLS-ALS"
      optimizer_reason <- paste(
        "linear=/smooth= uses PIRLS-ALS for binomial",
        "(LBFGS path not extended yet)"
      )
      optimizer <- "ALS"
    } else if (optimizer_used %in% unsupported) {
      stop(
        "`linear=` / `smooth=` are currently supported only with optimizer in ",
        "{'auto','ALS','PIRLS-ALS'} (structure-aware path).",
        call. = FALSE
      )
    }
  }
  optimizer_label <- optimizer_used
  backend_arg <- match.arg(backend)
  # Prefer explicit ttps(backend=...) over control when not auto
  if (!identical(backend_arg, "auto")) {
    control$backend <- backend_arg
  }
  # linear=/smooth= are implemented on the R ALS/PIRLS path
  if ((!is.null(linear) || !is.null(smooth)) &&
      (identical(backend_arg, "Rcpp") || identical(control$backend, "Rcpp"))) {
    warning("`linear=` / `smooth=` use the R backend; ignoring backend='Rcpp'.",
            call. = FALSE)
    control$backend <- "R"
    backend_arg <- "R"
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

  bs <- build_marginal_bases(X, k = k, degree = degree, knots = knots,
                             cyclic = cyclic, period = period)
  basis <- bs$basis
  p <- bs$k
  cyclic <- bs$cyclic
  if (!is.null(init) && ncol(basis[[1]]) != dim(init[[1]])[2]) {
    stop("init cores k does not match basis size from argument k/knots.",
         call. = FALSE)
  }

  # Optional profiled null-space (Gaussian Q0-orthogonalized ALS).
  null_info <- NULL
  max_ns <- control$null_space_max_npar %||% 4096L

  if (identical(null_space, "profiled")) {
    if (!identical(key, "gaussian")) {
      stop("null_space = 'profiled' is Gaussian-only for now.", call. = FALSE)
    }
    if (!is.null(linear) || !is.null(smooth)) {
      stop("null_space = 'profiled' does not support linear= / smooth= yet.",
           call. = FALSE)
    }
    if (!optimizer %in% c("auto", "ALS", "Damped-Newton-ALS", "LBFGS-ALS")) {
      stop(
        "null_space = 'profiled' requires ALS (got optimizer='", optimizer, "').",
        call. = FALSE
      )
    }
    if (!identical(backend, "R") && !identical(backend, "auto")) {
      warning("null_space = 'profiled' uses the R ALS backend.", call. = FALSE)
    }
    if (isTRUE(control$trace)) {
      cat(sprintf(
        "TTPsplines | null_space=profiled | q=%d | max_npar=%d\n",
        as.integer(penalty_order), as.integer(max_ns)
      ))
    }
    raw <- tt_als_fit_profiled_null(
      y, basis, ranks, lambda_spec, control,
      penalty_order = penalty_order, init_cores = init,
      offset = offset, weights = weights, max_npar = max_ns
    )
    null_info <- raw$null_space_info
  } else {
    # In array mode, attach marginal bases (B_list) to array_data now that basis is built.
    if (!is.null(array_data_out)) {
      # basis[[k]] has nrow = n_total (scattered rows); marginal basis is
      # bb$B[[k]] of size n_k x p. Reconstruct via glam_grid_bases.
      bb_arr <- glam_grid_bases(array_data_out$axes,
                                k = p, degree = as.integer(degree))
      array_data_out$B_marginal <- bb_arr$B
      # Y_centered will be updated inside ALS after intercept is known.
      # Pass Y (not centred) as Y_centered placeholder; ALS will subtract intercept.
      array_data_out$Y_centered <- array_data_out$Y
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
      init_cores = init,
      offset = offset,
      weights = weights,
      linear = linear,
      smooth = smooth,
      array_data = array_data_out
    )
  }

  npar_tt <- tt_npar(p, ranks)
  npar_full <- dense_npar(p, d)
  npar_tt_intrinsic <- npar_tt - tt_gauge_dim(ranks)
  residuals_resp <- y - raw$mu

  edf <- NA_real_
  edf_tt <- NA_real_
  edf_margin <- NULL
  edf_margin_cond <- NULL
  edf_note <- "not computed"
  if (isTRUE(control$compute_edf) && !is.null(raw$penalties)) {
    w_edf <- tt_edf_weights(raw, key, y, weights = weights)
    mnames <- colnames(X)
    null_proj <- if (identical(null_space, "profiled") &&
                     !is.null(null_info$projector)) {
      null_info$projector
    } else {
      NULL
    }
    parts <- tt_joint_edf_parts(
      raw$cores, basis, raw$penalties, as.numeric(raw$lambda),
      weight = w_edf, max_npar = control$edf_max_npar,
      null_proj = null_proj, names = mnames,
      penalty_order = as.integer(penalty_order),
      cyclic = cyclic,
      method = control$edf_method %||% "tedf"
    )
    edf_margin <- parts$edf_margin
    edf_margin_cond <- tt_margin_edf_cond(
      raw$cores, basis, raw$penalties, as.numeric(raw$lambda),
      weight = w_edf, names = mnames
    )
    if (identical(null_space, "profiled") && !is.null(null_proj)) {
      edf_tt <- parts$edf
      p0 <- as.integer(null_info$design_rank %||% null_info$projector$rank)
      edf <- if (is.finite(edf_tt)) p0 + edf_tt else NA_real_
      if (is.finite(edf)) {
        edf_note <- sprintf(
          paste0(
            "GDF_total = rank(X0)=%d + GDF_TT_perp=%.2f (Q0-residualized); ",
            "edf_margin = block tr(H_kk) of TT-perp (sums to edf_tt); ",
            "edf_margin_cond = ALS/cGCV diagnostic"
          ),
          p0, edf_tt
        )
      } else {
        edf_note <- "profiled GDF failed or skipped"
      }
    } else {
      edf <- parts$edf
      edf_tt <- edf
      if (is.finite(edf)) {
        edf_note <- paste0(
          "T-EDF: left-orthogonal cores + A'S_lambda A; ",
          "edf_margin = block tr(H_kk) (sums to edf); ",
          "edf_margin_cond = ALS/cGCV diagnostic"
        )
      } else if (npar_tt > control$edf_max_npar) {
        edf_note <- sprintf(
          "joint skipped (npar_TT=%d > edf_max_npar=%d); margin EDF may still be set",
          npar_tt, control$edf_max_npar
        )
      } else {
        edf_note <- "failed or skipped (see control$compute_edf / edf_max_npar)"
      }
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
      X = X,
      d = d,
      n = length(y),
      k = p,
      degree = bs$degree,
      knots = bs$knots,
      cyclic = cyclic,
      penalty_order = as.integer(penalty_order),
      penalty_mode = normalize_penalty_mode(
        raw$penalty_mode %||% control$penalty_mode %||% "global"
      ),
      cores = raw$cores,
      penalties = raw$penalties,
      rank = ranks,
      rank_internal = ranks[-c(1L, length(ranks))],
      rank_max = max(ranks),
      lambda = as.numeric(raw$lambda),
      lambda_method = raw$method_lambda,
      lambda_bounds = lam_info$lambda_bounds,
      lambda_boundary = lam_info$lambda_boundary,
      lambda_at_boundary = lam_info$lambda_at_boundary,
      intercept = raw$intercept %||% 0,
      beta = raw$beta %||% numeric(0),
      linear = linear,
      smooth = raw$smooth %||% smooth,
      offset = offset,
      null_space = null_space,
      null_space_info = null_info,
      weights = weights,
      fitted.values = raw$mu,
      linear.predictors = raw$eta,
      residuals = residuals_resp,
      deviance = raw$deviance,
      edf = edf,
      edf_tt = edf_tt,
      edf_margin = edf_margin,
      edf_margin_cond = edf_margin_cond,
      edf_note = edf_note,
      npar_tt = npar_tt,
      npar_tt_intrinsic = npar_tt_intrinsic,
      npar_dense = npar_full,
      compression_ratio = npar_full / max(npar_tt, 1),
      inference = NULL,
      ._inf = new.env(parent = emptyenv()),
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
      q_descent = raw$q_descent,
      cgcv = raw$cgcv,
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

#' @noRd
#' @export
ttpspline <- ttps

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
                                init_cores, offset = NULL, weights = NULL,
                                linear = NULL, smooth = NULL,
                                array_data = NULL) {
  offset <- normalize_offset(offset, length(y))
  weights <- normalize_weights(weights, length(y))
  linear <- normalize_linear(linear, length(y))
  smooth <- normalize_smooth(smooth, length(y))
  if ((!is.null(linear) || !is.null(smooth)) &&
      optimizer %in% c("Adam", "hybrid", "LBFGS", "GD",
                       "Damped-Newton-ALS", "LBFGS-ALS")) {
    stop(
      "`linear=` / `smooth=` are not supported with optimizer='", optimizer, "'.",
      call. = FALSE
    )
  }
  if (identical(optimizer, "Adam")) {
    return(tt_adam_fit(
      y, basis, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores, family = if (identical(key, "gaussian")) NULL else fam,
      offset = offset, weights = weights
    ))
  }

  if (identical(optimizer, "hybrid")) {
    return(tt_hybrid_fit(
      y, basis, fam, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores, offset = offset, weights = weights
    ))
  }

  if (identical(optimizer, "LBFGS")) {
    return(tt_lbfgs_fit(
      y, basis, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores,
      family = if (identical(key, "gaussian")) NULL else fam,
      offset = offset, weights = weights
    ))
  }

  if (identical(optimizer, "GD")) {
    return(tt_gd_fit(
      y, basis, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores,
      family = if (identical(key, "gaussian")) NULL else fam,
      offset = offset, weights = weights
    ))
  }

  if (identical(optimizer, "Damped-Newton-ALS")) {
    if (identical(key, "gaussian")) {
      # closed-form ALS is the natural Gaussian path; DN reduces to ridge ALS
      warning("Damped-Newton-ALS on Gaussian uses the same path as ALS.",
              call. = FALSE)
      out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                       init_cores = init_cores, offset = offset, weights = weights,
                       linear = linear, smooth = smooth)
      out$optimizer <- "Damped-Newton-ALS"
      out$backend <- "R"
      return(out)
    }
    return(tt_damped_newton_als_fit(
      y, basis, fam, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores, offset = offset, weights = weights
    ))
  }

  if (identical(optimizer, "LBFGS-ALS")) {
    if (identical(key, "gaussian")) {
      warning("LBFGS-ALS on Gaussian falls back to ALS (closed form).",
              call. = FALSE)
      out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                       init_cores = init_cores, offset = offset, weights = weights,
                       linear = linear, smooth = smooth)
      out$optimizer <- "LBFGS-ALS"
      out$backend <- "R"
      return(out)
    }
    return(tt_lbfgs_als_fit(
      y, basis, fam, ranks, lambda_spec, control, penalty_order,
      init_cores = init_cores, offset = offset, weights = weights
    ))
  }

  # ALS / PIRLS-ALS (default structure-aware path)
  if (identical(key, "gaussian")) {
    if (identical(backend, "Rcpp") && is.null(linear) && is.null(smooth) &&
        is.null(array_data)) {
      tt_als_fit_rcpp(y, basis, ranks, lambda_spec, control, penalty_order,
                      init_cores = init_cores, offset = offset, weights = weights,
                      linear = linear, smooth = smooth)
    } else {
      out <- tt_als_fit(y, basis, ranks, lambda_spec, control, penalty_order,
                       init_cores = init_cores, offset = offset, weights = weights,
                       linear = linear, smooth = smooth,
                       array_data = array_data)
      out$backend <- if (!is.null(array_data)) "R-array" else "R"
      out
    }
  } else {
    if (identical(backend, "Rcpp") && is.null(linear) && is.null(smooth)) {
      tt_pirls_fit_rcpp(y, basis, fam, ranks, lambda_spec, control, penalty_order,
                        init_cores = init_cores, offset = offset, weights = weights,
                        linear = linear, smooth = smooth)
    } else {
      tt_pirls_fit(y, basis, fam, ranks, lambda_spec, control, penalty_order,
                  init_cores = init_cores, offset = offset, weights = weights,
                  linear = linear, smooth = smooth)
    }
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
