#' Control parameters for [ttps()].
#'
#' Separates algorithmic knobs from the model API. Optimizer-specific fields
#' for L-BFGS / Adam are reserved; Adam/Keras is optional and not required
#' for the default ALS + cGCV workflows.
#'
#' @param max_sweeps Maximum ALS sweeps (Gaussian / inner GLM).
#' @param pirls_maxit Maximum outer PIRLS iterations (GLM).
#' @param tol Relative deviance / RSS change for stopping.
#' @param tol_lambda Relative log-λ change for cGCV outer convergence.
#' @param lambda_start Initial λ used when `lambda = "cGCV"`.
#' @param lambda_bounds Length-2 bounds for automatic λ search on log scale.
#' @param lambda_tol Alias of `tol_lambda`.
#' @param lambda_update Reserved (`"auto"`); future releases may extend.
#' @param backend `"auto"`, `"R"`, `"Rcpp"`, or `"keras"` (Adam only).
#' @param sparse `"auto"`, `TRUE`, or `FALSE` (v0: dense bases; hybrid reserved).
#' @param use_spectral_gcv Use spectral cache inside Brent cGCV when feasible.
#' @param outer_maxit Outer alternation iters (LBFGS/Adam + cGCV).
#' @param outer_tol Outer relative change tolerance.
#' @param lbfgs_maxit Max L-BFGS iterations.
#' @param gd_lr Initial / nominal GD step size (start of Armijo search when
#'   `gd_linesearch = TRUE`).
#' @param gd_maxit Maximum gradient-descent iterations.
#' @param gd_tol Infinity-norm gradient tolerance for GD stopping.
#' @param gd_linesearch Use Armijo backtracking (recommended; default TRUE).
#' @param gd_step_factor Backtracking factor for GD line search (default 0.5).
#' @param gd_step_min Smallest GD step before declaring line-search failure.
#' @param gd_armijo_c Armijo sufficient-decrease constant (default 1e-4).
#' @param adam_lr,adam_epochs,adam_batch_size,adam_patience Adam/Keras knobs
#'   (optional backend; `adam_batch_size = NULL` means full-batch).
#' @param trace Print iteration progress (`TRUE`/`FALSE`). Alias of `monitor`.
#' @param monitor Alias of `trace` — set `monitor = TRUE` to watch ALS / PIRLS /
#'   L-BFGS / GD sweeps. If both are supplied, either `TRUE` enables logging.
#' @param damping Legacy alias: for Bernoulli, enables `pirls_step_halving`
#'   when `pirls_step_halving` is left at default.
#' @param pirls_step_halving Bernoulli outer line search on the true penalized
#'   Bernoulli objective (parameter-space blend of TT cores).
#' @param step_factor Multiplicative factor for step-halving (default 0.5).
#' @param step_min Smallest accepted step size (default 1/128).
#' @param objective_tol Absolute tolerance for accepting a non-increasing step.
#' @param binomial_mu_eps Clip Bernoulli mean away from 0/1 in PIRLS.
#' @param binomial_weight_floor Shared floor for Bernoulli working variance
#'   used in both `W` and `z`.
#' @param seed RNG seed for TT core initialization.
#' @param init_sd SD for random TT initialization.
#' @param als_sweeps_per_pirls Inner ALS sweeps per PIRLS iteration (default
#'   `1`). Extra inner sweeps multiply `P_k^{full}` builds; on Chicago Poisson
#'   \(r=3\), `1` matched deviance of `3` at ~5–6× lower wall time.
#' @param als_sweeps_adaptive If `TRUE` (default), allow additional inner ALS
#'   sweeps up to `als_sweeps_per_pirls_max` when the working predictor is still
#'   moving after the base sweep count.
#' @param als_sweeps_per_pirls_max Cap on adaptive inner sweeps (default `3`).
#' @param hybrid_lbfgs_maxit Max L-BFGS iterations after ALS warm-start when
#'   `optimizer = "hybrid"` (experimental Bernoulli polish).
#' @param dn_max_sweeps Max outer sweeps for Damped-Newton-ALS.
#' @param dn_armijo_c Armijo constant for Damped-Newton-ALS.
#' @param dn_step_factor Backtracking factor for Damped-Newton-ALS.
#' @param dn_step_min Smallest Newton step before rejecting the update.
#' @param dn_ridge Optional explicit ridge added to the conditional Hessian
#'   (separate from the P-spline penalty; default 0).
#' @param block_lbfgs_maxit Max L-BFGS iterations per core for LBFGS-ALS.
#' @param block_lbfgs_sweeps Max outer ALS sweeps for LBFGS-ALS.
#' @param compute_edf Compute joint linearized EDF after fit (`TRUE`/`FALSE`).
#'   Skipped automatically when packed TT size exceeds `edf_max_npar`.
#' @param edf_max_npar Maximum packed TT parameters for joint EDF (memory guard).
#' @param warn_lambda_boundary Soft warning when cGCV λ hits search bounds.
#' @param cgcv_update cGCV dynamics: `"outer_simultaneous"` (default; fit all
#'   cores at fixed λ, freeze, Jacobi proposals, damped/trust update) or
#'   `"sequential"` (legacy Gauss–Seidel; can cascade to λ_max on Chicago
#'   Poisson under the global penalty).
#' @param cgcv_damping Log-scale mixing weight \(\rho\in(0,1]\) toward the
#'   raw proposal (default `0.25`).
#' @param cgcv_max_log10_step Trust region: max \(|\Delta\log_{10}\lambda|\)
#'   per update (default `1`). Use `Inf` to disable.
#' @param cgcv_parameterization `"free"` (default) or `"scale_anisotropy"`
#'   (\(\lambda_m=\lambda_0\omega_m\), \(\prod\omega_m=1\)).
#' @param cgcv_trace Store per-update cGCV diagnostics on the fit (`TRUE`).
#' @param cgcv_margin_order Optional permutation of `1:d` for sequential
#'   updates (order-sensitivity diagnostics).
#' @param cgcv_fit_sweeps ALS sweeps in each outer fit step (`NULL` →
#'   `max_sweeps`).
#' @param cgcv_lambda0_method For scale–anisotropy: `"fixed_start"` or
#'   `"log_grid"` overall-scale search.
#' @param cgcv_lambda0_grid Optional numeric grid for `"log_grid"`.
#' @return A list of class `"tt_control"`.
#' @export
tt_control <- function(max_sweeps = 50,
                       pirls_maxit = 50,
                       tol = 1e-8,
                       tol_lambda = 1e-4,
                       lambda_start = 1,
                       lambda_bounds = c(1e-4, 1e4),
                       lambda_tol = NULL,
                       lambda_update = "auto",
                       backend = c("auto", "R", "Rcpp", "keras"),
                       sparse = c("auto", TRUE, FALSE),
                       use_spectral_gcv = TRUE,
                       outer_maxit = 20,
                       outer_tol = 1e-5,
                       lbfgs_maxit = 500,
                       gd_lr = 1e-2,
                       gd_maxit = 5000L,
                       gd_tol = 1e-7,
                       gd_linesearch = TRUE,
                       gd_step_factor = 0.5,
                       gd_step_min = 1e-12,
                       gd_armijo_c = 1e-4,
                       adam_lr = 1e-3,
                       adam_epochs = 1000,
                       adam_batch_size = NULL,
                       adam_patience = 30,
                       trace = FALSE,
                       monitor = NULL,
                       damping = TRUE,
                       pirls_step_halving = NULL,
                       step_factor = 0.5,
                       step_min = 1 / 128,
                       objective_tol = 1e-10,
                       binomial_mu_eps = 1e-5,
                       binomial_weight_floor = 1e-4,
                       seed = 1,
                       init_sd = 0.15,
                       als_sweeps_per_pirls = 1,
                       als_sweeps_adaptive = TRUE,
                       als_sweeps_per_pirls_max = 3L,
                       hybrid_lbfgs_maxit = 50L,
                       dn_max_sweeps = 40L,
                       dn_armijo_c = 1e-4,
                       dn_step_factor = 0.5,
                       dn_step_min = 1e-12,
                       dn_ridge = 0,
                       block_lbfgs_maxit = 50L,
                       block_lbfgs_sweeps = 40L,
                       compute_edf = TRUE,
                       edf_max_npar = 2500L,
                       warn_lambda_boundary = TRUE,
                       cgcv_update = c("outer_simultaneous", "sequential"),
                       cgcv_damping = 0.25,
                       cgcv_max_log10_step = 1,
                       cgcv_parameterization = c("free", "scale_anisotropy"),
                       cgcv_trace = TRUE,
                       cgcv_margin_order = NULL,
                       cgcv_fit_sweeps = NULL,
                       cgcv_lambda0_method = c("fixed_start", "log_grid"),
                       cgcv_lambda0_grid = NULL) {
  backend <- match.arg(backend)
  cgcv_update <- match.arg(cgcv_update)
  cgcv_parameterization <- match.arg(cgcv_parameterization)
  cgcv_lambda0_method <- match.arg(cgcv_lambda0_method)
  if (is.character(sparse) && length(sparse) == 1L) {
    sparse <- match.arg(sparse, c("auto", "TRUE", "FALSE", "true", "false"))
    if (sparse %in% c("TRUE", "true")) sparse <- TRUE
    else if (sparse %in% c("FALSE", "false")) sparse <- FALSE
  } else if (!is.logical(sparse) && !identical(sparse, "auto")) {
    sparse <- "auto"
  }
  if (!is.null(lambda_tol)) tol_lambda <- lambda_tol
  if (is.null(pirls_step_halving)) pirls_step_halving <- isTRUE(damping)
  # monitor is a user-facing alias of trace (either TRUE enables logging)
  if (!is.null(monitor)) {
    trace <- isTRUE(trace) || isTRUE(monitor)
  }
  structure(
    list(
      max_sweeps = as.integer(max_sweeps),
      pirls_maxit = as.integer(pirls_maxit),
      tol = as.numeric(tol),
      tol_lambda = as.numeric(tol_lambda),
      lambda_start = as.numeric(lambda_start),
      lambda_bounds = as.numeric(lambda_bounds),
      lambda_update = as.character(lambda_update),
      backend = backend,
      sparse = sparse,
      use_spectral_gcv = isTRUE(use_spectral_gcv),
      outer_maxit = as.integer(outer_maxit),
      outer_tol = as.numeric(outer_tol),
      lbfgs_maxit = as.integer(lbfgs_maxit),
      gd_lr = as.numeric(gd_lr),
      gd_maxit = as.integer(gd_maxit),
      gd_tol = as.numeric(gd_tol),
      gd_linesearch = isTRUE(gd_linesearch),
      gd_step_factor = as.numeric(gd_step_factor),
      gd_step_min = as.numeric(gd_step_min),
      gd_armijo_c = as.numeric(gd_armijo_c),
      adam_lr = as.numeric(adam_lr),
      adam_epochs = as.integer(adam_epochs),
      adam_batch_size = if (is.null(adam_batch_size)) NULL else as.integer(adam_batch_size),
      adam_patience = as.integer(adam_patience),
      trace = isTRUE(trace),
      monitor = isTRUE(trace),
      damping = isTRUE(damping),
      pirls_step_halving = isTRUE(pirls_step_halving),
      step_factor = as.numeric(step_factor),
      step_min = as.numeric(step_min),
      objective_tol = as.numeric(objective_tol),
      binomial_mu_eps = as.numeric(binomial_mu_eps),
      binomial_weight_floor = as.numeric(binomial_weight_floor),
      seed = as.integer(seed),
      init_sd = as.numeric(init_sd),
      als_sweeps_per_pirls = as.integer(als_sweeps_per_pirls),
      als_sweeps_adaptive = isTRUE(als_sweeps_adaptive),
      als_sweeps_per_pirls_max = as.integer(als_sweeps_per_pirls_max),
      hybrid_lbfgs_maxit = as.integer(hybrid_lbfgs_maxit),
      dn_max_sweeps = as.integer(dn_max_sweeps),
      dn_armijo_c = as.numeric(dn_armijo_c),
      dn_step_factor = as.numeric(dn_step_factor),
      dn_step_min = as.numeric(dn_step_min),
      dn_ridge = as.numeric(dn_ridge),
      block_lbfgs_maxit = as.integer(block_lbfgs_maxit),
      block_lbfgs_sweeps = as.integer(block_lbfgs_sweeps),
      compute_edf = isTRUE(compute_edf),
      edf_max_npar = as.integer(edf_max_npar),
      warn_lambda_boundary = isTRUE(warn_lambda_boundary),
      cgcv_update = cgcv_update,
      cgcv_damping = as.numeric(cgcv_damping),
      cgcv_max_log10_step = as.numeric(cgcv_max_log10_step),
      cgcv_parameterization = cgcv_parameterization,
      cgcv_trace = isTRUE(cgcv_trace),
      cgcv_margin_order = if (is.null(cgcv_margin_order)) NULL else as.integer(cgcv_margin_order),
      cgcv_fit_sweeps = if (is.null(cgcv_fit_sweeps)) NULL else as.integer(cgcv_fit_sweeps),
      cgcv_lambda0_method = cgcv_lambda0_method,
      cgcv_lambda0_grid = if (is.null(cgcv_lambda0_grid)) NULL else as.numeric(cgcv_lambda0_grid),
      # Classical multidimensional P-spline penalty on Θ only (no surrogate).
      penalty_mode = "global"
    ),
    class = "tt_control"
  )
}

#' Package always uses the classical global penalty on Θ.
#' Removed surrogates (`own_margin` / `separable`) error if requested.
#' @keywords internal
#' @noRd
normalize_penalty_mode <- function(mode) {
  if (is.null(mode) || !length(mode)) return("global")
  mode <- as.character(mode)[[1L]]
  if (mode %in% c("global", "full")) return("global")
  if (mode %in% c("own_margin", "separable")) {
    stop(
      "penalty_mode '", mode, "' was removed: it is not the classical ",
      "multidimensional P-spline penalty on Θ. ",
      "TTPsplines always uses J_λ(Θ) = sum_m λ_m ||Θ ×_m Δ||_F^2.",
      call. = FALSE
    )
  }
  stop("Unknown penalty_mode '", mode, "'.", call. = FALSE)
}

#' @keywords internal
#' @noRd
is_global_penalty_mode <- function(mode) {
  identical(normalize_penalty_mode(mode %||% "global"), "global")
}

#' Resolve computational backend given optimizer preference.
#' @keywords internal
resolve_backend <- function(control, optimizer = "ALS") {
  be <- control$backend
  # Global P_k^full depends on other cores → ALS/PIRLS sweeps stay in R.
  if (optimizer %in% c("ALS", "PIRLS", "auto")) {
    if (identical(be, "Rcpp")) {
      warning(
        "ALS/PIRLS sweeps use the classical global penalty on Θ and run in R; ",
        "Rcpp accelerates P_k^full / global-penalty helpers when available.",
        call. = FALSE
      )
    }
    if (identical(be, "auto") || identical(be, "Rcpp")) {
      be <- "R"
    }
  }
  if (identical(be, "auto")) {
    if (identical(optimizer, "Adam")) {
      if (isTRUE(tt_has_keras())) return("keras")
      return("keras") # caller will error with install hint
    }
    # Direct / conditional-likelihood R paths
    if (optimizer %in% c("LBFGS", "GD", "hybrid",
                         "Damped-Newton-ALS", "LBFGS-ALS")) {
      return("R")
    }
    if (.ttpsplines_has_rcpp()) return("Rcpp")
    return("R")
  }
  if (identical(be, "Rcpp") && !.ttpsplines_has_rcpp()) {
    warning("Rcpp backend unavailable; falling back to R.", call. = FALSE)
    return("R")
  }
  if (identical(be, "keras") && !identical(optimizer, "Adam")) {
    warning("backend='keras' is only used with optimizer='Adam'; ignoring.",
            call. = FALSE)
    if (.ttpsplines_has_rcpp()) return("Rcpp")
    return("R")
  }
  be
}

.ttpsplines_has_rcpp <- function() {
  exists("tt_fit_d_cpp", mode = "function") ||
    exists("tt_cgcv_fit_cpp", mode = "function")
}
