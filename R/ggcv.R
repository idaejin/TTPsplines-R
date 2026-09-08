# Public TT-gGCV (joint global GCV) API.
#
# Default product selector remains lambda = "cGCV". Joint gGCV is opt-in:
# expensive Sobol + SAA Monte Carlo GDF search (lab stack), Gaussian only.

#' Joint TT-gGCV smoothing-parameter search
#'
#' Optimizes all directional \(\boldsymbol\lambda\) jointly by a Monte Carlo
#' estimate of global GCV (Hutchinson-style GDF of the fixed-rank TT ALS map).
#' This is the experimental **TT-gGCV** oracle; it does **not** replace the
#' default conditional selector `lambda = "cGCV"`.
#'
#' @param y,X Gaussian response and covariates.
#' @param formula,data Optional formula interface; ignored if `y`/`X` supplied.
#' @param rank Fixed TT rank (scalar).
#' @param family Currently only [stats::gaussian()].
#' @param theta_lower,theta_upper Box in \(\log_{10}\lambda\) (default from
#'   `control$lambda_bounds`).
#' @param n_global Sobol sample size (`NULL` → dimension default).
#' @param n_refine Number of local `nlminb` refinements.
#' @param n_diverse Diverse Sobol elites before refine.
#' @param M_search,M_final Monte Carlo GDF probe counts (search / final).
#' @param core_starts_search,core_starts_final ALS inits per evaluation.
#' @param min_dist,boundary_tol Diversity / boundary diagnostics.
#' @param seed Master seed.
#' @param k,degree,penalty_order Basis / penalty.
#' @param control [tt_control()] for ALS budgets.
#' @param include_cgcv_anchor Seed with a cGCV fit.
#' @param extra_theta Optional extra \(\theta=\log_{10}\lambda\) starts.
#' @param epsilon_rel,scheme Forwarded to Monte Carlo GDF.
#' @param fit_backend Lab ALS backend (`"R"` / `"Rcpp_fixed"`).
#' @param adaptive_fidelity,fidelity,gdf_init Lab fidelity / GDF init knobs.
#' @param verbose Print stage progress.
#' @param ... Unused (error if supplied).
#' @return A list with `lambda`, `theta`, `gcv`, `gdf`, `fit`, `boundary`, and
#'   search diagnostics. The embedded `fit` is a [ttps()] object.
#'
#' @seealso [ttps()] with `lambda = "gGCV"`, [tt_control()] (`ggcv_*` knobs).
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 200
#' X <- cbind(runif(n), runif(n))
#' y <- sin(2 * pi * X[, 1]) + 0.3 * sin(pi * X[, 2]) + rnorm(n, sd = 0.25)
#' opt <- tt_ggcv(y, X, rank = 2, k = 6,
#'                n_global = 16, M_search = 4, M_final = 8,
#'                control = tt_control(max_sweeps = 8))
#' opt$lambda
#' }
#' @export
tt_ggcv <- function(y = NULL,
                    X = NULL,
                    formula = NULL,
                    data = NULL,
                    rank,
                    family = stats::gaussian(),
                    theta_lower = NULL,
                    theta_upper = NULL,
                    n_global = NULL,
                    n_refine = 5L,
                    n_diverse = 10L,
                    M_search = 15L,
                    M_final = 40L,
                    core_starts_search = 1L,
                    core_starts_final = 3L,
                    min_dist = 0.75,
                    boundary_tol = 0.05,
                    seed = 1L,
                    k = 8L,
                    degree = 3L,
                    penalty_order = 2L,
                    control = tt_control(max_sweeps = 40L,
                                         compute_edf = FALSE,
                                         seed = 1L),
                    include_cgcv_anchor = TRUE,
                    extra_theta = NULL,
                    epsilon_rel = 1e-3,
                    scheme = c("forward", "central"),
                    fit_backend = c("R", "Rcpp_fixed"),
                    adaptive_fidelity = TRUE,
                    fidelity = NULL,
                    gdf_init = c("probe_warm", "cold"),
                    verbose = FALSE,
                    ...) {
  if (length(list(...))) {
    stop("Unused arguments in tt_ggcv(): ",
         paste(names(list(...)), collapse = ", "), call. = FALSE)
  }
  # Map lambda_bounds -> log10 box when caller omits theta_* .
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  if (is.null(theta_lower)) theta_lower <- log10(bounds[1L])
  if (is.null(theta_upper)) theta_upper <- log10(bounds[2L])

  out <- tt_global_lambda_optimize(
    y = y,
    X = X,
    formula = formula,
    data = data,
    rank = rank,
    family = family,
    theta_lower = theta_lower,
    theta_upper = theta_upper,
    n_global = n_global,
    n_refine = n_refine,
    n_diverse = n_diverse,
    M_search = M_search,
    M_final = M_final,
    core_starts_search = core_starts_search,
    core_starts_final = core_starts_final,
    min_dist = min_dist,
    boundary_tol = boundary_tol,
    seed = seed,
    k = k,
    degree = degree,
    penalty_order = penalty_order,
    control = control,
    include_cgcv_anchor = include_cgcv_anchor,
    extra_theta = extra_theta,
    epsilon_rel = epsilon_rel,
    scheme = scheme,
    fit_backend = fit_backend,
    adaptive_fidelity = adaptive_fidelity,
    fidelity = fidelity,
    gdf_init = gdf_init,
    verbose = verbose
  )
  out$method <- "gGCV"
  if (!is.null(out$fit) && inherits(out$fit, c("ttps", "ttpspline"))) {
    out$fit$lambda_method <- "gGCV"
    out$fit$ggcv <- list(
      gcv = out$gcv,
      gdf = out$gdf,
      gdf_mc_se = out$gdf_mc_se,
      theta = out$theta,
      boundary = out$boundary,
      winner_source = out$winner_source,
      cost = out$cost,
      elapsed = out$elapsed
    )
  }
  out
}

#' Resolve `lambda = "gGCV"` inside [ttps()] (Gaussian scattered only).
#' @keywords internal
#' @noRd
.ttps_dispatch_ggcv <- function(y, X, rank, k, degree, penalty_order,
                                control, init = NULL, cl = NULL) {
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  n_global <- control$ggcv_n_global
  n_refine <- as.integer(control$ggcv_n_refine %||% 5L)
  M_search <- as.integer(control$ggcv_M_search %||% 15L)
  M_final <- as.integer(control$ggcv_M_final %||% 40L)
  seed <- as.integer(control$seed %||% 1L)

  ctrl <- control
  ctrl$compute_edf <- isTRUE(control$compute_edf)
  # Optimizer path uses its own ALS budgets; keep user's max_sweeps as a floor.
  if (is.null(ctrl$max_sweeps) || ctrl$max_sweeps < 8L) {
    ctrl$max_sweeps <- 40L
  }

  opt <- tt_ggcv(
    y = y,
    X = X,
    rank = rank,
    k = k,
    degree = degree,
    penalty_order = penalty_order,
    theta_lower = log10(bounds[1L]),
    theta_upper = log10(bounds[2L]),
    n_global = n_global,
    n_refine = n_refine,
    M_search = M_search,
    M_final = M_final,
    seed = seed,
    control = ctrl,
    include_cgcv_anchor = isTRUE(control$ggcv_include_cgcv_anchor %||% TRUE),
    verbose = isTRUE(control$trace)
  )
  fit <- opt$fit
  if (is.null(fit) || !inherits(fit, c("ttps", "ttpspline"))) {
    stop("lambda = \"gGCV\" failed to produce a TT fit.", call. = FALSE)
  }
  if (!is.null(cl)) fit$call <- cl
  fit$lambda_method <- "gGCV"
  fit$ggcv <- list(
    gcv = opt$gcv,
    gdf = opt$gdf,
    gdf_mc_se = opt$gdf_mc_se,
    theta = opt$theta,
    boundary = opt$boundary,
    winner_source = opt$winner_source,
    cost = opt$cost,
    elapsed = opt$elapsed,
    diagnostics = opt$diagnostics
  )
  fit
}
