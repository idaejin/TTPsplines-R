# Poisson TT-gGCV for scattered rows. tt_ggcv_poisson() is a thin wrapper
# over tt_ggcv(family = poisson()): the shared engine (ggcv_engine.R) on the
# exact fixed-lambda PIRLS map, UBRE = D + 2 GDF (scale 1) by default.
#
# The former modes are gone. "working" (Gaussian gGCV on the row-scaled
# working response; it scaled the covariates instead of the basis rows, so
# it was not the WLS problem, and dropped cyclic margins) now warns and is
# ignored. "algorithmic" (probes clipped at 0, a truncated search budget,
# mixed-fidelity comparisons, Sobol-only anisotropy) is replaced by the
# engine, which is what ggcv_glm_mode = "algorithmic" now means; its
# evaluator was removed (the archival lab review scripts of 2026-10-01 that
# call it reproduce at package commit 77d5fc6).

#' Poisson TT-gGCV smoothing selection for scattered data
#'
#' Thin wrapper over [tt_ggcv()] with `family = poisson()`: global UBRE
#' (deviance + 2 GDF, scale 1) of the fixed-\eqn{\lambda} PIRLS-ALS map, with
#' a Monte Carlo GDF from non-negative probes (no clipping at zero counts)
#' and the search of [tt_ggcv_array()]. Search knobs left `NULL` here are
#' read from the `ggcv_*` fields of `control` (see [tt_control()]), then
#' from the defaults of [tt_ggcv()].
#'
#' @inheritParams tt_ggcv
#' @param family Must be [stats::poisson()].
#' @param mode Deprecated. The exact fixed-\eqn{\lambda} map (`"algorithmic"`)
#'   is always used; `"working"` (here or in `control$ggcv_glm_mode`) gives a
#'   warning and is ignored: the working-response proxy was removed.
#' @param anisotropic If `TRUE`, one smoothing group per margin
#'   (`groups = seq_len(d)`); if `FALSE`, isotropic. `NULL` (default) reads
#'   `control$ggcv_groups`, then `control$ggcv_poisson_anisotropic` (one
#'   group per margin when `TRUE`), in the order of `ttps(lambda = "gGCV")`;
#'   isotropic when both are `NULL`.
#' @param control [tt_control()] for every fixed-\eqn{\lambda} fit (its
#'   `pirls_maxit` is the cap of the adaptive budget or the fixed budget;
#'   default cap 60 PIRLS iterations) and source of the `ggcv_*` knobs.
#' @param n_global,n_refine,M_search,M_final Search settings of [tt_ggcv()];
#'   `NULL` reads `control$ggcv_<name>`, then the [tt_ggcv()] default.
#' @param budget,fit_tol Iteration budget of [tt_ggcv()] (`"adaptive"` or
#'   `"fixed"`) and the convergence tolerance of the adaptive budget (see
#'   [tt_ggcv_array()]); `NULL` reads `control$ggcv_budget` /
#'   `control$ggcv_fit_tol`, then the defaults `"adaptive"` / `1e-7`.
#' @param seed Probe seed; `NULL` = `control$seed`.
#' @return The value of [tt_ggcv()] plus `ggcv_glm_mode = "algorithmic"`.
#' @export
tt_ggcv_poisson <- function(y,
                            X,
                            rank,
                            family = stats::poisson(),
                            offset = NULL,
                            cyclic = NULL,
                            period = NULL,
                            knots = NULL,
                            mode = NULL,
                            anisotropic = NULL,
                            k = 8L,
                            degree = 3L,
                            penalty_order = 2L,
                            control = tt_control(max_sweeps = 400L,
                                                 pirls_maxit = 60L,
                                                 compute_edf = FALSE,
                                                 seed = 1L),
                            n_global = NULL,
                            n_refine = NULL,
                            M_search = NULL,
                            M_final = NULL,
                            seed = NULL,
                            budget = NULL,
                            fit_tol = NULL,
                            gdf_method = NULL,
                            verbose = FALSE) {
  fam <- normalize_family(family)
  if (!identical(family_key(fam), "poisson")) {
    stop("tt_ggcv_poisson() currently supports family = poisson() only.",
         call. = FALSE)
  }
  mode <- match.arg(mode %||% control$ggcv_glm_mode %||% "algorithmic",
                    c("algorithmic", "working"))
  if (identical(mode, "working")) .ggcv_warn_working("tt_ggcv_poisson()")
  # resolved here (an explicit `mode` overrides the knob): tt_ggcv() must
  # not warn a second time
  control$ggcv_glm_mode <- "algorithmic"
  # same order as ttps(lambda = "gGCV"): anisotropic, ggcv_groups, then
  # ggcv_poisson_anisotropic
  groups <- .ggcv_poisson_groups(control, ncol(as.matrix(X)), anisotropic)
  kb <- .ggcv_control_knobs(control)
  out <- tt_ggcv(
    y = y, X = X, rank = rank, family = fam, offset = offset, cyclic = cyclic,
    period = period, knots = knots, k = k, degree = degree,
    penalty_order = penalty_order, groups = groups, n_grid = kb$n_grid,
    n_global = n_global %||% kb$n_global, n_refine = n_refine %||% kb$n_refine,
    n_final = kb$n_final, M_search = M_search %||% kb$M_search,
    M_final = M_final %||% kb$M_final, tol = kb$tol,
    criterion = kb$criterion, probe_init = kb$probe_init,
    n_cores = kb$n_cores, seed = seed %||% control$seed %||% 1L,
    control = control, budget_check = kb$budget_check,
    budget_fallback = kb$budget_fallback, budget = budget %||% kb$budget,
    fit_tol = fit_tol %||% kb$fit_tol,
    gdf_method = gdf_method %||% kb$gdf_method, verbose = verbose
  )
  out$ggcv_glm_mode <- "algorithmic"
  out
}

#' Deprecation warning for the removed working-response proxy.
#' @keywords internal
#' @noRd
.ggcv_warn_working <- function(where) {
  warning(
    where, ": mode \"working\" (ggcv_glm_mode) is deprecated and ignored. ",
    "The working-response proxy was removed (it row-scaled the covariates ",
    "instead of the basis rows); gGCV uses the exact fixed-lambda PIRLS map.",
    call. = FALSE
  )
}
