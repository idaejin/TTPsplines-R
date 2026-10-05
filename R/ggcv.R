# Public TT-gGCV for scattered rows (Gaussian and Poisson) and the ttps()
# dispatch of lambda = "gGCV" (scattered and array data).
#
# tt_ggcv() builds the fixed-lambda map of ttps() on the given rows
# (.ggcv_map_scattered()) and runs the engine of ggcv_engine.R, the one used
# by tt_ggcv_array(): non-negative Rademacher probes against a shared
# reference, common random numbers, a per-lambda iteration budget (adaptive:
# a pilot finds where the penalized objective settles, n_star, and the
# estimator, reference and probes run B = ceiling(1.5 n_star) + 2 iterations)
# from a cold common initialization, probe validation, an isotropic-first
# search, paired-SE thresholds and the noise-aware final choice with a
# decision-aware budget check. tt_gdf()
# scores one lambda with the same evaluator (scattered counterpart of
# tt_gdf_array()). criterion = "auto" gives GCV for Gaussian and
# UBRE = D + 2 GDF (scale 1) for Poisson. Rows carry unit weight;
# Berman-Turner quadrature weights enter through the offset. Extra starts
# come only from the caller (extra_theta) or from Sobol points (n_global).
#
# ttps(lambda = "gGCV") reads the ggcv_* knobs of tt_control() (array mode:
# ggcv_array_* first), selects with the budget ggcv_budget (caps
# ggcv_max_sweeps / ggcv_pirls_maxit) and the cold common initialization,
# then refits the model at the selected lambda with the caller's own
# arguments and control (ggcv_refit = TRUE). The default selector of ttps()
# remains "cGCV".

#' Joint TT-gGCV smoothing-parameter search for scattered data
#'
#' Selects the smoothing vector \eqn{\lambda} of a scattered-data TT
#' P-spline ([ttps()], Gaussian or Poisson) by minimising a global criterion
#' of the fixed-\eqn{\lambda} fit: GCV \eqn{n D / (n - \mathrm{GDF})^2} or
#' UBRE \eqn{D / n + 2 \phi\,\mathrm{GDF} / n - \phi}, where \eqn{D} is the
#' deviance and GDF the Monte Carlo trace of the derivative of the fitting
#' map \eqn{y \mapsto \hat\mu(y; \lambda)}. `criterion = "auto"` takes GCV
#' for Gaussian data and UBRE with \eqn{\phi = 1} (deviance + 2 GDF) for
#' Poisson data. Offsets (e.g. Berman-Turner quadrature rows) are
#' supported; observation weights are not.
#'
#' The evaluator and the search are those of [tt_ggcv_array()] (see its
#' sections *Search*, *Perturbations* and *Iteration budget*): an isotropic
#' grid with halvings and one parabolic step, an optional grouped pattern
#' search (`groups`), and a final stage at `M_final` probes. There the
#' winner is the most regular converged candidate with a stable GDF
#' (isotropic before grouped, then the largest \eqn{\lambda}) that the
#' minimum does not beat by more than `max(tol, 2 * SE)` (paired SE over the
#' common probes). Its GDF is then re-estimated at twice its iteration
#' budget; when it moves by more than `budget_tol`, the candidate is
#' rescored with the less favourable of its two scores and kept if it still
#' wins, or when no other converged candidate with a stable GDF is left to
#' check (`budget_fallback`). On the Dette d = 8 example with
#' Berman-Turner rows the plain argmin sat on the grid boundary, where the
#' Monte Carlo GDF was 432 \eqn{\pm} 334, against 145 \eqn{\pm} 10 at the
#' \eqn{\lambda} chosen by thinning CV. Iteration budgets are adaptive by
#' default (`budget`). This is not a replacement for the default selector
#' `lambda = "cGCV"`.
#'
#' @param y,X Response and covariates (`n x d`, `d >= 2`).
#' @param formula,data Optional formula interface (`y ~ x1 + x2`); used when
#'   `y` or `X` is `NULL`. The intercept column is dropped.
#' @param rank TT rank, scalar or vector as in [ttps()] (not collapsed).
#' @param family [stats::gaussian()] (default) or [stats::poisson()].
#' @param offset Optional offset (length `n` or scalar, linear-predictor
#'   scale), e.g. `log(exposure)` or Berman-Turner quadrature log-weights.
#' @param cyclic,period,knots Passed to [ttps()].
#' @param k,degree,penalty_order Basis / penalty, passed to [ttps()].
#' @param extra_theta Optional extra starts in \eqn{\log_{10}\lambda}: a
#'   vector is one start; a matrix has one row per start and one column per
#'   group or per margin (margins are averaged within groups). Evaluated with
#'   the isotropic grid when there is one group, as starts of the grouped
#'   stage otherwise.
#' @param n_global Sobol starts for the grouped stage (default `0`; `NULL`
#'   is `0`). Ignored with one group.
#' @param seed Probe seed (any finite number; probe `j` uses
#'   `seed * 100003 + j`, wrapped into the integer range when it leaves it).
#'   The TT cores of every fit start from `control$seed`.
#' @param control [tt_control()] for every fixed-\eqn{\lambda} fit; its
#'   `max_sweeps` (Gaussian) / `pirls_maxit` (Poisson) is the cap of the
#'   adaptive budget, or the fixed budget with `budget = "fixed"` (default
#'   caps: 400 sweeps / 60 PIRLS iterations). Its `tol` is not used (see
#'   `fit_tol`).
#' @param ... Not used (error if supplied).
#' @inheritParams tt_ggcv_array
#' @return The list returned by [tt_ggcv_array()] (`lambda`, `theta`,
#'   `score`, `score_dev`, `gcv`, `ubre`, `gdf`, `gdf_se`, `gdf_cv`,
#'   `stable`, `deviance`, `n_eff`, `criterion`, `groups`, `decision`,
#'   `rule`, `best_theta`, `boundary`, `budget`, `budget_checks`,
#'   `budget_verified`, `paired`, `fit` with `$ggcv`, `search`, `n_evals`,
#'   `n_fits`, `method`, `elapsed`) plus `gdf_mc_se` (alias of `gdf_se`).
#'
#' @seealso [ttps()] with `lambda = "gGCV"` (knobs `ggcv_*` of
#'   [tt_control()]), [tt_gdf()] to score a given \eqn{\lambda},
#'   [tt_ggcv_array()] for product grids, [tt_ggcv_poisson()].
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 200
#' X <- cbind(runif(n), runif(n))
#' y <- sin(2 * pi * X[, 1]) + 0.3 * sin(pi * X[, 2]) + rnorm(n, sd = 0.25)
#' opt <- tt_ggcv(y, X, rank = 2, k = 6, M_search = 4, M_final = 8)
#' opt$lambda
#' opt$paired   # final candidates and the choice
#' head(opt$search)
#'
#' ## Poisson counts with an exposure offset, one smoothing group per margin
#' expo <- runif(n, 0.5, 2)
#' yp <- rpois(n, expo * exp(sin(2 * pi * X[, 1])))
#' opt_p <- tt_ggcv(yp, X, rank = 2, k = 6, family = poisson(),
#'                  offset = log(expo), groups = 1:2)
#' opt_p$decision
#' opt_p$budget_verified
#' }
#' @export
tt_ggcv <- function(y = NULL,
                    X = NULL,
                    formula = NULL,
                    data = NULL,
                    rank,
                    family = stats::gaussian(),
                    offset = NULL,
                    cyclic = NULL,
                    period = NULL,
                    knots = NULL,
                    k = 8L,
                    degree = 3L,
                    penalty_order = 2L,
                    groups = NULL,
                    theta_lower = NULL,
                    theta_upper = NULL,
                    extra_theta = NULL,
                    n_grid = 9L,
                    n_global = 0L,
                    n_refine = 2L,
                    n_final = 3L,
                    M_search = 4L,
                    M_final = 16L,
                    tol = 1,
                    criterion = c("auto", "gcv", "ubre"),
                    scale = NULL,
                    probe_init = c("cold", "warm"),
                    probe_budget = NULL,
                    epsilon_rel = 1e-3,
                    n_cores = 1L,
                    seed = 1L,
                    control = tt_control(max_sweeps = 400L,
                                         pirls_maxit = 60L,
                                         compute_edf = FALSE,
                                         seed = 1L),
                    budget_check = TRUE,
                    budget_tol = 0.01,
                    budget_fallback = TRUE,
                    budget = c("adaptive", "fixed"),
                    fit_tol = 1e-7,
                    gdf_method = c("mc", "exact", "auto"),
                    verbose = FALSE,
                    ...) {
  t0 <- proc.time()[["elapsed"]]
  if (length(list(...))) {
    stop("Unused arguments in tt_ggcv(): ",
         paste(names(list(...)), collapse = ", "), call. = FALSE)
  }
  if (is.null(y) || is.null(X)) {
    if (is.null(formula) || is.null(data)) {
      stop("Supply y + X or formula + data.", call. = FALSE)
    }
    mf <- stats::model.frame(formula, data = data)
    y <- stats::model.response(mf)
    X <- stats::model.matrix(formula, data = data)
    if ("(Intercept)" %in% colnames(X)) {
      X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
    }
  }
  criterion <- match.arg(criterion)
  probe_init <- match.arg(probe_init)
  budget <- match.arg(budget)
  map <- .ggcv_map_scattered(
    y = y, X = X, family = family, rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, cyclic = cyclic, period = period,
    knots = knots, offset = offset
  )
  # the removed working-response proxy (tt_ggcv_poisson() and ttps() warn
  # themselves and pass the knob on as "algorithmic")
  if (identical(map$key, "poisson") && identical(control$ggcv_glm_mode, "working")) {
    .ggcv_warn_working("tt_ggcv()")
  }
  st <- .ggcv_setup(
    map, groups = groups, criterion = criterion, scale = scale,
    probe_init = probe_init, probe_budget = probe_budget,
    epsilon_rel = epsilon_rel, n_cores = n_cores, seed = seed,
    control = control, budget = budget, fit_tol = fit_tol,
    gdf_method = match.arg(gdf_method)
  )
  starts <- .ggcv_group_starts(extra_theta, st$groups)
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  out <- .ggcv_search(
    st, lo = theta_lower %||% log10(bounds[1L]),
    hi = theta_upper %||% log10(bounds[2L]), n_grid = n_grid,
    n_refine = n_refine, n_global = n_global %||% 0L, n_final = n_final,
    M_search = M_search, M_final = M_final, tol = tol, theta_start = starts,
    budget_check = budget_check, budget_tol = budget_tol,
    budget_fallback = budget_fallback, verbose = verbose
  )
  out$gdf_mc_se <- out$gdf_se
  out$elapsed <- proc.time()[["elapsed"]] - t0
  out
}

#' Global GDF and gGCV score of a scattered-data TT fit at a fixed lambda
#'
#' Scattered-row counterpart of [tt_gdf_array()]: scores one smoothing vector
#' with the evaluator of [tt_ggcv()], i.e. the fixed-\eqn{\lambda} fit of
#' [ttps()] on these rows (ALS for Gaussian, PIRLS-ALS for Poisson), its
#' deviance, GCV / UBRE and a Monte Carlo GDF
#' \eqn{\mathrm{tr}(\partial\hat\mu / \partial y)} from non-negative
#' perturbations against a shared reference fit (see [tt_ggcv_array()],
#' sections *Perturbations* and *Iteration budget*).
#'
#' @param y,X Response and covariates (`n x d`, `d >= 2`).
#' @param lambda Smoothing vector (scalar or length d, non-negative).
#' @param family [stats::gaussian()] (default) or [stats::poisson()].
#' @param rank TT rank, scalar or vector as in [ttps()].
#' @param k,degree,penalty_order,cyclic,period,knots Passed to [ttps()].
#' @param offset Optional offset (length `n` or scalar, linear-predictor
#'   scale), e.g. Berman-Turner quadrature log-weights.
#' @param M Number of probes (ignored when `probes` is a matrix).
#' @param probes `"rademacher"` (default; probe `j` drawn from
#'   `seed * 100003 + j`, the probes of [tt_ggcv()] and [tt_gdf_array()]),
#'   `"unit"` (one perturbation per row: the exact finite-difference trace,
#'   one fit per row, small problems only), or a numeric matrix of -1 / +1
#'   values with one row per observation and one column per probe, a user
#'   probe bank (`M = ncol(probes)`).
#' @param seed Probe seed (any finite number; probe `j` uses
#'   `seed * 100003 + j`, wrapped into the integer range when it leaves it).
#'   The TT cores of every fit start from `control$seed`.
#' @param control [tt_control()] for every fit: `max_sweeps` (Gaussian) /
#'   `pirls_maxit` (Poisson) is the cap of the adaptive budget or the fixed
#'   budget (default caps: 400 sweeps, 60 PIRLS iterations); `control$seed`
#'   fixes the common cold initialization. Its `tol` is not used. The
#'   estimator is the cold [ttps()] fit with this control, `tol = 0` and the
#'   iteration budget set to B (`budget`), the same map as the probes; with
#'   `budget = "fixed"`, B is the control budget, so a `ttps()` fit that
#'   stopped after `n` iterations is reproduced by `budget = "fixed"` and
#'   that budget set to `n`.
#' @param scale Known scale for `ubre` (Poisson default `1`; Gaussian `NA`
#'   unless given).
#' @inheritParams tt_ggcv_array
#' @return A list with `gdf`, `gdf_se`, `deviance`, `gcv`, `ubre`, `n_eff`,
#'   `contrib` (per-probe contributions; per row for `"unit"`), `M_ok`,
#'   `converged` and `last_rel_change` (convergence of the estimator),
#'   `budget` (iteration budget B of the estimator, reference and probe
#'   fits), `n_iter` (adaptive: the iteration where the pilot's objective
#'   settled, the cap when it did not; fixed: iterations of the estimator),
#'   `gdf_cv` (per-probe coefficient of variation; `0` for `"unit"`),
#'   `stable` (the contributions form a usable GDF, see [tt_ggcv_array()];
#'   treat `FALSE` as a failed GDF), `reevals` (probe-validation rounds, each
#'   at twice the budget), `fit` (the estimator), `epsilon`, `n_fits` and
#'   `time_s`.
#' @seealso [tt_ggcv()] for the selection, [tt_gdf_array()] for product
#'   grids.
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 200
#' X <- cbind(runif(n), runif(n))
#' y <- sin(2 * pi * X[, 1]) + 0.3 * sin(pi * X[, 2]) + rnorm(n, sd = 0.25)
#' g <- tt_gdf(y, X, lambda = c(0.1, 1), rank = 2, k = 6, M = 16)
#' c(g$gdf, g$gdf_se, g$budget)
#'
#' ## a user probe bank (one column per probe)
#' R <- matrix(sample(c(-1, 1), n * 8, replace = TRUE), n, 8)
#' tt_gdf(y, X, lambda = c(0.1, 1), rank = 2, k = 6, probes = R)$gdf
#' }
#' @export
tt_gdf <- function(y,
                   X,
                   lambda,
                   family = stats::gaussian(),
                   rank,
                   k = 8L,
                   degree = 3L,
                   penalty_order = 2L,
                   cyclic = NULL,
                   period = NULL,
                   knots = NULL,
                   offset = NULL,
                   M = 16L,
                   probes = "rademacher",
                   probe_init = c("cold", "warm"),
                   probe_budget = NULL,
                   epsilon_rel = 1e-3,
                   n_cores = 1L,
                   seed = 1L,
                   budget = c("adaptive", "fixed"),
                   fit_tol = 1e-7,
                   control = tt_control(max_sweeps = 400L,
                                        pirls_maxit = 60L,
                                        compute_edf = FALSE,
                                        seed = 1L),
                   scale = NULL,
                   gdf_method = c("mc", "exact")) {
  probe_init <- match.arg(probe_init)
  budget <- match.arg(budget)
  gdf_method <- match.arg(gdf_method)
  map <- .ggcv_map_scattered(
    y = y, X = X, family = family, rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, cyclic = cyclic, period = period,
    knots = knots, offset = offset
  )
  .ggcv_gdf_point(
    map, lambda = lambda, M = M, probes = probes, scale = scale,
    probe_init = probe_init, probe_budget = probe_budget,
    epsilon_rel = epsilon_rel, n_cores = n_cores, seed = seed,
    budget = budget, fit_tol = fit_tol, control = control,
    gdf_method = gdf_method
  )
}

#' Extra starts in group space: a vector is one start; a matrix has one row
#' per start and one column per group or per margin (margins averaged within
#' groups).
#' @keywords internal
#' @noRd
.ggcv_group_starts <- function(theta, groups, what = "extra_theta") {
  if (is.null(theta) || !length(theta)) return(NULL)
  G <- max(groups)
  d <- length(groups)
  m <- if (is.null(dim(theta))) {
    matrix(as.numeric(theta), nrow = 1L)
  } else {
    as.matrix(theta)
  }
  storage.mode(m) <- "double"
  if (ncol(m) == d && G < d) {
    m <- do.call(rbind, lapply(seq_len(nrow(m)), function(i) {
      as.numeric(tapply(m[i, ], groups, mean))
    }))
  }
  if (ncol(m) != G) {
    stop(sprintf("`%s` must have one column per group (%d) or per margin (%d).",
                 what, G, d), call. = FALSE)
  }
  m
}

# ---------------------------------------------------------------------------
# ttps() dispatch
# ---------------------------------------------------------------------------

#' gGCV knobs of a [tt_control()] object: `<prefix><name>` first (array
#' mode: `ggcv_array_<name>`, set with `$<-`), then the unified
#' `ggcv_<name>`, then the engine defaults. Controls without
#' `ggcv_max_sweeps` / `ggcv_pirls_maxit` select with their own budget as
#' the cap (or fixed budget).
#' @keywords internal
#' @noRd
.ggcv_control_knobs <- function(control, prefix = NULL) {
  knob <- function(name, default) {
    v <- if (is.null(prefix)) NULL else control[[paste0(prefix, name)]]
    v %||% control[[paste0("ggcv_", name)]] %||% default
  }
  list(
    groups = knob("groups", NULL),
    n_grid = as.integer(knob("n_grid", 9L)),
    n_global = as.integer(knob("n_global", 0L)),
    n_refine = as.integer(knob("n_refine", 2L)),
    n_final = as.integer(knob("n_final", 3L)),
    M_search = as.integer(knob("M_search", 4L)),
    M_final = as.integer(knob("M_final", 16L)),
    tol = as.numeric(knob("tol", 1)),
    criterion = knob("criterion", "auto"),
    probe_init = knob("probe_init", "cold"),
    n_cores = as.integer(knob("n_cores", 1L)),
    budget_check = !isFALSE(knob("budget_check", TRUE)),
    budget_fallback = !isFALSE(knob("budget_fallback", TRUE)),
    budget = as.character(knob("budget", "adaptive"))[1L],
    fit_tol = as.numeric(knob("fit_tol", 1e-7)),
    gdf_method = as.character(knob("gdf_method", "auto"))[1L],
    refit = !isFALSE(knob("refit", TRUE)),
    max_sweeps = as.integer(knob("max_sweeps", control$max_sweeps)),
    pirls_maxit = as.integer(knob("pirls_maxit", control$pirls_maxit))
  )
}

#' Smoothing groups of a Poisson gGCV selection, in the same order for
#' [tt_ggcv_poisson()] and both `ttps(lambda = "gGCV")` dispatchers
#' (scattered rows and arrays): the explicit `anisotropic` argument (`TRUE`:
#' one group per margin; `FALSE`: isotropic), then the groups knob (`groups`:
#' `control$ggcv_groups`, in array mode `ggcv_array_groups` first), then
#' `control$ggcv_poisson_anisotropic` (one group per margin when `TRUE`);
#' `NULL` (isotropic) otherwise.
#' @keywords internal
#' @noRd
.ggcv_poisson_groups <- function(control, d, anisotropic = NULL,
                                 groups = control$ggcv_groups) {
  if (!is.null(anisotropic)) return(if (isTRUE(anisotropic)) seq_len(d))
  groups %||% (if (isTRUE(control$ggcv_poisson_anisotropic)) seq_len(d))
}

#' Resolve `lambda = "gGCV"` inside [ttps()] for scattered rows (Gaussian or
#' Poisson): engine selection on the call's rows with the unified `ggcv_*`
#' knobs, then the refit of `.ttps_ggcv_finish()`.
#' @keywords internal
#' @noRd
.ttps_dispatch_ggcv_scattered <- function(y, X, family, rank, k, degree,
                                          penalty_order, cyclic, period,
                                          knots, offset, control,
                                          refit = NULL, cl = NULL) {
  fam <- normalize_family(family)
  is_pois <- identical(family_key(fam), "poisson")
  ctrl <- control
  if (is_pois && identical(control$ggcv_glm_mode, "working")) {
    .ggcv_warn_working("ttps(lambda = \"gGCV\")")
    ctrl$ggcv_glm_mode <- "algorithmic"  # warned once; tt_ggcv() stays quiet
  }
  kb <- .ggcv_control_knobs(control)
  groups <- if (is_pois) .ggcv_poisson_groups(control, ncol(X)) else kb$groups
  off <- if (is.null(offset) || all(as.numeric(offset) == 0)) NULL else offset
  ctrl$max_sweeps <- kb$max_sweeps
  ctrl$pirls_maxit <- kb$pirls_maxit
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  opt <- tt_ggcv(
    y = y, X = X, rank = rank, family = fam, offset = off, cyclic = cyclic,
    period = period, knots = knots, k = k, degree = degree,
    penalty_order = penalty_order, groups = groups,
    theta_lower = log10(bounds[1L]), theta_upper = log10(bounds[2L]),
    n_grid = kb$n_grid, n_global = kb$n_global, n_refine = kb$n_refine,
    n_final = kb$n_final, M_search = kb$M_search, M_final = kb$M_final,
    tol = kb$tol, criterion = kb$criterion, probe_init = kb$probe_init,
    n_cores = kb$n_cores, seed = as.integer(control$seed %||% 1L),
    control = ctrl, budget_check = kb$budget_check,
    budget_fallback = kb$budget_fallback, budget = kb$budget,
    fit_tol = kb$fit_tol, gdf_method = kb$gdf_method,
    verbose = isTRUE(control$trace)
  )
  .ttps_ggcv_finish(opt, refit = if (kb$refit) refit, cl = cl)
}

#' Fit returned by `ttps(lambda = "gGCV")`: `refit(lambda)` at the selected
#' lambda (the caller's own arguments and control, RNG stream preserved) or,
#' when `refit` is `NULL` (`ggcv_refit = FALSE`), the engine's estimator fit
#' at the winner (the cold fit at the winner's budget B with `tol = 0`).
#' Attaches the selection as `fit$ggcv`.
#' @keywords internal
#' @noRd
.ttps_ggcv_finish <- function(opt, refit = NULL, cl = NULL) {
  do_refit <- is.function(refit)
  fit <- if (do_refit) .tt_with_preserved_seed(refit(opt$lambda)) else opt$fit
  if (!inherits(fit, "ttpspline")) {
    stop("lambda = \"gGCV\" failed to produce a TT fit.", call. = FALSE)
  }
  fit$lambda_method <- "gGCV"
  fit$lambda_at_boundary <- isTRUE(opt$boundary)
  fit$ggcv <- c(opt$fit$ggcv, list(
    refit = do_refit, n_evals = opt$n_evals, n_fits = opt$n_fits,
    elapsed = opt$elapsed
  ))
  if (!is.null(cl)) fit$call <- cl
  fit
}
