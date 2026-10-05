# Global TT-gGCV for array (product-grid) data, ttps(array = TRUE).
# Fitting regimes (naming per 2026-08-19 DECISION): K-ALS for unweighted
# Gaussian (Kronecker Gram); weighted GLAM row-tensor PIRLS for Poisson or
# observation weights. The selector treats both as one fixed-lambda map.
#
# Criterion at a fixed smoothing vector lambda:
#   gGCV(lambda) = n * D(y, mu_hat) / (n - GDF)^2,
#   GDF = tr(d mu_hat / d y) of the fixed-lambda fitting map
#         y |-> ttps(Y, array = TRUE, lambda = lambda),
# with n = number of cells carrying positive weight. Optional UBRE form
# D / n + 2 * scale * GDF / n - scale (Poisson: scale = 1 known).
# criterion = "auto" (default) follows mgcv's GCV.Cp rule: UBRE when the
# scale is known (Poisson), GCV when it is estimated (Gaussian). GCV on a
# Poisson deviance estimates the scale from D / n, which goes wrong when the
# response is under- or near-degenerate (STATS19 occupied-only pilot).
#
# GDF is estimated by Hutchinson / Monte Carlo finite differences with
# NON-NEGATIVE perturbations, so Poisson responses stay valid on zero-heavy
# counts (no clipping at 0):
#   r ~ Rademacher on the weighted cells, z = (r + 1) / 2 in {0, 1},
#   J r = 2 J (z - 1/2) ~= (2 / eps) * (mu(y + eps z) - mu(y + eps / 2)),
#   E[r' J r] = tr J.
# One shared reference fit at y + eps / 2 replaces the base fit in the
# difference, so each probe costs one fit and the base fit's own
# convergence error cancels.
#
# Determinism: the estimator, the reference and the probes of a point run the
# SAME fixed iteration budget B (control$tol is forced to 0, so neither PIRLS
# nor ALS stops early) from the SAME initialization ("cold", default). The
# GDF is then the derivative of one smooth map, which is the regime in which
# Hutchinson GDF tracks the joint EDF (manuscript appendix on cGCV versus
# gGCV). With budget = "adaptive" (default) B follows lambda: a pilot at y
# finds the first iteration n_star after which the penalized objective has
# settled (relative change <= fit_tol at two consecutive iterations) under
# the cap control$pirls_maxit / control$max_sweeps, and
# B = min(cap, ceiling(1.5 * n_star) + 2); budget = "fixed" runs every fit at
# the control budget. probe_init = "warm" restarts probes from the
# estimator's cores with a shorter budget; cheaper, but it under-responds
# when ALS mixes slowly.
#
# The evaluator, the cache and the search live in ggcv_engine.R (shared with
# the scattered routes); this file builds the array map and the public API.

#' Global TT-gGCV smoothing selection for array data (TT-GLAM, `array = TRUE`)
#'
#' Selects \eqn{\lambda} for a product-grid TT P-spline fitted with
#' [ttps()]`(array = TRUE)` by minimising a global criterion of the
#' fixed-\eqn{\lambda} array fit (K-ALS for unweighted Gaussian data,
#' weighted GLAM row-tensor PIRLS for Poisson responses or observation
#' weights): GCV \eqn{n D / (n - \mathrm{GDF})^2} or UBRE
#' \eqn{D / n + 2 \phi\,\mathrm{GDF} / n - \phi}, where GDF is the Monte Carlo
#' trace of the derivative of the fitting map. Offset arrays are supported.
#'
#' @section Search:
#' 1. *Isotropic stage* (always, also with several groups): `n_grid` points
#'    on the shared \eqn{\log_{10}\lambda}, `n_refine` halvings of the grid
#'    step around the best point and one parabolic step.
#' 2. *Grouped stage* (several groups only): optional starts (`theta_start`,
#'    `n_global` Sobol points), then a coordinate pattern search from the
#'    isotropic optimum with steps 0.5 and 0.25 in \eqn{\log_{10}\lambda},
#'    extending accepted moves with doubling steps.
#' 3. *Final stage*: the isotropic and grouped optima, the nearest isotropic
#'    evaluation with a stable GDF on the smoother side (larger
#'    \eqn{\lambda}) and the next best distinct stable points (`n_final` in
#'    all, plus the smoother neighbour with one group) are re-scored with
#'    `M_final` probes.
#'
#' Comparisons are made on the scale of `score_dev`:
#' \eqn{D / \phi + 2\,\mathrm{GDF}} for UBRE (scaled deviance plus twice the
#' GDF, the AIC scale) or \eqn{n \log \mathrm{GCV}} for GCV, neither of which
#' depends on the units of the response: a move or a candidate improves
#' on another only if it lowers this score by more than `max(tol, 2 * SE)`,
#' where SE is the standard error of the paired difference over the common
#' probes. The winner is the most regular converged final candidate with a
#' stable GDF that the minimum does not improve on by this rule: isotropic
#' before grouped (`decision`), then the largest \eqn{\lambda}, in the spirit
#' of the one-standard-error rule (`tol = 0` gives a purely Monte Carlo
#' rule). On the Dette d = 8 example with Berman-Turner rows the plain argmin
#' sat on the grid boundary, where the Monte Carlo GDF was 432 \eqn{\pm} 334,
#' against 145 \eqn{\pm} 10 at the \eqn{\lambda} chosen by thinning CV.
#'
#' The winner's GDF is then re-estimated at twice its iteration budget with
#' the same probes (`budget_check`). If it moves by more than `budget_tol`,
#' the fixed-\eqn{\lambda} map is not converged there. With
#' `budget_fallback = TRUE` the candidate is rescored with the less
#' favourable of its scores at B and 2B and the choice is redone with that
#' conservative score: if it still wins it is kept, with a warning that gives
#' the GDF change and the margin, and `budget_verified = FALSE`; otherwise
#' the new winner gets its own check (at most three checks). The GDF change
#' thus matters only when it could reverse the choice: dropping such a
#' candidate outright once returned a model about 19 score units worse by
#' the exact criterion over a GDF change worth about 1 unit. An unconverged
#' or unstable winner gets one diagnostic check and no fallback, and no
#' fallback check is spent on such a candidate. The returned candidate is
#' always one whose GDF was checked: when the redone choice would return an
#' unconverged or unstable candidate (e.g. the failing winner's score at 2B
#' is not finite and no other converged candidate with a stable GDF is
#' left), or a candidate that the limit of three checks leaves unchecked,
#' the first winner is kept, with `budget_verified = FALSE` and a warning.
#'
#' @section Perturbations:
#' Probes are non-negative (\eqn{y + \epsilon (r+1)/2} against a shared
#' reference \eqn{y + \epsilon/2}), so Poisson counts are never clipped at
#' zero and each probe costs one refit. The step \eqn{\epsilon} is
#' `epsilon_rel` times the RMS of the response (Gaussian) or times
#' \eqn{\max(\mathrm{RMS}, 1)} (Poisson). Probes are drawn on cells with
#' positive weight only and are identical across \eqn{\lambda} (common random
#' numbers); the first `M_search` probes of the final bank equal the search
#' bank, so a final-stage evaluation only fits the estimator, the reference
#' and the `M_final - M_search` new probes. With `budget = "fixed"` and cold
#' probes, all fits of all points of a search step form one parallel batch.
#' With the adaptive budget (default), the pilots of the step's new points
#' run first as one batch, because they fix B, and the estimators,
#' references and probes form a second batch; warm probes also run the
#' estimators first.
#'
#' Probe validation: a probe response can need more iterations than the
#' response itself from the same cold start, so an evaluation whose probe
#' contributions contain a non-finite value or an outlier (a contribution
#' above `n_eff` in absolute value, a GDF \eqn{\le 0}, or a per-probe
#' coefficient of variation `gdf_cv` above
#' \eqn{\max(1, 3\sqrt{2 / \max(\mathrm{GDF}, 1)})}) is re-evaluated whole,
#' with the same probes, at 2B and then 4B (one more parallel batch each).
#' If it stays so, it is flagged unstable (`stable = FALSE`) and treated
#' like an unconverged one: it is never an accepted pattern-search move, the
#' isotropic optimum (unless no stable isotropic point exists) or an
#' eligible final candidate (unless nothing else is eligible, with a
#' warning). The CV limit depends on the GDF because the per-probe CV of a
#' symmetric smoother with eigenvalues in \eqn{[0, 1]} is at most
#' \eqn{\sqrt{2 / \mathrm{GDF}}}: a limit of 1 flagged heavily smoothed,
#' converged evaluations (GDF of about 4) with four probes by chance. The
#' limit is 1 from GDF 18 up.
#'
#' @section Iteration budget:
#' The estimator, reference and probe fits of one \eqn{\lambda} run a fixed
#' budget B (`control$tol` is set to `0`, so ALS / PIRLS never stop early)
#' from the same cold initialization, so that all probes differentiate the
#' same deterministic map and the deviance and the returned fit come from it
#' too. With `budget = "adaptive"` (default), `control$max_sweeps` (ALS
#' sweeps, Gaussian) or `control$pirls_maxit` (PIRLS iterations, Poisson) is
#' a cap. A pilot at \eqn{y} runs cold for 20 ALS sweeps (10 PIRLS
#' iterations), then twice as many, and so on up to the cap, until the
#' penalized objective has changed by at most `fit_tol` relative to its
#' previous value at two consecutive iterations; the first iteration where
#' that holds is `n_iter`, the evaluation is converged when it exists, and
#' \eqn{B = \min(\mathrm{cap}, \lceil 1.5\, n_{iter} \rceil + 2)} (the cap
#' when not converged), because the GDF settles more slowly than the
#' objective. The rule reads the fit history and has no floor, so B does not
#' depend on the units of the response (the stopping rule of [ttps()] tests
#' the RSS, which can stall at a turning point while the fit still moves,
#' against \eqn{\max(1, \mathrm{RSS})}). A per-iteration rule cannot see a
#' fit that drifts very slowly: the budget check below is the safeguard. A
#' fixed budget that suits one \eqn{\lambda} need not suit another: for
#' low-rank Gaussian ALS on scattered rows (n = 400, d = 2, rank 2, k = 8) a
#' 20-sweep budget gave GDF 14.4 against 21.5 to 21.9 at 80 to 320 sweeps,
#' and a GDF that was not monotone in \eqn{\lambda}. `budget = "fixed"` runs
#' every fit at the control budget and flags an evaluation as converged when
#' the relative change of its estimator's last iteration is at most `1e-5`.
#' The search log records B (`budget`), `n_iter`, `converged`, `stable` and
#' `reevals` (probe-validation rounds); the winner is chosen among converged
#' and stable final candidates, and a warning is issued when the selected
#' fit is not converged.
#'
#' @param Y d-way response array (d >= 2).
#' @param axes List of grid coordinates, one per margin (default: unit grids).
#' @param family `gaussian()` or `poisson()`.
#' @param rank,k,degree,penalty_order,cyclic,period,knots Passed to [ttps()].
#' @param weights Optional non-negative weight array (same dim as `Y`).
#' @param offset Optional offset array (same dim as `Y`, linear-predictor scale).
#' @param groups Integer vector of length d mapping margins to smoothing
#'   groups that share one \eqn{\lambda}. `NULL` = isotropic (one group);
#'   `seq_len(d)` = fully anisotropic.
#' @param theta_lower,theta_upper Search box in \eqn{\log_{10}\lambda}
#'   (scalar or one value per group; default from `control$lambda_bounds`).
#'   The isotropic stage searches the common range
#'   `[max(theta_lower), min(theta_upper)]`.
#' @param theta_start Optional matrix (or vector) of extra starting points in
#'   group space (\eqn{\log_{10}\lambda}): evaluated with the grid when there
#'   is one group, as starts of the grouped stage otherwise.
#' @param n_grid Grid size of the isotropic stage (at least 3).
#' @param n_global Sobol starts for the grouped stage (default `0`: the
#'   grouped stage starts from the isotropic optimum). Ignored with one
#'   group.
#' @param n_refine Halvings of the isotropic grid step around the best point,
#'   followed by one parabolic step. `0` = grid only.
#' @param n_final Number of candidates re-scored with `M_final` probes (with
#'   several groups, at least the isotropic and the grouped optima and the
#'   smoother isotropic neighbour; with one group the smoother neighbour is
#'   added on top of `n_final`).
#' @param M_search,M_final Rademacher probes in the search / final stages.
#' @param tol Minimum improvement of `score_dev` (default `1`; UBRE:
#'   \eqn{D / \phi + 2\,\mathrm{GDF}}, the AIC scale; GCV:
#'   \eqn{n \log \mathrm{GCV}}) for a move to be accepted, or for the minimum
#'   to beat a more regular final candidate; the paired threshold `2 * SE`
#'   applies when larger.
#' @param criterion `"auto"` (default): UBRE with known scale for Poisson
#'   and GCV for Gaussian, as `mgcv`'s `"GCV.Cp"`; or force `"gcv"` /
#'   `"ubre"`. For overdispersed counts use `"gcv"` or pass an estimated
#'   `scale` with `"ubre"`.
#' @param scale Scale for UBRE. Poisson default `1`; required for Gaussian
#'   UBRE.
#' @param probe_init `"cold"` (default) or `"warm"`.
#' @param probe_budget Iteration budget for warm probes (PIRLS iterations for
#'   Poisson, ALS sweeps for Gaussian); ignored when `probe_init = "cold"`.
#' @param epsilon_rel Finite-difference step relative to the RMS of the
#'   response on the cells (rows) with positive weight: `epsilon_rel * RMS`
#'   for Gaussian data (`epsilon_rel` when the RMS is 0), so the GDF does not
#'   depend on the units of the response; `epsilon_rel * max(RMS, 1)` for
#'   Poisson counts.
#' @param n_cores Parallel workers over the fits of a batch
#'   (`parallel::mclapply`; forked processes, so `1` on Windows).
#' @param seed Probe seed (any finite number; probe `j` uses
#'   `seed * 100003 + j`, wrapped into the integer range when it leaves it).
#' @param control [tt_control()] for every fixed-\eqn{\lambda} fit; its
#'   `max_sweeps` (Gaussian) / `pirls_maxit` (Poisson) is the cap of the
#'   adaptive budget, or the fixed budget (default caps: 400 sweeps, 60
#'   PIRLS iterations). Its `tol` is not used (see `fit_tol`).
#' @param budget_check If `TRUE` (default), re-estimate the winner's GDF with
#'   twice its iteration budget and the same probes. A relative change above
#'   `budget_tol` means the fixed-\eqn{\lambda} map is not converged at that
#'   budget (seen for low-rank Poisson fits at 10 PIRLS iterations). The
#'   doubling change is a lower bound on the distance to the converged map
#'   when ALS converges slowly (2 to 4 PIRLS iterations moved GDF by 1.5%
#'   while the converged value was 10% away), hence the tight default.
#' @param budget_tol Relative GDF change that fails the check (`0.01`).
#' @param budget_fallback If `TRUE` (default), a converged and stable winner
#'   that fails the budget check is rescored with the less favourable of its
#'   scores at B and 2B and the choice is redone; the new winner, if any, is
#'   checked in turn (at most three checks). The first winner is kept when
#'   the redone choice would return an unconverged or unstable candidate, or
#'   one left unchecked by the limit of checks. `FALSE` only warns.
#' @param budget `"adaptive"` (default) or `"fixed"` iteration budget; see
#'   section *Iteration budget*.
#' @param fit_tol Convergence tolerance of the adaptive budget (default
#'   `1e-7`): the pilot at \eqn{y} has converged at the first iteration after
#'   which the relative change of the penalized objective has been at most
#'   `fit_tol` at two consecutive iterations (ALS sweeps or PIRLS
#'   iterations; no floor, so free of the units of the response).
#' @param gdf_method `"mc"` (default): Monte Carlo GDF of the fixed-budget
#'   map. `"exact"`: one fit per point, polished to a stationary point by
#'   trust-region Newton, GDF by [tt_edf_exact()] (no probes, `gdf_se = 0`;
#'   a point whose fit does not reach a local minimum is unconverged).
#'   `"auto"`: exact for at most 2000 TT parameters and rows x parameters at
#'   most 5e7, else Monte Carlo.
#' @param verbose Print one line per computed evaluation and per stage.
#' @return A list with `lambda` (length d), `theta` (group scale), `score`
#'   (`ubre` or `gcv`), `score_dev` (the comparison scale,
#'   \eqn{D / \phi + 2\,\mathrm{GDF}} or \eqn{n \log \mathrm{GCV}}), `gcv`,
#'   `ubre`, `gdf`, `gdf_se`, `gdf_cv` (per-probe coefficient of variation of
#'   the GDF) and `stable` (the probe contributions form a usable GDF: all
#'   finite, none above `n_eff` in absolute value, GDF \eqn{> 0} and
#'   `gdf_cv` within the limit of section *Perturbations*), `deviance`,
#'   `n_eff`, `criterion`,
#'   `groups`, `decision` (`"isotropic"` or `"grouped"`), `rule` (one-line
#'   description of the final choice), `best_theta` (the minimum among the
#'   eligible final candidates), `boundary`, `budget` (check of the winner:
#'   GDF at twice the budget, its relative change, budget mode, B and
#'   `n_iter`),
#'   `budget_checks` (one row per check, with the GDF and `score_dev` at B
#'   and 2B and `ok`), `budget_verified` (`TRUE` when the returned winner
#'   passed its check, `NA` without a check), `paired` (final candidates with
#'   `score_dev`, `score_dev_used` (the score in the final choice: the less
#'   favourable of B and 2B after a failed check), `gdf_cv`, `stable`,
#'   `converged`, the paired difference and SE against the minimum,
#'   `qualifies`, `best`, `winner`, `budget_ok` and the difference against the
#'   winner), `fit` (the [ttps()] fit at the winner, with `$ggcv`), `search`
#'   (one row per computed evaluation, with `gdf_cv`, `stable`, convergence
#'   flags, the budget B, `n_iter`, `reevals` (probe-validation rounds) and
#'   the number of fits run), `n_evals` (computed evaluations, budget checks
#'   excluded), `n_fits` (fits run, budget checks included), `method` and
#'   `elapsed`.
#' @seealso [tt_gdf_array()] to score a given \eqn{\lambda}; [tt_ggcv()] for
#'   scattered rows.
#' @export
tt_ggcv_array <- function(Y,
                          axes = NULL,
                          family = stats::gaussian(),
                          rank = 3L,
                          k = 8L,
                          degree = 3L,
                          penalty_order = 2L,
                          cyclic = NULL,
                          period = NULL,
                          knots = NULL,
                          weights = NULL,
                          offset = NULL,
                          groups = NULL,
                          theta_lower = NULL,
                          theta_upper = NULL,
                          theta_start = NULL,
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
                          verbose = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  criterion <- match.arg(criterion)
  probe_init <- match.arg(probe_init)
  budget <- match.arg(budget)
  map <- .ggcv_map_array(
    Y = Y, axes = axes, family = family, rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, cyclic = cyclic, period = period,
    knots = knots, weights = weights, offset = offset
  )
  # the removed working-response proxy, as in tt_ggcv() (the ttps() array
  # route warns itself and passes the knob on as "algorithmic")
  if (identical(map$key, "poisson") && identical(control$ggcv_glm_mode, "working")) {
    .ggcv_warn_working("tt_ggcv_array()")
  }
  st <- .ggcv_setup(
    map, groups = groups, criterion = criterion, scale = scale,
    probe_init = probe_init, probe_budget = probe_budget,
    epsilon_rel = epsilon_rel, n_cores = n_cores, seed = seed,
    control = control, budget = budget, fit_tol = fit_tol,
    gdf_method = match.arg(gdf_method)
  )
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  out <- .ggcv_search(
    st, lo = theta_lower %||% log10(bounds[1L]),
    hi = theta_upper %||% log10(bounds[2L]), n_grid = n_grid,
    n_refine = n_refine, n_global = n_global, n_final = n_final,
    M_search = M_search, M_final = M_final, tol = tol,
    theta_start = theta_start, budget_check = budget_check,
    budget_tol = budget_tol, budget_fallback = budget_fallback,
    verbose = verbose
  )
  out$elapsed <- proc.time()[["elapsed"]] - t0
  out
}

#' Global GDF and gGCV score of an array TT fit at a fixed lambda
#'
#' Scores one smoothing vector with the same estimator as [tt_ggcv_array()]:
#' the fixed-\eqn{\lambda} array fit (K-ALS or weighted PIRLS), its deviance,
#' and a Monte Carlo GDF from non-negative perturbations (see the sections
#' *Perturbations* and *Iteration budget* of [tt_ggcv_array()]). [tt_gdf()]
#' is the scattered-row counterpart.
#'
#' @inheritParams tt_ggcv_array
#' @param lambda Smoothing vector (scalar or length d, non-negative).
#' @param M Number of probes (ignored when `probes` is a matrix).
#' @param probes `"rademacher"` (default; probe `j` drawn from
#'   `seed * 100003 + j`), `"unit"` (one perturbation per weighted cell: the
#'   exact finite-difference trace, one fit per cell, small arrays only), or a
#'   numeric matrix of -1 / +1 values with one row per cell (in
#'   `as.numeric(Y)` order) and one column per probe, a user probe bank
#'   (`M = ncol(probes)`; rows of zero-weight cells are ignored).
#' @param control [tt_control()] for every fit: `max_sweeps` (Gaussian) /
#'   `pirls_maxit` (Poisson) is the cap of the adaptive budget or the fixed
#'   budget (default caps: 400 sweeps, 60 PIRLS iterations). The defaults
#'   changed: up to package commit 77d5fc6 `tt_gdf_array()` ran a fixed
#'   budget of 10 sweeps / 20 PIRLS iterations; `budget = "fixed"` with
#'   `control = tt_control(max_sweeps = 10, pirls_maxit = 20,
#'   compute_edf = FALSE, seed = 1)` gives that map again.
#' @return A list with `gdf`, `gdf_se`, `deviance`, `gcv`, `ubre`, `n_eff`,
#'   `contrib` (per-probe contributions; per cell for `"unit"`), `M_ok`,
#'   `converged` and `last_rel_change` (convergence of the estimator),
#'   `budget` (iteration budget B of the estimator, reference and probe
#'   fits), `n_iter` (adaptive: the iteration where the pilot's objective
#'   settled, the cap when it did not; fixed: iterations of the estimator),
#'   `gdf_cv` (per-probe coefficient of variation; `0` for `"unit"`),
#'   `stable` (see [tt_ggcv_array()]), `reevals` (probe-validation rounds,
#'   each at twice the budget), `fit` (the estimator), `epsilon`, `n_fits`
#'   and `time_s`.
#' @seealso [tt_gdf()] for scattered rows.
#' @export
tt_gdf_array <- function(Y,
                         lambda,
                         axes = NULL,
                         family = stats::gaussian(),
                         rank = 3L,
                         k = 8L,
                         degree = 3L,
                         penalty_order = 2L,
                         cyclic = NULL,
                         period = NULL,
                         knots = NULL,
                         weights = NULL,
                         offset = NULL,
                         M = 16L,
                         probes = "rademacher",
                         scale = NULL,
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
                         gdf_method = c("mc", "exact")) {
  probe_init <- match.arg(probe_init)
  budget <- match.arg(budget)
  gdf_method <- match.arg(gdf_method)
  map <- .ggcv_map_array(
    Y = Y, axes = axes, family = family, rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, cyclic = cyclic, period = period,
    knots = knots, weights = weights, offset = offset
  )
  .ggcv_gdf_point(
    map, lambda = lambda, M = M, probes = probes, scale = scale,
    probe_init = probe_init, probe_budget = probe_budget,
    epsilon_rel = epsilon_rel, n_cores = n_cores, seed = seed,
    budget = budget, fit_tol = fit_tol, control = control,
    gdf_method = gdf_method
  )
}

# ---------------------------------------------------------------------------
# ttps() dispatch
# ---------------------------------------------------------------------------

#' Resolve `lambda = "gGCV"` inside [ttps()] when `array = TRUE`.
#'
#' Search knobs are read from `control`: first the array-specific names
#' `ggcv_array_<knob>` (set on a [tt_control()] object with `$<-`), then the
#' unified `ggcv_<knob>`, then the defaults of [tt_ggcv_array()]. Knobs:
#' `groups`, `n_grid`, `n_global`, `n_refine`, `n_final`, `M_search`,
#' `M_final`, `tol`, `criterion`, `probe_init`, `n_cores`, `budget_check`,
#' `budget_fallback`, `budget`, `fit_tol`, `refit`, and the iteration cap
#' (adaptive budget) or fixed budget of every selection fit, `max_sweeps` /
#' `pirls_maxit` (`ggcv_max_sweeps` / `ggcv_pirls_maxit`). Poisson groups
#' follow `.ggcv_poisson_groups()` (the groups knob, then
#' `ggcv_poisson_anisotropic`), as on scattered rows. `refit(lambda)` refits
#' the caller's model at the selected lambda (see `.ttps_ggcv_finish()`);
#' without it the selection fit is returned.
#' @keywords internal
#' @noRd
.ttps_dispatch_ggcv_array <- function(Y, axes, family, rank, k, degree,
                                      penalty_order, cyclic, period, knots,
                                      weights, offset, control, cl = NULL,
                                      refit = NULL) {
  w <- if (is.null(weights) || all(abs(as.numeric(weights) - 1) < 1e-12)) {
    NULL
  } else {
    weights
  }
  off <- if (is.null(offset) || all(as.numeric(offset) == 0)) NULL else offset
  kb <- .ggcv_control_knobs(control, prefix = "ggcv_array_")
  # Poisson: ggcv_array_groups / ggcv_groups, then ggcv_poisson_anisotropic,
  # as on scattered rows and in tt_ggcv_poisson()
  is_pois <- identical(family_key(normalize_family(family)), "poisson")
  groups <- if (is_pois) {
    .ggcv_poisson_groups(control, length(dim(Y)), groups = kb$groups)
  } else {
    kb$groups
  }
  ctrl <- control
  # the removed working-response proxy: warned once here, as on scattered
  # rows (the array route never used it, so nothing else changes)
  if (is_pois && identical(control$ggcv_glm_mode, "working")) {
    .ggcv_warn_working("ttps(array = TRUE, lambda = \"gGCV\")")
    ctrl$ggcv_glm_mode <- "algorithmic"  # warned once; tt_ggcv_array() stays quiet
  }
  ctrl$max_sweeps <- kb$max_sweeps
  ctrl$pirls_maxit <- kb$pirls_maxit
  opt <- tt_ggcv_array(
    Y = Y, axes = axes, family = family, rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, cyclic = cyclic, period = period,
    knots = knots, weights = w, offset = off, groups = groups,
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
