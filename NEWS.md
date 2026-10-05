# TTPsplines NEWS

## Development (0.0.0.9001)

### gGCV: adaptive pilot on the objective, probe validation, decision-aware budget check

* Adaptive budget, new pilot rule: a pilot at y runs cold with `tol = 0`
  for 20 ALS sweeps (10 PIRLS iterations), then twice as many, up to the
  cap, until the relative change of the penalized objective, read from the
  fit history without a floor, has been at most `fit_tol` at two
  consecutive iterations. That iteration is `n_iter`; the evaluation is
  converged when it exists, and B = min(cap, ceiling(1.5 n_iter) + 2) (the
  cap otherwise). The estimator (deviance, returned fit, `tt_gdf()$fit`,
  `ggcv_refit = FALSE`) is now the cold fit at B with `tol = 0`, the same
  map as the probes, not the tolerance-stopped run. The previous pilot
  stopped where `ttps()` stops, on the RSS against max(1, RSS): on a
  low-rank fit (n = 400, rank 3, k = 10, lambda = 1e-3) it flagged
  convergence after 13 sweeps at a turning point of the RSS while the
  objective still changed by 2.5e-5 per sweep (after 3 sweeps for y in
  units of 1e-3), and B = 15 gave GDF 54 against 84 at 2000 sweeps; the new
  rule settles at sweep 69 in any units (B = 106). A per-iteration rule
  cannot see a very slow drift; the budget check remains the safeguard.
* Probe validation: an evaluation whose probe contributions contain a
  non-finite value, a value above `n_eff` in absolute value, a GDF <= 0 or a
  per-probe coefficient of variation above max(1, 3 sqrt(2 / max(GDF, 1)))
  is re-evaluated with the same probes at 2B and then 4B. The CV limit is 1
  from GDF 18 up and larger below: the per-probe CV of a symmetric smoother
  with eigenvalues between 0 and 1 is at most sqrt(2 / GDF), and a limit of
  1 flagged converged, heavily smoothed evaluations (GDF 3.4 to 4.7, four
  probes) whose contributions equalled the exact ones (11 of 410 search
  rows on smooth truths at rank = k = 6; none with the new limit). If it
  stays so, it is flagged unstable
  (`stable = FALSE`) and is never an accepted pattern move, the isotropic
  optimum (unless no stable isotropic point exists) or an eligible final
  candidate (unless nothing else is eligible, with a warning); the smoother
  and next-best final candidates are stable search points. A fixed
  60-sweep rank-2 fit gave contributions (-0.6, -18.9, 53.8, 24.4), and
  (6.2, 10.8, 12.1, 6.2) at 120 sweeps. `stable` now covers all these
  criteria (it was `gdf_cv <= 1`); the search log gains `reevals`, and
  `tt_gdf()` / `tt_gdf_array()` return `stable` and `reevals`.
* Decision-aware budget check, replacing the rule that dropped a failing
  winner: when the winner's GDF moves by more than `budget_tol` at twice
  its budget, the candidate is rescored with the less favourable of its
  scores at B and 2B and the choice is redone with that conservative score.
  It is kept if it still wins (`budget_verified = FALSE`, with a warning
  that gives the GDF change and the margin); otherwise the new winner gets
  its own check (at most three checks). On the anisotropic test data
  (seed 61) at fixed budgets of 80 to 400 sweeps, before the core-solve
  change below, the old rule dropped the grouped optimum over a GDF change
  worth about 1 score unit, although it led the returned isotropic
  candidate by about 20 units (19 by the exact criterion); a test now forces
  that GDF change (-5.6%) and checks that the grouped optimum is kept. The
  first winner is always checked, but no fallback check is
  spent on an unconverged or unstable candidate: when no final candidate
  converged, one diagnostic check runs instead of three (7 x 6 x 6 x 5
  Poisson array: 18 instead of 54 fits in budget checks, same selection).
  The returned candidate is always a checked one: when the redone choice
  would return an unconverged or unstable candidate (the failing winner's
  score at 2B is not finite and no other converged, stable candidate is
  left; a mocked case returned an unconverged, unchecked candidate 4.5
  score units worse), or a candidate that the limit of three checks leaves
  unchecked, the first winner is kept (`budget_verified = FALSE`, with a
  warning). A kept winner's warning names a runner-up only when the rule
  could choose it (converged, stable GDF) and otherwise calls the winner
  the only eligible final candidate; a GDF change that is not finite reads
  "not finite" (it printed "NA%").
  `budget_checks` gains `score_dev` and `score_dev_2x`; `paired` gains
  `score_dev_used`.
* Fit errors are kept: `tt_gdf()` and `tt_gdf_array()` stop with the fit's
  message when the estimator fails (they returned NA silently), failed
  perturbation fits give a warning, and selections that end with no finite
  criterion report the first fit error (e.g. "`period` must be a list of
  length d").
* Poisson smoothing groups resolve in one order in `tt_ggcv_poisson()` and
  in `ttps(lambda = "gGCV")` on scattered rows and on arrays: the
  `anisotropic` argument, then `ggcv_groups` (array mode:
  `ggcv_array_groups` first), then `ggcv_poisson_anisotropic`.
  `tt_ggcv_poisson()` read the last two in the opposite order, and the
  array route ignored `ggcv_poisson_anisotropic`.
* `tt_ggcv(family = poisson())`, `tt_ggcv_array(family = poisson())` and
  `ttps(array = TRUE, family = poisson(), lambda = "gGCV")` warn on
  `ggcv_glm_mode = "working"` like the other routes (once per call). The
  array route never used the working proxy, so its selections do not
  change.
* The internal RSS-and-objective stop that the old pilot added to
  `tt_als_fit_sequential()` during this cycle is removed again.

### gGCV numerics: core solve, probe seeds, UBRE scale, probe step

* Global-mode TT core systems (every fixed-lambda ALS / PIRLS-ALS core
  update of `ttps()`) are solved by one symmetric eigendecomposition, the
  minimum-norm solve that drops eigenvalues below 1e-12 times the largest,
  instead of Cholesky, then LU, then a ridge. At saturating rank under heavy
  smoothing these systems are singular along gauge directions (reciprocal
  condition down to 1e-21); the old cascade switched with roundoff, so the
  fixed-lambda map jumped under tiny changes of y and Monte Carlo GDFs were
  off by up to 1e4 (Exp B smoke mesh, rank 10: 6 of 81 cells off by more
  than 1 df). On that mesh the paired error against the dense smoother is
  now at most 0.035 df. Fitted values of well-conditioned fits change by
  less than 1e-8 (relative to sd(y)). Speed cost: the eigendecomposition
  costs about 9 Cholesky factorizations, and every fixed-lambda global-mode
  fit of `ttps()` pays it, so every fit of the gGCV engine does. On four
  benchmarks of fixed-lambda `ttps()` fits against the old solve (Gaussian
  d = 3, rank 6 and 10, k = 10; Gaussian d = 4, rank 6, k = 8; Poisson
  d = 3, rank 8, k = 10; n = 2000 to 3000) a fit took 1.3 to 1.8 times as
  long; the largest cores (rank 10) slowed down most.
* Probe seeds: probe `j` uses `seed * 100003 + j` in double precision,
  wrapped into the integer range only when it leaves it. Seeds above 21474
  (e.g. date-like `control$seed`) no longer overflow and stop gGCV, and
  every seed that worked before keeps its probes. The probes of a batch are
  drawn once in the parent process, so a failed draw is an error, not NA
  contributions from forked workers; `seed` must be one finite number.
* UBRE comparisons use the AIC scale, `score_dev = D / scale + 2 * GDF`, so
  `tol` and the paired SE mean the same for any units of the response. With
  `D + 2 * scale * GDF`, `tol` was in units of y^2 and Gaussian UBRE
  selections changed with the units of y. Poisson (scale 1) is unchanged;
  the `score` / `ubre` fields keep their definition.
* Gaussian probe step `epsilon_rel * RMS(y)` without the floor at 1
  (`epsilon_rel` when y = 0), so the GDF does not depend on the units of y
  (for a response in units of 1e-3 the old step was 1.5 x RMS and overstated
  a low-rank GDF by up to 3 df). Poisson keeps `epsilon_rel * max(RMS, 1)`.

### gGCV: noise-aware final choice, adaptive iteration budget, `tt_gdf()`

* Final choice: among the converged final candidates, the most regular one
  that the minimum does not beat by more than `max(tol, 2 * paired SE)`
  wins (isotropic before grouped, then the largest lambda), in the spirit of
  the one-standard-error rule; this replaces the plain argmin with an
  isotropic preference. On the Dette d = 8 example with Berman-Turner rows
  the plain argmin sat on the grid boundary, where the Monte Carlo GDF was
  432 +- 334 (145 +- 10 at the lambda chosen by thinning CV). The final
  stage adds the nearest isotropic evaluation on the smoother side; the
  value gains `rule` and `best_theta`, and `paired` gains `gdf_cv`,
  `stable`, `diff_vs_best`, `se_vs_best`, `qualifies`, `best`, `winner` and
  `budget_ok`.
* Budget check of the winner (`budget_check`, `budget_fallback = TRUE`,
  knob `ggcv_budget_fallback`): the winner's GDF is re-estimated at twice
  its iteration budget, with the decision-aware fallback described above.
  New fields `budget_checks` and `budget_verified`.
* Adaptive iteration budget (`budget = "adaptive"`, the default; knob
  `ggcv_budget`; tolerance `fit_tol`, default `1e-7`, knob `ggcv_fit_tol`):
  `pirls_maxit` (Poisson) / `max_sweeps` (Gaussian) is a cap, and B follows
  lambda (pilot rule above). At a fixed 20-sweep budget the low-rank
  Gaussian ALS GDF was 14.4 against 21.5 to 21.9 at 80 to 320 sweeps. The
  search log gains `budget`, `n_iter`, `gdf_cv` and `stable`.
* `budget = "fixed"` runs every fit at the control budget. With the control
  of commit 77d5fc6 (`tt_control(max_sweeps = 10, pirls_maxit = 20)`) it
  reproduces the fixed-lambda GDF evaluations of `tt_gdf_array()` at that
  commit up to roundoff (about 1e-12 relative on two test arrays; the core
  solve changed, see above), provided the probe step is the same (Poisson;
  Gaussian responses with RMS(y) >= 1) and the probe validation does not
  re-evaluate the point. It does not reproduce selections: the search, the
  smoother final candidate, the noise-aware rule, the converged-only
  eligibility and the budget check all changed, so `tt_ggcv_array()`,
  `tt_ggcv()`, `tt_ggcv_poisson()` and `ttps(lambda = "gGCV")` select
  differently in either budget mode (lab callers such as `dette_ggcv.R`,
  `11_ggcv_fullgrid.R` and `08_lambda_selectors.R` included). Old
  selections need the package at commit 77d5fc6.
* New defaults (caps): `tt_ggcv()`, `tt_ggcv_array()`, `tt_ggcv_poisson()`,
  `tt_gdf_array()` and the new `tt_gdf()` use
  `tt_control(max_sweeps = 400, pirls_maxit = 60)`; `tt_control()` knobs
  `ggcv_max_sweeps` 20 -> 400 and `ggcv_pirls_maxit` 25 -> 60. This changes
  the default results of `tt_gdf_array()`, which ran a fixed budget of 10
  sweeps / 20 PIRLS iterations: `budget = "fixed"` with
  `control = tt_control(max_sweeps = 10, pirls_maxit = 20,
  compute_edf = FALSE, seed = 1)` gives the old map (see the previous
  item). With `ggcv_refit = FALSE`, `ttps(lambda = "gGCV")` returns the
  selection's estimator (the cold fit at the winner's budget B, `tol = 0`).
* New exported `tt_gdf()`: GDF, deviance, GCV and UBRE of a scattered-row
  fit at a fixed lambda (counterpart of `tt_gdf_array()`). Both accept
  `probes = <matrix>`, a user bank of -1 / +1 probes (one row per
  observation or cell, one column per probe).
* Warning when the Monte Carlo GDF at the selected lambda is unstable
  (`stable = FALSE`, see the probe validation above).
* `inst/benchmarks/ggcv_array/validate_gdf.R` and
  `check_tedf_saturating.R` pin `budget = "fixed"`, the budgets they study.

### One gGCV engine for scattered and array data

* `tt_ggcv()` (Gaussian and Poisson scattered rows, offsets allowed),
  `tt_ggcv_poisson()` and `tt_ggcv_array()` share one evaluator and one
  search (`R/ggcv_engine.R`): non-negative probes against a shared
  reference, common random numbers, an iteration budget per lambda (see
  above) from a cold common initialization; isotropic grid with halvings
  and a parabolic step, optional grouped pattern search, final candidates
  compared with the paired Monte Carlo SE (`tol`, default 1 score unit;
  final choice above). `tt_ggcv()` no longer calls
  the lab optimizer `tt_global_lambda_optimize()` (since removed, see
  "Lab oracle stack removed" below). `tt_ggcv()` gains
  `offset`, `cyclic`, `period`, `knots`, `groups`, `criterion`, `scale`,
  `tol`, `n_cores` and the budget check.
* `ttps(lambda = "gGCV")` selects with the `ggcv_*` knobs of `tt_control()`
  and the iteration budget `ggcv_budget` (caps `ggcv_max_sweeps` /
  `ggcv_pirls_maxit`), then refits
  at the selected lambda with the caller's own arguments and control
  (`ggcv_refit = TRUE`); `fit$ggcv` keeps the selection diagnostics. The
  selection now uses the call's `rank` as given (not `max(rank)`) and, for
  Gaussian data, its `offset`, `cyclic`, `period` and `knots`.
* New `tt_control()` knobs: `ggcv_groups`, `ggcv_n_grid`, `ggcv_n_final`,
  `ggcv_tol`, `ggcv_criterion`, `ggcv_probe_init`, `ggcv_n_cores`,
  `ggcv_budget_check`, `ggcv_max_sweeps`, `ggcv_pirls_maxit`, `ggcv_refit`.
  Changed defaults: `ggcv_n_refine` 5 -> 2, `ggcv_M_search` 15 -> 4,
  `ggcv_M_final` 40 -> 16, `ggcv_glm_mode` `"working"` -> `"algorithmic"`.
* Deprecated: the Poisson working-response proxy (`mode = "working"` /
  `ggcv_glm_mode = "working"`) warns and is ignored; it row-scaled the
  covariates instead of the basis rows.
* Removed old gGCV steps that no selection path uses:
  * the cGCV anchor: argument `include_cgcv_anchor` of `tt_ggcv()` and
    `tt_ggcv_poisson()` and knob `ggcv_include_cgcv_anchor` of
    `tt_control()`. Extra starts remain available (`extra_theta`,
    `theta_start`, `n_global`); a control object that still carries the
    field (set with `$<-`) runs unchanged and the field is ignored.
  * the lab-only arguments of `tt_ggcv()` (`n_diverse`,
    `core_starts_search`, `core_starts_final`, `min_dist`, `boundary_tol`,
    `scheme`, `fit_backend`, `adaptive_fidelity`, `fidelity`, `gdf_init`);
    they now fail as unused arguments.
  * the old algorithmic Poisson evaluator (`.tt_pois_fit_fixed()`,
    `.tt_pois_ggcv_eval()`, probes clipped at 0). Archival scripts that
    call it reproduce with the package at commit 77d5fc6.
* `tt_ggcv_array()`: `n_global` default 0 (was `8 * G` Sobol points),
  `n_refine` counts isotropic halvings, new `tol`; the value gains
  `score_dev`, `decision`, `paired`, `n_evals` and per-evaluation
  convergence flags in `search`.

### Lab oracle stack removed

* Removed the internal (never exported) lab gGCV oracle: `tt_global_gcv()`,
  `tt_global_gdf_mc()`, `tt_global_lambda_optimize()`,
  `tt_global_lambda_optimize_v1()` and their `.tt_lab_*` helpers (files
  `R/global_gcv_lab.R`, `R/global_lambda_optimize.R`,
  `R/global_lambda_optimize_v1.R`), their 42 tests in 7 files and the
  benchmark directory `inst/benchmarks/global_gcv/`. gGCV selection runs
  only on the engine (`tt_ggcv()`, `tt_ggcv_poisson()`, `tt_ggcv_array()`,
  `ttps(lambda = "gGCV")`). Scripts that call the oracle (for example
  `tt_global_gdf_mc()` for a Monte Carlo GDF) reproduce with the package at
  commit 77d5fc6.
* Helpers still in use were moved first, code unchanged:
  `.tt_with_preserved_seed()`, `.tt_clone_cores()` and
  `.tt_lab_rademacher_probes()` (same name; lab scripts call it through
  `TTPsplines:::`) now live in `R/utils_internal.R`; the Sobol generator of
  the optional `n_global` starts moved into the engine as
  `.ggcv_sobol_params()`, `.ggcv_sobol_unit()` and `.ggcv_sobol_box()`
  (same points).
* Kept: the fixed-lambda ALS building blocks under the global penalty
  (`tt_als_core_update_global()`, `tt_als_sweep_global()`,
  `tt_als_fit_fixed_global()`, C++ `tt_als_fit_fixed_global_cpp`) and their
  tests.

### gGCV for array data (`array = TRUE`)

* New exported `tt_ggcv_array()` (smoothing selection) and `tt_gdf_array()`
  (GDF and score at a fixed lambda); `ttps(array = TRUE, lambda = "gGCV")`
  dispatches there (knobs `control$ggcv_array_*`).
* GDF by Monte Carlo finite differences with non-negative perturbations
  against a shared reference (no clipping of Poisson counts at 0), fixed
  iteration budget (now adaptive by default, see the first section) and
  common cold initialization; optional exact unit-vector trace for small
  arrays.
* `criterion = "auto"` (default): UBRE with scale 1 for Poisson, GCV for
  Gaussian (as `mgcv` `"GCV.Cp"`).
* `budget_check = TRUE`: re-estimates the winner's GDF at twice the iteration
  budget and warns if it moves by more than 1% (now with a fallback to the
  next candidate, see the first section).
* Validation: `inst/benchmarks/ggcv_array/validate_gdf.R` (low rank, d = 3 and 5).

### Fixes

* `glam_fit_poisson(fit_weights = )`: the intercept update and the deviance now
  use the fit weights. Before, cells with weight 0 (e.g. held-out test cells)
  still entered the intercept, so the fit was not the weighted MLE.
* `ttps()` no longer reseeds the caller's global RNG stream; the core
  initialization is unchanged for a given `control$seed`. Scripts that drew
  random numbers after a `ttps()` call without their own `set.seed()` now get
  different draws than before.

### Removed Adam/Keras stub

* Dropped unused `optimizer = "Adam"` / `backend = "keras"` stub, `tt_has_keras()`,
  `tt_keras_status()`, `reticulate` Suggests, and the optional Keras benchmark script.
  Direct-likelihood path remains ALS / PIRLS-ALS / LBFGS / GD (+ experimental hybrids).

### Information criteria

* New exported `tt_ic(fit, "AIC"|"BIC")`: working AIC/BIC from joint linearized `fit$edf` (+1 for intercept).

### Array input mode (`array = TRUE`)

* `ttps(Y, axes = list(...), rank = r, k = k, lambda = ..., array = TRUE)`
  accepts a d-way array `Y` (dim1 fastest, as from `array(y, dim = c(n1,...,nd))`)
  plus a list of marginal coordinate vectors `axes`.
* When `array = TRUE` the conditional Gram and RHS for each TT core are computed
  **without materialising** the n × q_k design matrix, using the Kronecker
  structure of the complete grid:
    - `S_k = kron(R'R, kron(Bk'Bk, L'L))` (Gaussian, unweighted)
    - `b_k` via triple-mode contraction of `Y_centered` over `(L_uniq, Bk, R_uniq)`
* `axes` can be omitted; defaults to a unit-interval grid for each margin.
* Gaussian unweighted grids use the Kronecker Gram path; Poisson / weighted
  grids use the weighted array Gram (still without forming \(X_k\)).
* Restrictions: no `linear=` / `smooth=`; `null_space = "joint"` only;
  no `lambda = "CV"` / `"gGCV"`. Walkthrough: `vignette("array-mode")`.
* Result is numerically identical to `array = FALSE` on the same grid data
  (fitted values differ by ≤ machine epsilon in all tests).
* New internal functions: `tt_gram_rhs_array()`, `.tt_array_extract_interfaces()`,
  `.tt_array_rhs()` in `R/tt_gram_array.R`.

### Null-space modes (TODO-SC-NULL)

* `ttps(..., null_space = )`:
  - `"joint"` (default) — single TT on full \(\Theta\).
  - `"profiled"` — **experimental** Gaussian NSP: profile \(\beta_0\);
    ALS on \(Q_0 y\) and \(Q_0 Z_k\); GDF \(= \mathrm{rank}(X_0)+\mathrm{GDF}_{TT,\perp}\).
    No GLM / cyclic / `linear` / `smooth` yet.
* Removed: `"sequential"` / `"separate"` (fixed-offset OLS→TT; not joint NSP).
* Cap: `tt_control(null_space_max_npar = ...)`.
* Gate: `inst/benchmarks/null_space/GATE_PROFILED.md`.

### EDF API (joint + per margin)

* `fit$edf_margin` = diagonal blocks \(\operatorname{tr}(H_{kk})\) of the joint
  influence \(H=(J'J+P)^{-1}J'J\); **sums to** `fit$edf` (TT parameter-block
  partition).
* `fit$edf_margin_cond` = conditional ALS core EDFs (cGCV diagnostic; not
  additive).
* `tt_edf(fit)` returns joint + both margin summaries.

### Margin drop test

* New `tt_margin_drop_test()`: leave-one-margin-out deviance comparison with
  approximate F (`method = "nested"`) or permutation (`method = "permute"`)
  p-values; soft `drop_candidates` when `p > alpha`.

### Margin Activity Path

* New helper `tt_margin_activity_path()`: screen TT margins by partial-range
  activity along an isotropic \(\lambda\) path, rank covariates, and choose a
  nested top-\(m\) subset by \(K\)-fold CV (`select = "1se"` default, or
  `"min"` / `"none"`). Refits [ttps()] with `lambda = "cGCV"` on the selected
  columns. S3 `print` / `plot` methods.
* Vignette: `vignette("margin-activity-path", package = "TTPsplines")`.
* Protocol / extended demo:
  `inst/benchmarks/margin_path/PROTOCOL_MAP.md`,
  `inst/benchmarks/margin_path/run_map_example.R`.
* Prefer the full name **Margin Activity Path** (avoid acronym "MAP").

## Development (0.0.0.9000)

### TT-DLNM (distributed lag)

* New API `ttps_dlnm(y, list(temp=, pm10=, …), lag = L, …)`: joint exposure×lag
  P-spline surface in TT form,
  \(\eta_t=\ldots+\sum_{\ell=0}^{L}h(x_{t-\ell},\ell)\), **without** building the
  dense cross-basis \(W\). ALS/PIRLS uses summed conditional designs over lags.
* Calendar confounding via existing `linear=` / `smooth=` (recommended:
  `s(year)+s(month, bs="cc")` + DOW).
* `predict_dlnm(fit, var=, at=, cen=, type="overall"|"slice")` for
  Gasparrini-style overall (lag basis summed) and fixed-lag slices.
* Tests: `tests/testthat/test-dlnm.R`.

### Additive smooths + parametric linear

* `ttps(..., smooth = )` — additive 1D P-splines jointly with the TT surface:
  `bs = "ps"` (open) or `"cc"` (circular), basis size `k`, penalty order `m`
  (alias `penalty_order`), and per-term `lambda` (numeric or `"cGCV"`;
  default `lambda_smooth = "cGCV"`) **or** `target_edf` (root-find \(\lambda\)
  so \(\mathrm{edf}(\lambda)\approx\) target; useful for epi time trends).
  Example:
  `smooth = list(time = list(x = d$time, bs = "ps", k = 80, m = 2, target_edf = 40))`.
* Estimated by backfitting inside ALS / PIRLS-ALS; smooth \(\lambda\) via
  conditional cGCV, target EDF, or fixed (same spectral helpers as TT cores).
* `summary()` prints a **Smooth terms** table (`bs`, `k`, `m`, `edf`,
  `target_edf`, `lambda`, method) plus the glm-style parametric coefficient
  table for `linear=`.
* `predict(..., se.fit=TRUE)` now works with `linear=` / `smooth=` treating
  those terms as a fixed offset (TT Level-1 SE only). Optional
  `contrast_row=` gives SE of a link contrast for centered RR curves.
* Threaded through `predict` (matching `smooth=` / `linear=` newdata),
  `tt_rank_select` / `tt_rank_refit`, `ttps_multistart`.
* Unsupported (error): LBFGS / GD / hybrid / DN-ALS / LBFGS-ALS.

### Methodological

* **DECISION:** the package always uses the classical global discrete
  P-spline penalty on \(\Theta\),
  \(J_{\boldsymbol\lambda}(\Theta)=\sum_m\lambda_m\|\Theta\times_m\Delta\|_F^2\),
  via the exact conditional restriction
  \(P_k^{\mathrm{full}}=A_k^\top S_{\boldsymbol\lambda} A_k\).
  The former own-margin / separable surrogate is **removed** (not a
  classical multidimensional P-spline criterion).
* Fixed-\(\lambda\) Gaussian ALS records a per-core \(Q\) non-increase
  diagnostic (`fit$q_descent`).
* cGCV searches \(\lambda_k\) with fixed cross-margin offset \(P_{k,-k}\).
* L-BFGS / GD / PIRLS / Damped-Newton-ALS / LBFGS-ALS use the same global
  \(J_{\boldsymbol\lambda}\).
* Rcpp helpers: `tt_conditional_penalty_full_cpp`,
  `tt_global_penalty_value_cpp` (ALS/PIRLS sweeps remain in R).

### Diagnostics

* New `ttps_multistart()`: several random TT inits, best fit by penalized
  objective (or deviance), start table, and per-margin boundary fractions.
  Use when cGCV λ hits search bounds or low-\(r\) ALS looks init-sensitive.
  `summary()` of cGCV fits with boundary hits points here.
* Rank CV multi-start remains `tt_rank_select(..., n_starts = ...)`.

### Bases / prediction

* Open B-spline knots for covariates in \([0,1]\) now span the **unit
  interval** (same heuristic as cyclic margins), so `predict()` on
  `seq(0,1)` no longer collapses to the intercept outside `range(X)`.
* Soft warning when `newdata` falls outside a non-cyclic knot span.

### cGCV dynamics (experimental)

* `tt_control(cgcv_update=)`: **`"outer_simultaneous"` (default)** —
  fit all cores at fixed \(\lambda\) → freeze → Jacobi proposals →
  damped / trust-region update — or `"sequential"` (legacy Gauss–Seidel;
  can oversmooth-cascade on Chicago Poisson).
* Defaults after Chicago validation: `cgcv_damping = 0.25`,
  `cgcv_max_log10_step = 1`. Also `"scale_anisotropy"` parameterization
  and `tt_cgcv_frozen_curves()`.
* Fit objects store `fit$cgcv` with proposals / traces. Global
  \(J_{\boldsymbol\lambda}\) unchanged; own-margin not restored.
