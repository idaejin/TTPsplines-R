# TTPsplines NEWS

## Development (0.0.0.9001)

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
* Unsupported (error): LBFGS / GD / hybrid / Adam / DN-ALS / LBFGS-ALS.

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
