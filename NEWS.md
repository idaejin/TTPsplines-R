# TTPsplines NEWS

## Development (0.0.0.9000)

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
