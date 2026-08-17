# Protocol: fixed-λ ALS C++ path (global \(P_k^{\mathrm{full}}\))

**Package:** `TTPsplines` (`ttpsplines-pkg`)  
**Updated:** 2026-08-17

## Architecture (DECISION)

```text
outer λ (R: Sobol / BO / global-GCV / cGCV)
  +
fixed-λ ALS (C++, global P_k^full, parity vs R)
```

Do **not** port outer cGCV / global-GCV search to C++ first.

## P0 — Documentation — **PASS** (2026-08-17)

- README / `?tt_control` / `?ttps` / scalability vignette: ALS/PIRLS sweeps stay in **R**.
- `backend = "Rcpp"` warns; kernels only (Gram / \(P^{\mathrm{full}}\)).
- `tt_cgcv_fit_cpp` / `tt_glm_pirls_cgcv_cpp` marked **legacy own-margin**.
- Commit: `b27e12a`.

## P1 — Single core update — **PASS** (2026-08-17)

Gate: `tests/testthat/test-als-core-update-global.R` — **352 PASS** (R vs Rcpp).

### Scope

| In | Out |
|----|-----|
| Gaussian | GLM / PIRLS |
| Numeric fixed \(\boldsymbol\lambda\) | cGCV |
| `null_space = "joint"` | profiled NSP |
| One margin \(k\) | Full sweep / multi-sweep fitter |
| Global \(P_k^{\mathrm{full}}\) | Own-margin surrogate |

### API (internal)

- C++: `tt_als_core_update_global_cpp(...)`
- R: `tt_als_core_update_global(..., backend = "R"|"Rcpp")`

### Gate quantities (R vs C++)

\[
S=Z_k^\top W Z_k,\quad
b=Z_k^\top W y_c,\quad
P_k^{\mathrm{full}},\quad
g_k^{\mathrm{new}},\quad
\mathrm{RSS},\quad
J_\lambda,\quad
L_{\mathrm{pen}}=Q
\]

Relative objective gate on well-conditioned problems:

\[
\frac{|L_{\mathrm{Cpp}}-L_R|}{1+|L_R|}<10^{-8}
\]

(also matrix / \(g\) absolute tolerances in `tests/testthat/test-als-core-update-global.R`).

### Cases

- \(d \in \{2,3,5\}\), rank \(\in \{1,2,3\}\)
- isotropic and anisotropic \(\lambda\)
- edge and interior cores
- extreme \(\lambda\) (near-singular stress)
- observation weights

### Next (not P1)

- **P2** one full fixed-λ sweep (same core order as R) — see below
- **P3** multi-sweep fitter + warm start / convergence
- **P4** outer cGCV calling C++ fixed-λ (optional)

## P2 — One fixed-λ sweep — **PASS** (2026-08-17)

Gate: `tests/testthat/test-als-sweep-global.R` — **135 PASS**.

### Scope

| In | Out |
|----|-----|
| One Gauss–Seidel sweep over all margins | Multi-sweep / tol / history |
| Fixed numeric \(\boldsymbol\lambda\) | cGCV / λ updates |
| LTR / RTL / arbitrary `margin_order` | Intercept refresh mid-sweep |
| Global \(P_k^{\mathrm{full}}\) | Own-margin |

### API (internal)

- C++: `tt_als_sweep_global_cpp(...)`
- R: `tt_als_sweep_global(..., backend = "R"|"Rcpp")`

Interfaces are rebuilt each core (≡ R `design_interface_cache = FALSE`). Gate compares \(\widehat f\), \(\eta\), RSS, \(J_\lambda\), \(L_{\mathrm{pen}}\), and cores (same init ⇒ unique conditional solves).

### Next

- **P3** multi-sweep fitter (`max_sweeps`, `tol`, warm start, no input mutation, convergence flags) — see below
- **P4** outer cGCV calling C++ fixed-λ (optional)

## P3 — Multi-sweep fixed-λ fitter — **PASS** (2026-08-17)

Gate: `tests/testthat/test-als-fit-fixed-global.R` — **84 PASS**.

### Scope

| In | Out |
|----|-----|
| `max_sweeps`, RSS `tol` (stop after sweep `> 2`) | cGCV / λ search |
| Intercept refresh each sweep | `linear=` / `smooth=` |
| Warm start via `init` cores | Profiled null-space |
| History / converged / reason | Public `ttps(backend="Rcpp")` switch |

### API (internal)

- C++: `tt_als_fit_fixed_global_cpp(...)`
- R: `tt_als_fit_fixed_global(..., backend = "R"|"Rcpp")`

This is the piece global-GCV should call for each numeric \(\boldsymbol\lambda\).

### Next

- **P4** (optional): outer cGCV / Sobol / GDF MC calling this fitter
- Wire `backend = "Rcpp"` for ALS only after production parity + speedup evidence
- Microbenchmark Gram vs sweep overhead before claiming speedups

## P4A — Wire global-GCV lab to P3 — **PASS** (2026-08-17)

Gate: `tests/testthat/test-global-gcv-p4a-backend.R`.

### API (lab / internal)

```r
fit_backend = c("R", "Rcpp_fixed")
```

Threaded through:

- `.tt_lab_fit_fixed()` — **single dispatcher** for all fixed-λ ALS
- `tt_global_gcv()`, `tt_global_gdf_mc()`, `.tt_lab_refit_from_base()`
- `.tt_lab_make_theta_evaluator()`, `.tt_lab_reeval_multistart()`
- `tt_global_lambda_optimize()`, `tt_global_lambda_optimize_v1()`

`R` = public [ttps()] ALS (unchanged algorithm).  
`Rcpp_fixed` = P3 `tt_als_fit_fixed_global(..., backend = "Rcpp")` wrapped as `"ttpspline"`.

**Not** in P4A: warm starts between λ, adaptive sweep budgets, fused_blocked, interface cache.

### Next

- End-to-end velocity script comparing `fit_backend` on Sobol/GDF workloads
- **P4B** adaptive fidelity + controlled GDF warm starts
- Then fused_blocked + interface cache inside P3

## Benchmark note

Script: `run_fixed_lambda_velocity.R` · report: `VELOCITY_FIXED_LAMBDA.md` (2026-08-17).

Measured end-to-end fixed-λ fit speedup **~1.1–1.9×** (Darwin arm64) vs R
reference; single-core micro ~1× (kernels already shared). Profile before
claiming larger gains from further C++ work.
