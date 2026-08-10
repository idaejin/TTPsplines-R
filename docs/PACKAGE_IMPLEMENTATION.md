# Package implementation status

**Package:** `TTPsplines` (`01_PROJECTS/ttpsplines-pkg/`, GitHub `idaejin/TTPsplines-R`)  
**Research lab (validation reference):** `01_PROJECTS/TTPsplines/`  
**Last updated:** 2026-08-09

## Complexity layers (terminology)

`tt_complexity()` reports **distinct** objects — never a single “complexity”:

1. **Full tensor coefficients** \(N_{\mathrm{full}}=\prod p_k\)
2. **TT stored parameters** \(N_{\mathrm{TT}}=\sum r_{k-1}p_k r_k\)
3. **Compression ratio** \(CR=N_{\mathrm{full}}/N_{\mathrm{TT}}\)
4. **TT intrinsic dimension** \(\dim(\mathcal{M}_r)=N_{\mathrm{TT}}-\sum_{k=1}^{d-1}r_k^2\) (gauge)
5. **Joint EDF** (statistical; after \(\lambda\); **≠** \(N_{\mathrm{TT}}\))

Reading: \(r\) = structural capacity; \(\lambda\) = smoothness; EDF = effective fit complexity.

## Architecture

Public API separates three decisions:

| Axis | Argument | Options |
|------|----------|---------|
| Optimizer | `optimizer` | `ALS` (default), `LBFGS`, `Adam` (stub) |
| Smoothing | `lambda` | scalar / length-`d` fixed, `"cGCV"` |
| Backend | `backend` / `tt_control(backend=)` | `auto`, `R`, `Rcpp`, `keras` |

Internal λ dispatch:

```
parse_lambda_spec(lambda, d, control)
  → list(method, values, automatic)
update_lambda(method, workspace)
  → fixed | cGCV   (cFS / cREML not implemented yet)
```

Same statistical model for all optimizers: non-additive TT P-spline surface on scattered \(X \in \mathbb{R}^{n\times d}\).

## Files touched / added

| Path | Role |
|------|------|
| `R/lambda.R` | Parser + cGCV Brent + optional spectral cache |
| `R/control.R` | Extended `tt_control()` + `resolve_backend()` |
| `R/initialization.R` | `tt_initialize()` |
| `R/als.R` | ALS Gaussian; `init`; `values`; spectral flag |
| `R/pirls.R` | PIRLS GLM; `init`; weighted cGCV |
| `R/optimizer_lbfgs.R` | Global L-BFGS + outer cGCV |
| `R/optimizer_adam.R` | `tt_has_keras()` / `tt_keras_status()` + stub |
| `R/ttpspline.R` | `optimizer` / `backend` / `init` dispatch |
| `R/methods.R` | Summary shows optimizer / outer / λ evals |
| `DESCRIPTION` / `NAMESPACE` | Collate + exports |
| `tests/testthat/test-api.R` | λ / init / LBFGS / Adam stub tests |

## Public API

```r
ttps(y, X, family, rank, k, degree, penalty_order,
          lambda, optimizer, backend, init, control, knots)
tt_control(...); tt_rank(); tt_initialize(); tt_complexity(); tt_rank_profile()
tt_has_keras(); tt_keras_status()
# S3: print, summary, predict, fitted, residuals, coef, deviance, plot
```

## Family status

| Family | Status | Notes |
|--------|--------|-------|
| Gaussian | **Working** | ALS (R/Rcpp), LBFGS, fixed + cGCV |
| Poisson | **Working (ALS/PIRLS)** | Fixed + cGCV; LBFGS NLL path available |
| Bernoulli | **Working with caveats** | PIRLS damping; high-rank `|η|` can blow up |

## Optimizer status

| Optimizer | Status |
|-----------|--------|
| ALS | **Primary / default** — Gaussian + PIRLS GLM |
| LBFGS | **Implemented (R)** — joint cores; outer cGCV alternation |
| Adam/Keras | **Stub only** — clear error; no TF hard dependency |

## cGCV status

- Conditional GCV per core with Brent/`optimize()` on \(\eta=\log\lambda\)
- Workspace caches `S`, `b`, `P`, weighted `Xw`/`yw`
- Optional spectral cache (`use_spectral_gcv = TRUE`)
- ALS: update λ on each core visit; LBFGS/Adam: outer freeze-cores then conditional update
- **Not** global GCV; **not** cFS/cREML

## Rcpp / sparse

| Item | Status |
|------|--------|
| Rcpp ALS / cGCV / GLM PIRLS | Present via `src/tt_pspline_nd.cpp` when compiled |
| Sparse `Matrix` hybrid | Hook (`sparse=` in control); dense default in v0 |
| Matrix-free \(X_k\) | Experimental / not default |

## Tests

- Gaussian fixed/cGCV, poisson, bernoulli, predict=fitted, reserved λ methods
- λ expansion/validation, `tt_initialize` reproducibility, ALS vs LBFGS surface proximity, Adam stub message, anisotropic λ
- Complexity layers (storage / intrinsic / EDF)

## Benchmarks

Lab + `inst/benchmarks/` remain the reference. Optimizer comparison benchmarks still to be fleshed out; do not block ALS/cGCV correctness.

## Speedups / parity

- Prefer end-to-end timings in the research lab audit notes
- Numerical parity vs lab scripts: compare **η / μ / objective**, never raw cores (gauge)

## Blockers / deferred

1. **Adam/Keras full implementation** — deferred until ALS + cGCV + Rcpp + GLM are solid  
2. Joint EDF — **implemented** (linearized stacked Jacobian; skip via `edf_max_npar`)  
3. Sparse / matrix-free backends — experimental  
4. Automated lab↔package parity suite  
5. Bernoulli stability at high rank  

## Not implemented yet (reserved API hooks)

- `lambda = "cFS"`, `"cREML"`
- TT-cSOP, Schall, DMRG, Riemannian, MALS  

---

## Minimal examples (same model; only optimizer / λ / backend change)

```r
# Fixed isotropic
fit <- ttps(y, X, rank = 3, lambda = 1)

# Fixed anisotropic
fit <- ttps(y, X, rank = 3, lambda = c(0.1, 1, 10, 1, 0.5))

# Automatic cGCV
fit <- ttps(y, X, rank = 3, lambda = "cGCV")

# L-BFGS + cGCV
fit <- ttps(y, X, rank = 3, optimizer = "LBFGS", lambda = "cGCV")

# Adam — errors with install/status guidance until implemented
# fit <- ttps(y, X, rank = 3, optimizer = "Adam", lambda = "cGCV")

# Fair optimizer comparison
init <- tt_initialize(X, rank = 3, seed = 123)
fit_als   <- ttps(y, X, rank = 3, lambda = 1, optimizer = "ALS",   init = init)
fit_lbfgs <- ttps(y, X, rank = 3, lambda = 1, optimizer = "LBFGS", init = init)
```

**Priority:** correctness > reproducibility > stability > efficiency > API polish.
