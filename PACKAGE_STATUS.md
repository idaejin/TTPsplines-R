# TTPsplines package — prototype status

**Canonical status report:** [`docs/PACKAGE_IMPLEMENTATION.md`](docs/PACKAGE_IMPLEMENTATION.md)

**Folder:** `01_PROJECTS/ttpsplines-pkg/` (avoids macOS case-clash with lab `TTPsplines/`)  
**Package name:** `TTPsplines`  
**Version:** `0.0.0.9000`  
**Lab (unchanged source of truth for sims):** `01_PROJECTS/TTPsplines/`

## API axes

```text
optimizer ∈ {ALS, LBFGS, Adam*}
lambda    ∈ {scalar, length-d, "cGCV"}   # cFS/cREML not implemented yet
backend   ∈ {auto, R, Rcpp, keras*}
* Adam/keras = optional stub
```

Public entry: `ttpspline(..., optimizer=, backend=, init=)` plus `tt_initialize()`.

## 1. Files created (high level)

```text
DESCRIPTION, NAMESPACE, LICENSE, README.md, PORTING.md
R/  TTPsplines-package, linalg, basis, penalties, rank, control,
    initialization, families, tt_geometry, lambda, als, pirls,
    optimizer_lbfgs, optimizer_adam, ttpspline,
    complexity, methods, rank_profile, glam, RcppExports
src/ tt_pspline_nd.cpp (+ Makevars, RcppExports.cpp)
tests/testthat/test-ttpspline.R, test-api.R, test-complexity-layers.R
docs/PACKAGE_IMPLEMENTATION.md
vignettes/ introduction, generalized, rank-compression
inst/benchmarks/ *.R + python/benchmark_keras_tt.py (stubs)
```

## 2. Existing lab code reused

| Lab | Package use |
|---|---|
| `src/tt_pspline_nd.cpp` | copied → package Rcpp backend |
| TT ALS / interfaces / designs | reimplemented cleanly in `R/tt_geometry.R`, `R/als.R` |
| GLM PIRLS + cGCV ideas | `R/pirls.R`, `R/lambda.R` |
| GLAM RH | `R/glam.R` (`glam_fit_gaussian`) |
| Bases / difference penalties | `R/basis.R`, `R/penalties.R` |

Lab scripts/docs/outputs **not moved or deleted**.

## 3. API implemented

- `ttpspline(y, X, family, rank, k, lambda, control, …)`
- `tt_control()`, `tt_rank()`, `tt_initialize()`, `tt_complexity()`, `tt_rank_profile()`
- S3: `print`, `summary`, `predict`, `fitted`, `residuals`, `coef`, `deviance`, `plot`
- `glam_fit_gaussian()` for grid compression benchmarks
- Modular `update_lambda()` with `"fixed"` / `"cGCV"`; `"cFS"` / `"cREML"` reserved (not implemented)

## 4–8. Status by feature

| Feature | Status |
|---|---|
| Gaussian ALS fixed λ | **working** (R + Rcpp) |
| Gaussian cGCV | **working** (R + Rcpp) |
| Poisson PIRLS | **working** (R; Rcpp path available) |
| Bernoulli PIRLS | **working** (R; soft damping) |
| cGCV for GLM | **working** via modular λ / Rcpp |
| Joint EDF | **working** (linearized; size-guarded) |
| Complexity layers | **working** (`tt_complexity`) |
| Rcpp backend | **working** (`backend = "auto"/"Rcpp"`) |
| Sparse Matrix path | **hook only** (`sparse` in control; dense v0) |
| Benchmarks | **runnable** — `inst/benchmarks/` |
| Tests | **testthat** |
| Vignettes | **scaffolded** (`eval: false`) |

## 9. Remaining blockers

1. Stronger Bernoulli line-search; document λ bounds for separation  
2. Optional: sparse B-spline storage; matrix-free \(S_k\)  
3. Parity suite: lab scripts vs package on shared seeds  
4. `devtools::document()` for Rd pages  
5. Full Adam/Keras backend (currently stub)  
6. Split `src/tt_pspline_nd.cpp` into modules (cosmetic)

## 10. Reserved λ hooks (not implemented)

```r
lambda = "cFS"   # stop with clear message
lambda = "cREML"
update_lambda(method = ...)  # extend switch
```

## Minimal working examples

```r
devtools::load_all("01_PROJECTS/ttpsplines-pkg")

# Gaussian
fit <- ttpspline(y, X, family = gaussian(), rank = 3, k = 8, lambda = "cGCV")

# Poisson
fit <- ttpspline(y, X, family = poisson(),  rank = 3, k = 8, lambda = 1)

# Bernoulli
fit <- ttpspline(y, X, family = binomial(), rank = 3, k = 8, lambda = 1)
```
