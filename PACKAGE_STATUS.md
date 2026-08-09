# TTPsplines package — prototype status (Paper 1)

**Folder:** `01_PROJECTS/ttpsplines-pkg/` (avoids macOS case-clash with lab `TTPsplines/`)  
**Package name:** `TTPsplines`  
**Version:** `0.0.0.9000`  
**Lab (unchanged source of truth for sims):** `01_PROJECTS/TTPsplines/`

## 1. Files created (high level)

```text
DESCRIPTION, NAMESPACE, LICENSE, README.md, PORTING.md
R/  TTPsplines-package, linalg, basis, penalties, rank, control,
    families, tt_geometry, lambda, als, pirls, ttpspline,
    complexity, methods, rank_profile, glam, RcppExports
src/ tt_pspline_nd.cpp (+ Makevars, RcppExports.cpp)
tests/testthat/test-ttpspline.R
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
- `tt_control()`, `tt_rank()`, `tt_complexity()`, `tt_rank_profile()`
- S3: `print`, `summary`, `predict`, `fitted`, `residuals`, `coef`, `deviance`, `plot`
- `glam_fit_gaussian()` for grid compression benchmarks
- Modular `update_lambda()` with `"fixed"` / `"cGCV"`; `"cFS"` / `"cREML"` reserved

## 4–8. Status by feature

| Feature | Status |
|---|---|
| Gaussian ALS fixed λ | **working** (R + Rcpp) |
| Gaussian cGCV | **working** (R + Rcpp) |
| Poisson PIRLS | **working** (R; Rcpp path available) |
| Bernoulli PIRLS | **working** (R; soft damping) |
| cGCV for GLM | **working** via modular λ / Rcpp |
| Rcpp backend | **working** (`backend = "auto"/"Rcpp"`) |
| Sparse Matrix path | **hook only** (`sparse` in control; dense v0) |
| Benchmarks | **runnable** — `inst/benchmarks/` (CSV/PNG in `results/`) |
| Tests | **15 PASS** (`testthat`) |
| Vignettes | **scaffolded** (`eval: false`) |

## 9. Remaining blockers for Paper 1 manuscript software

1. Replace `OWNER` GitHub URLs; `use_mit_license()` / push repo  
2. Joint EDF / AIC (optional; currently `edf = NA`)  
3. Fill `inst/benchmarks/` from lab smoke studies  
4. Optional: sparse B-spline storage; matrix-free \(S_k\)  
5. Stronger Bernoulli line-search; document λ bounds for separation  
6. Parity suite: lab scripts vs package on shared seeds  
7. `devtools::document()` for Rd pages  
8. Split `src/tt_pspline_nd.cpp` into `tt_core` / `tt_weighted` (cosmetic)

## 10. Paper 2 hooks

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
