# TTPsplines package — prototype status

**Canonical implementation notes:** [`docs/PACKAGE_IMPLEMENTATION.md`](docs/PACKAGE_IMPLEMENTATION.md)

**Folder:** `01_PROJECTS/ttpsplines-pkg/` (avoids macOS case-clash with lab `TTPsplines/`)  
**Package name:** `TTPsplines`  
**Version:** `0.0.0.9001`  
**GitHub:** https://github.com/idaejin/TTPsplines-R  
**Lab:** `01_PROJECTS/TTPsplines/` — consumes this package for the engine.

**DECISION (2026-08-10):** package = sole home of engine unit tests.

This file is a **living snapshot** (last refresh 2026-09-09). Prefer NEWS +
vignettes for user-facing truth; prefer `docs/PACKAGE_IMPLEMENTATION.md` for
deep engine notes.

## API axes

```text
optimizer ∈ {auto*, ALS, PIRLS-ALS, Damped-Newton-ALS*, LBFGS-ALS*, GD*, LBFGS, hybrid*}
lambda    ∈ {scalar, length-d, "cGCV", "CV", "gGCV"*}
backend   ∈ {auto, R, Rcpp}
array     ∈ {FALSE (scattered), TRUE (product grid Y + axes)}
* auto: Gaussian→ALS, Poisson→PIRLS-ALS, binomial→LBFGS
* gGCV: experimental, Gaussian scattered only
```

Public entry: `ttps(...)` (+ `ttps_dlnm()` for exposure×lag; vignette deferred).

## Vignettes (all knit `ok` in `vignettes/_build_check.csv`)

| Vignette | Topic |
|----------|--------|
| `getting-started` | Scattered TT, datasets |
| `array-mode` | Product grids `array = TRUE` |
| `generalized` | GLM families; `linear=` / `smooth=` |
| `cgcv` | λ selectors |
| `rank-selection` | CV + 1-SE for \(r\) |
| `aic-bic` | `tt_ic()` AIC/BIC from EDF |
| `margin-activity-path` | Margin screening |
| `glam-vs-tt` / `glam-poisson` | Dense GLAM baselines |
| `scalability` | Storage / timings / backends |
| `uncertainty` | Conditional SE / bands |

Manual smoke: `tests/manual/smoke_test_apis.Rmd`.

## Status by feature (short)

| Feature | Status |
|---|---|
| Gaussian ALS / cGCV | working |
| Poisson PIRLS / Bernoulli LBFGS (`auto`) | working |
| `array = TRUE` (Gaussian Kronecker; Poisson weighted) | working; no linear/smooth/CV/gGCV |
| `tt_ic()` AIC/BIC | working (`AIC()` S3 not yet) |
| Rank helpers `tt_rank_select` / `tt_rank_refit` | working |
| `linear=` / `smooth=` (ALS/PIRLS) | working; vignette section in `generalized` |
| Dense GLAM Gaussian / Poisson | working |
| `ttps_dlnm` / `predict_dlnm` | working API; **no vignette yet** (deferred) |
| Rcpp | kernels only; ALS sweeps stay in R |
| Joint gGCV | experimental; Gaussian scattered |

## Minimal examples

```r
devtools::load_all("01_PROJECTS/ttpsplines-pkg")

fit <- ttps(y, X, family = gaussian(), rank = 3, k = 8, lambda = "cGCV")
fit <- ttps(Y, axes = list(...), array = TRUE, rank = 2, k = 8, lambda = 1)
tt_ic(fit, "AIC")
```
