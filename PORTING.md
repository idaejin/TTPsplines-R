# Porting map (lab → package)

**Lab:** `../TTPsplines/` (unchanged)  
**Package folder:** `ttpsplines-pkg/` · **R package:** `TTPsplines`

## Done in v0 prototype

- Public API `ttpspline()` / controls / rank / complexity / S3  
- Gaussian ALS + cGCV (R + Rcpp)  
- Poisson / Bernoulli PIRLS  
- Modular `update_lambda()` (`fixed`, `cGCV`; stubs for `cFS`/`cREML`)  
- GLAM Gaussian grid helper  
- testthat smoke (15 tests)  
- vignette + benchmark stubs  

See `PACKAGE_STATUS.md`.

## Still copy/adapt as needed

| Lab | Note |
|---|---|
| `experiment_*.R` | → fill `inst/benchmarks/` |
| `R/glam_array.R` GLM grid | optional extend `glam_fit_*` |
| docs | keep in lab; cite from vignettes |

## Before GitHub

1. Replace `OWNER` in DESCRIPTION/README  
2. `devtools::document()`  
3. `devtools::check()`  
4. `usethis::use_git(); usethis::use_github()` on **this folder only**
