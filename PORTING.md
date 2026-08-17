# Porting map (lab ↔ package)

**Lab:** `../TTPsplines/` (private `idaejin/TTPsplines`) — papers, lit, manuscript, Paper-2 prototypes.  
**Package folder:** `ttpsplines-pkg/` · **R package:** `TTPsplines` · **GitHub:** `idaejin/TTPsplines-R`

## Architecture (DECISION 2026-08-10)

- **Package** owns the fitting engine and unit tests.
- **Lab consumes** the installed package for Paper-1 engine demos (`library(TTPsplines)` / `R/00_load_package.R`).
- Lab retains Paper-2 prototypes (TT-cFS / cREML / SA-CAB / SOP) until they are ported here.
- Lab `R/tt_pspline_nd.R`, `tt_glm_pirls.R`, … → `R/_archive_engine/` (historical; Paper-2 scaffolding)
- Lab engine unit tests removed; package `tests/testthat/` is authoritative

## Done in v0 prototype

- Public API `ttps()` / controls / rank / complexity / S3  
- Gaussian ALS + cGCV (R sweeps; Rcpp Gram / \(P_k^{\mathrm{full}}\) helpers)  
- Poisson / Bernoulli PIRLS (+ family-aware `optimizer="auto"`)  
- Conditional experimental solvers (Damped-Newton-ALS, LBFGS-ALS, GD)  
- Example surfaces: Ishigami / Sobol-g / Friedman (`simulate_*`, `data(...)`)  
- Modular `update_lambda()` (`fixed`, `cGCV`; stubs for `cFS`/`cREML`)  
- GLAM Gaussian grid helper  
- testthat suite  
- vignette + benchmark stubs  

See `PACKAGE_STATUS.md`.

## Still copy/adapt as needed

| Lab | Note |
|---|---|
| `tt_sop.R` / `tt_cstar_nd.R` / `experiment_tt_cfs_*` | → port TT-cFS / cREML when Paper 2 freezes API |
| `sa_cab_lambda.R` / SOP/CSOP | stay lab-only until justified |
| `experiment_*.R` (Paper 1) | prefer package API; fill `inst/benchmarks/` as needed |
| `R/glam_array.R` GLM grid | **done** — `glam_fit_poisson` + `glam_grid_bases` |
| docs / lit / manuscript | remain in lab |

## Install (collaborators)

```r
pak::pak("idaejin/TTPsplines-R")
library(TTPsplines)
```
