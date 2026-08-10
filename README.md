# TTPsplines

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**Experimental** R package for **Tensor-Train P-splines**: non-additive multidimensional smooth / GLM regression on **scattered** continuous covariates.

The TT factorization compresses the **coefficient tensor** of a tensor-product P-spline — observations need **not** lie on a grid.

## Install

```r
# install.packages("pak")
pak::pak("idaejin/TTPsplines-R")
# or: remotes::install_github("idaejin/TTPsplines-R")
```

To register vignettes in the installed library (needed for `vignette(...)`):

```r
# from a clone of this repo:
devtools::install(build_vignettes = TRUE)
# or from GitHub:
remotes::install_github("idaejin/TTPsplines-R", build_vignettes = TRUE)
```

Then open with an explicit package (or `library(TTPsplines)` first):

```r
vignette("getting-started", package = "TTPsplines")
browseVignettes("TTPsplines")
```

`devtools::load_all()` does **not** install vignettes into the library, so
`vignette("getting-started")` will still warn until you `install(..., build_vignettes = TRUE)`.

Local (development folder `ttpsplines-pkg/`):

```r
devtools::load_all("path/to/ttpsplines-pkg")
# or: pak::local_install("path/to/ttpsplines-pkg")
```

## Quick start

### Gaussian (`auto` → ALS)

```r
library(TTPsplines)

set.seed(1)
n <- 800
X <- matrix(runif(n * 3), n, 3)
f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.5 * X[, 3]
y <- f + rnorm(n, sd = 0.25)

fit <- ttps(
  y, X,
  family = gaussian(),
  rank = 2,
  k = 10,
  lambda = "cGCV",
  control = tt_control(max_sweeps = 12, backend = "auto")
)

summary(fit)
tt_complexity(fit)
predict(fit, X[1:5, ], type = "response")
```

With `optimizer = "auto"` (default), `summary(fit)` reports the family rule explicitly, e.g.:

```text
Requested optimizer:    auto
Selected optimizer:     ALS
Reason:                 gaussian family default
```

### Poisson (`auto` → PIRLS-ALS)

```r
library(TTPsplines)

set.seed(2)
n <- 800
X <- matrix(runif(n * 3), n, 3)
f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.4 * X[, 3]
eta <- f - mean(f) + log(3)
y <- rpois(n, exp(eta))

fit_p <- ttps(
  y, X,
  family = poisson(),
  rank = 2,
  k = 8,
  lambda = 1,
  control = tt_control(pirls_maxit = 25, als_sweeps_per_pirls = 4,
                       backend = "auto")
)

summary(fit_p)
# Requested optimizer: auto | Selected: PIRLS-ALS | poisson family default
predict(fit_p, X[1:5, ], type = "response")
```

### Bernoulli (`auto` → LBFGS)

```r
library(TTPsplines)

set.seed(3)
n <- 800
X <- matrix(runif(n * 3), n, 3)
f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.4 * X[, 3]
eta <- 1.5 * (f - mean(f))
y <- rbinom(n, 1, plogis(eta))

fit_b <- ttps(
  y, X,
  family = binomial(),
  rank = 2,
  k = 8,
  lambda = 5,
  control = tt_control(lbfgs_maxit = 300, backend = "auto")
)

summary(fit_b)
# Requested optimizer: auto | Selected: LBFGS | binomial family default
predict(fit_b, X[1:5, ], type = "response")
```

## Family-aware `auto` (v1)

| Family | Selected optimizer |
|---|---|
| Gaussian | `ALS` |
| Poisson | `PIRLS-ALS` |
| binomial | `LBFGS` |

Always overridable:

```r
# Explicit overrides (research / benchmarking)
fit_b_als <- ttps(y, X, family = binomial(), rank = 2, k = 8, lambda = 5,
                       optimizer = "PIRLS-ALS")
fit_g_lb  <- ttps(yg, X, family = gaussian(), rank = 2, k = 8, lambda = 1,
                       optimizer = "LBFGS")
```

Fixed anisotropic λ: `lambda = c(1, 10, 0.5)`.

Fit fields: `optimizer_requested`, `optimizer_used`, `optimizer_reason`
(and `optimizer` ≡ `optimizer_used` for compatibility).

### Monitor progress

```r
# either form enables the same iteration logs
fit <- ttps(y, X, family = gaussian(), rank = 2, k = 8, lambda = 1,
                 monitor = TRUE)
fit <- ttps(y, X, family = gaussian(), rank = 2, k = 8, lambda = 1,
                 control = tt_control(monitor = TRUE, max_sweeps = 12))
# equivalent: tt_control(trace = TRUE)
```

With `monitor = TRUE` and `backend = "auto"`, the package uses the **R** path so ALS/PIRLS sweep lines are printed. Pass `backend = "Rcpp"` explicitly if you prefer the fast path without per-sweep logs.

### cGCV λ near search bounds

Default search interval is `lambda_bounds = c(1e-4, 1e4)`. If selected λ hug the edges, `summary(fit)` reports:

```text
Lambda:                 9922.5, 0.0431, 0.000100
Lambda boundary:        upper, interior, lower
Lambda search bounds:   [0.0001, 10000]
```

and (unless `warn_lambda_boundary = FALSE`) emits a soft warning. Fields: `fit$lambda_boundary`, `fit$lambda_bounds`, `fit$lambda_at_boundary`.

## Example datasets (Ishigami / Sobol-g / Friedman)

```r
library(TTPsplines)
data(ishigami)   # d=3
data(sobol_g)    # d=4
data(friedman)   # d=5

X <- as.matrix(ishigami[, c("x1", "x2", "x3")])
fit <- ttps(ishigami$y, X, rank = 2, k = 8, lambda = 1)
summary(fit)
```

On-the-fly redraws: `simulate_ishigami()`, `simulate_sobol_g()`,
`simulate_friedman()`. Truth helpers: `f_ishigami()`, `f_sobol_g()`,
`f_friedman()`.

Vignette: `vignette("getting-started", package = "TTPsplines")`.
Script: `Rscript inst/examples/example_test_functions.R`.

## Choosing the TT rank (CV + 1-SE)

`ttps(..., rank = r)` always uses that exact rank (no auto-selection).
For predictive choice of \(r\):

```r
sel <- tt_rank_select(y, X, ranks = 1:5, lambda = 1, folds = 5, rule = "1se")
sel
plot(sel)
fit <- tt_rank_refit(sel)   # full-data refit at selected_rank
```

Vignette: `vignette("rank-selection", package = "TTPsplines")`.

## Choosing λ (cGCV)

```r
fit <- ttps(y, X, rank = 2, k = 8, lambda = "cGCV")
fit$lambda
fit$lambda_boundary   # check for bound hits
```

Vignette: `vignette("cgcv", package = "TTPsplines")`.

## AIC / BIC (linearized EDF)

In-sample information criteria use joint linearized `fit$edf` (not `npar_tt`):

```r
# Gaussian: n * log(RSS/n) + 2 * (edf + 1)
# Poisson / Bernoulli: deviance + 2 * (edf + 1)
```

Vignette: `vignette("aic-bic", package = "TTPsplines")`.
There is no `AIC()` / `BIC()` method yet.

## GLAM Poisson (Currie–Durbán–Eilers)

On a regular age × year grid with exposures (Currie, Durbán & Eilers, 2006):

```r
data(glam_poisson)
bb <- glam_grid_bases(list(age = glam_poisson$age, year = glam_poisson$year), k = 10)
fit <- glam_fit_poisson(
  glam_poisson$Y, bb$B, lambda = c(10, 1),
  offset = log(glam_poisson$exposure)
)
```

Also: `glam_fit_gaussian()`, `simulate_glam_poisson()`,
`compare_glam_tt_gaussian()` / `compare_glam_tt_scale()` (Gaussian GLAM vs TT
on \(d=3,5,7\) grids, including \(n\times k\) scale at \(d=7\)).
Vignettes: `vignette("glam-vs-tt")`, `vignette("glam-poisson")`,
`vignette("scalability")`.
Scripts: `inst/examples/example_glam_poisson.R`,
`inst/examples/example_glam_gaussian_vs_tt.R`.

| In v0 | Not yet |
|---|---|
| TT-ALS / PIRLS / global L-BFGS | TT-cFS, cREML |
| Family-aware `auto` optimizer | rank-/λ-adaptive `auto` |
| Gaussian / Poisson / Bernoulli | SA-CAB, SOP, DMRG |
| `lambda` fixed / `"cGCV"` | automatic rank inside `ttps()` |
| `tt_rank_select()` + `tt_rank_refit()` | LRT / bootstrap rank tests |
| Experimental: `GD`, `Damped-Newton-ALS`, `LBFGS-ALS` | mixed effects / TMB |
| GLAM Gaussian + Poisson (fixed λ, d≤3) | higher-d GLAM / REML |

## License

MIT

## Benchmarks

Not run in `R CMD check`. From the package root:

```r
devtools::load_all()
source("inst/benchmarks/run_all.R")
```

Or one family:

```bash
Rscript inst/benchmarks/benchmark_gaussian.R
TTPSPLINES_BENCH_WHICH=poisson,bernoulli Rscript inst/benchmarks/run_all.R
```

Results: `inst/benchmarks/results/*.csv` (+ PNG). See `inst/benchmarks/README.md`.

Quick three-family smoke (or per-family scripts):

```bash
Rscript inst/examples/example_three_families.R
Rscript inst/examples/example_poisson.R
Rscript inst/examples/example_bernoulli.R
Rscript inst/examples/example_test_functions.R
```
