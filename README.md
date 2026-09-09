# TTPsplines

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**Experimental** R package for **Tensor-Train P-splines**: non-additive multidimensional **statistical** smooth / GLM regression on **scattered** continuous covariates.

The TT factorization compresses the **coefficient tensor** \(\Theta\) of a tensor-product P-spline (no observation grid required) and uses classical **directional discrete-difference (P-spline) penalties**. That coefficient geometry already appears in tensor-network B-splines for system identification; this package does **not** claim priority for “TT of B-spline weights + difference penalties.” Its focus is **smoother practice**: GLM families, TT-aware conditional GCV (default), optional k-fold `lambda = "CV"`, experimental joint TT-gGCV, predictive rank selection helpers, array-mode / GLAM baselines when feasible, and an open fitting API (`ttps()`).

Within that geometry the package keeps \(r \neq \lambda \neq \mathrm{EDF}\) conceptually distinct (rank is structural capacity; \(\lambda\) is roughness; EDF is a post-penalty diagnostic).

## Install

**Package only** (fast; no vignettes):

```r
# install.packages("pak")
pak::pak("idaejin/TTPsplines-R")
```

**Package + vignettes** (needed for `vignette(...)` / `browseVignettes()`):

```r
# install.packages("remotes")
remotes::install_github(
  "idaejin/TTPsplines-R",
  force = TRUE,              # required if remotes skips an unchanged SHA
  build_vignettes = TRUE,
  dependencies = TRUE        # knitr, rmarkdown, ...
)
```

From a local clone:

```r
devtools::install(build_vignettes = TRUE)
# or without vignettes: devtools::load_all()
```

`pak` and `devtools::load_all()` do **not** register vignettes in the library.

### Vignettes

After a vignette-enabled install:

```r
browseVignettes("TTPsplines")
vignette("getting-started", package = "TTPsplines")
```

| Vignette | Topic |
|----------|--------|
| `getting-started` | Scattered TT fit, families, datasets |
| `array-mode` | Product grids with `array = TRUE` |
| `cgcv` | λ: fixed / cGCV / CV / gGCV |
| `rank-selection` | TT rank via CV + 1-SE |
| `margin-activity-path` | Margin screening |
| `generalized` | Poisson / Bernoulli; `linear=` / `smooth=` |
| `aic-bic` | In-sample AIC / BIC from EDF |
| `glam-vs-tt` | Full tensor vs TT |
| `glam-poisson` | GLAM Poisson + exposure |
| `scalability` | Storage / timing notes |
| `uncertainty` | SE / bands (Level-1) |

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
yg <- rnorm(nrow(X))  # any Gaussian response
fit_g_lb <- ttps(yg, X, family = gaussian(), rank = 2, k = 8, lambda = 1,
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

ALS / PIRLS sweeps always run in **R** under the classical global penalty
\(P_k^{\mathrm{full}}\). Compiled **Rcpp kernels** (Gram/RHS, \(P_k^{\mathrm{full}}\)
assembly, etc.) are used from that R loop when available. Requesting
`backend = "Rcpp"` for ALS/PIRLS does **not** switch to a full C++ fitter; it
emits a warning and still uses the R sweep (with Rcpp helpers). Prefer
`monitor = TRUE` (or `tt_control(trace = TRUE)`) to print per-sweep progress.

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

## Complete grids (`array = TRUE`)

For a full rectangular response array (Gaussian ALS), pass the array and
marginal axes so Gram/RHS use Kronecker structure without materialising the
scattered \(n\times q_k\) design:

```r
set.seed(1)
n1 <- 12; n2 <- 10; n3 <- 8
x1 <- seq(0, 1, length.out = n1)
x2 <- seq(0, 1, length.out = n2)
x3 <- seq(0, 1, length.out = n3)
g <- expand.grid(x1 = x1, x2 = x2, x3 = x3)
mu <- sin(2 * pi * g$x1) * cos(2 * pi * g$x2) + 0.5 * g$x3
# dim1-fastest layout (same as array(y, dim = c(n1, n2, n3)))
Y <- array(mu + rnorm(n1 * n2 * n3, sd = 0.2), dim = c(n1, n2, n3))

fit_a <- ttps(
  Y, axes = list(x1, x2, x3),
  rank = 2, k = 8, lambda = 1,
  array = TRUE,
  control = tt_control(max_sweeps = 8, backend = "auto")
)
summary(fit_a)
```

Restrictions in this version: no `linear=` / `smooth=` /
`null_space = "profiled"`; no `lambda = "CV"` / `"gGCV"`. Numerically matches
scattered `ttps()` on the same grid. Gaussian unweighted grids use Kronecker
Gram; Poisson / weighted grids use the weighted array path. For
exposure-weighted Poisson with **dense** \(\Theta\) prefer
`glam_fit_poisson()`. Full walkthrough:
`vignette("array-mode", package = "TTPsplines")`.

## Choosing the TT rank (CV + 1-SE)

`ttps(..., rank = r)` always uses that exact rank (no auto-selection).
For predictive choice of \(r\):

```r
sel <- tt_rank_select(y, X, ranks = 1:5, lambda = 1, folds = 5, rule = "1se")
# init-sensitive problems (e.g. Ishigami at low r):
# sel <- tt_rank_select(..., n_starts = 5)
sel
plot(sel)
fit <- tt_rank_refit(sel)   # full-data refit at selected_rank
```

Vignette: `vignette("rank-selection", package = "TTPsplines")`.
Warm-start from a neighbouring rank: `tt_truncate_rank(fit$cores, rank = 2)`.

## Margin Activity Path (which covariates)

When many margins may be null, screen columns before a full-\(d\) cGCV fit:

```r
path <- tt_margin_activity_path(
  y, X, rank = 2, k = 5, select = "1se", folds = 5, seed = 1
)
path$selected_names
plot(path)
fit <- path$fit   # TT + cGCV on selected margins
```

Vignette: `vignette("margin-activity-path", package = "TTPsplines")`.
Example (reproducible): `Rscript inst/examples/reprex_margin_activity_path.R`.
Extended demo: `Rscript inst/examples/example_margin_activity_path.R`.
This chooses the **margin set**; it does not replace `tt_rank_select()` (\(r\))
or `lambda = "cGCV"` (smoothness). Prefer the full name *Margin Activity Path*
(avoid the acronym "MAP").

Complementary leave-one-out / permutation drop diagnostic:

```r
tst <- tt_margin_drop_test(y, X, rank = 2, k = 5, lambda = 1, method = "nested")
tst$drop_candidate_names
```

## Choosing λ (fixed / cGCV / CV / gGCV)

| Spec | Role | Notes |
|------|------|-------|
| numeric / length-`d` | Fixed isotropic or anisotropic \(\lambda\) | Always available |
| `"cGCV"` (default) | Conditional / product GCV inside ALS / PIRLS | Default dynamics: `cgcv_update = "outer_simultaneous"` |
| `"CV"` | K-fold CV of each \(\lambda_k\) | ALS / PIRLS only; default `cv_sweeps = 1` (tune then freeze) |
| `"gGCV"` | Joint TT-gGCV (Monte Carlo GDF) | Experimental; Gaussian scattered only; expensive — prefer `tt_ggcv()` |


```r
# Default product selector
fit <- ttps(y, X, rank = 2, k = 8, lambda = "cGCV")
fit$lambda
fit$lambda_boundary   # check for bound hits

# K-fold CV (ALS / PIRLS only). Default: tune on sweep 1, then freeze.
fit_cv <- ttps(
  y, X, rank = 2, k = 8, lambda = "CV",
  control = tt_control(cv_folds = 5, cv_rule = "min", cv_sweeps = 1)
)
fit_cv$cv          # folds, grid, rule, trace

# Joint TT-gGCV (opt-in oracle; prefer tt_ggcv() for search diagnostics)
# fit_g <- ttps(y, X, rank = 2, k = 6, lambda = "gGCV",
#               control = tt_control(ggcv_n_global = 16, ggcv_M_search = 4))
# opt <- tt_ggcv(y, X, rank = 2, k = 6, n_global = 16, M_search = 4)
```

Useful knobs: `tt_control(cv_folds, cv_ngrid, cv_grid, cv_sweeps, cv_rule)` and
`tt_control(ggcv_n_global, ggcv_n_refine, ggcv_M_search, ggcv_M_final, ggcv_include_cgcv_anchor)`.

If cGCV λ sits on a search bound, diagnose with multi-start (stable vs unstable
hits) via `ttps_multistart()` — see `vignette("cgcv")`.

Vignette: `vignette("cgcv", package = "TTPsplines")`.

## EDF (joint and per margin)

```r
fit <- ttps(y, X, rank = 2, k = 8, lambda = "cGCV")  # compute_edf = TRUE by default
fit$edf          # joint linearized EDF
fit$edf_margin   # block tr(H_kk); sums to fit$edf
fit$edf_margin_cond  # ALS/cGCV conditional traces (diagnostic; not additive)
tt_edf(fit)      # tidy extract / recompute
```

`sum(fit$edf_margin)` equals `fit$edf` (parameter-block partition). Use `edf_margin_cond` only to diagnose cGCV core flexibility.

## AIC / BIC (linearized EDF)

In-sample information criteria use joint linearized `fit$edf` (not `npar_tt`)
via exported `tt_ic()`:

```r
tt_ic(fit, "AIC")   # Gaussian: n * log(RSS/n) + 2 * (edf + 1)
tt_ic(fit, "BIC")   # Poisson / Bernoulli: deviance + pen * (edf + 1)
```

Vignette: `vignette("aic-bic", package = "TTPsplines")`.
There is no `AIC()` / `BIC()` S3 method yet (`tt_ic()` is the API).

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
Vignettes: `vignette("array-mode")`, `vignette("glam-vs-tt")`,
`vignette("glam-poisson")`, `vignette("scalability")`.
Scripts: `inst/examples/example_glam_poisson.R`,
`inst/examples/example_glam_gaussian_vs_tt.R`.

| In v0 | Not yet |
|---|---|
| TT-ALS / PIRLS / global L-BFGS | rank-/λ-adaptive `auto` |
| Family-aware `auto` optimizer | rank-/λ-adaptive `auto` |
| Gaussian / Poisson / Bernoulli | SA-CAB, SOP, DMRG |
| `lambda` fixed / `"cGCV"` / `"CV"`; experimental `"gGCV"` / `tt_ggcv()` | automatic rank inside `ttps()` |
| `tt_rank_select()` + `tt_rank_refit()` | LRT / bootstrap rank tests |
| `tt_margin_activity_path()` (margin screening) | group-lasso on TT cores |
| Experimental: `GD`, `Damped-Newton-ALS`, `LBFGS-ALS` | mixed effects / TMB |
| `array = TRUE` (Gaussian / Poisson grids) + GLAM Poisson | higher-d GLAM / REML; DLNM vignette |

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
