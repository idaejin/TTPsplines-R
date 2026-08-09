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

Local (development folder `ttpsplines-pkg/`):

```r
devtools::load_all("path/to/ttpsplines-pkg")
# or: pak::local_install("path/to/ttpsplines-pkg")
```

## Quick start

```r
library(TTPsplines)

set.seed(1)
n <- 800
X <- matrix(runif(n * 3), n, 3)
f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.5 * X[, 3]
y <- f + rnorm(n, sd = 0.25)

fit <- ttpspline(
  y, X,
  family = gaussian(),
  rank = 3,
  k = 8,
  lambda = "cGCV",
  control = tt_control(max_sweeps = 12, backend = "auto")
)

summary(fit)
tt_complexity(fit)
predict(fit, X[1:5, ], type = "response")
```

Poisson / Bernoulli:

```r
# y ~ Poisson(exp(f));  y ~ Bernoulli(plogis(f))
fit_p <- ttpspline(y, X, family = poisson(),  rank = 3, k = 8, lambda = 1)
fit_b <- ttpspline(y, X, family = binomial(), rank = 3, k = 8, lambda = 1)
```

Fixed anisotropic λ: `lambda = c(1, 10, 0.5)`.

## Current scope

| In v0 | Not yet |
|---|---|
| TT-ALS / PIRLS | TT-cFS, cREML |
| Gaussian / Poisson / Bernoulli | SA-CAB, SOP, DMRG |
| `lambda` fixed / `"cGCV"` | automatic rank |
| GLAM grid baseline helper | mixed effects / TMB |

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

Quick three-family smoke:

```bash
Rscript inst/examples/example_three_families.R
```
