# Self-check: public lambda = "gGCV" returns a fit (smoke budgets).
# From ttpsplines-pkg root: Rscript tests/manual/check_ggcv_api.R
pkg_root <- "/Users/daejin/Dropbox/IE/research/01_PROJECTS/ttpsplines-pkg"
suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(pkg_root, quiet = TRUE)
  } else {
    library(TTPsplines)
  }
})

set.seed(1)
n <- 120L
X <- cbind(runif(n), runif(n))
y <- sin(2 * pi * X[, 1]) + 0.3 * sin(pi * X[, 2]) + rnorm(n, sd = 0.25)
ctrl <- tt_control(
  max_sweeps = 6L,
  ggcv_n_global = 8L,
  ggcv_n_refine = 1L,
  ggcv_M_search = 3L,
  ggcv_M_final = 6L,
  compute_edf = FALSE,
  warn_lambda_boundary = FALSE
)

fit <- ttps(y, X, rank = 2L, k = 6L, lambda = "gGCV", control = ctrl)
stopifnot(inherits(fit, "ttpspline"))
stopifnot(identical(fit$lambda_method, "gGCV"))
stopifnot(length(fit$lambda) == 2L, all(fit$lambda > 0))
stopifnot(!is.null(fit$ggcv$gcv), is.finite(fit$ggcv$gcv))
cat("check_ggcv_api: OK\n")
