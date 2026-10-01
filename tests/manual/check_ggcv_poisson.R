# Self-check: Poisson lambda = "gGCV" (working + tiny algorithmic).
# From ttpsplines-pkg root: Rscript tests/manual/check_ggcv_poisson.R
pkg_root <- Sys.getenv(
  "TTP_PKG",
  unset = "/Users/daejin/Dropbox/IE/research/01_PROJECTS/ttpsplines-pkg"
)
suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(pkg_root, quiet = TRUE)
  } else {
    library(TTPsplines)
  }
})

set.seed(2)
n <- 80L
X <- cbind(runif(n), runif(n))
eta <- 0.8 + sin(2 * pi * X[, 1]) + 0.4 * cos(2 * pi * X[, 2])
y <- rpois(n, exp(eta))

ctrl_w <- tt_control(
  max_sweeps = 6L,
  pirls_maxit = 12L,
  ggcv_glm_mode = "working",
  ggcv_n_global = 6L,
  ggcv_n_refine = 1L,
  ggcv_M_search = 3L,
  ggcv_M_final = 4L,
  compute_edf = FALSE,
  warn_lambda_boundary = FALSE,
  seed = 2L
)
fit_w <- ttps(y, X, family = poisson(), rank = 2L, k = 6L,
              lambda = "gGCV", control = ctrl_w)
stopifnot(inherits(fit_w, "ttpspline"))
stopifnot(identical(fit_w$lambda_method, "gGCV"))
stopifnot(identical(fit_w$ggcv$mode, "working"))
stopifnot(length(fit_w$lambda) == 2L, all(fit_w$lambda > 0))
cat("check_ggcv_poisson working: OK  lambda=",
    paste(sprintf("%.3g", fit_w$lambda), collapse = ","), "\n")

ctrl_a <- tt_control(
  max_sweeps = 4L,
  pirls_maxit = 8L,
  ggcv_glm_mode = "algorithmic",
  ggcv_poisson_anisotropic = FALSE,
  ggcv_n_global = 5L,
  ggcv_n_refine = 1L,
  ggcv_M_search = 2L,
  ggcv_M_final = 3L,
  ggcv_include_cgcv_anchor = FALSE,
  compute_edf = FALSE,
  warn_lambda_boundary = FALSE,
  seed = 2L
)
fit_a <- ttps(y, X, family = poisson(), rank = 2L, k = 6L,
              lambda = "gGCV", control = ctrl_a)
stopifnot(identical(fit_a$ggcv$mode, "algorithmic"))
stopifnot(is.finite(fit_a$ggcv$gcv), is.finite(fit_a$ggcv$gdf))
cat("check_ggcv_poisson algorithmic: OK  lambda=",
    paste(sprintf("%.3g", fit_a$lambda), collapse = ","),
    " gdf=", sprintf("%.2f", fit_a$ggcv$gdf), "\n")
