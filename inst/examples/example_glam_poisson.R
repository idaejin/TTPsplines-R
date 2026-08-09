# Currie–Durbán–Eilers GLAM Poisson on an age × year grid
#   Rscript inst/examples/example_glam_poisson.R
# Run from the package root (or with TTPsplines installed).
#
# Reference: Currie, Durbán & Eilers (2006), JRSS-B — generalized linear
# array models (GLAM) with rotated-H algebra; Poisson PIRLS on a grid.

root <- if (file.exists("DESCRIPTION") &&
             identical(unname(read.dcf("DESCRIPTION")[, "Package"]), "TTPsplines")) {
  normalizePath(".")
} else if (file.exists("../../DESCRIPTION")) {
  normalizePath("../..")
} else {
  NULL
}

if (!is.null(root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, quiet = TRUE)
} else {
  library(TTPsplines)
}

cat("Currie–Durbán–Eilers GLAM Poisson (age × year grid)\n")
data(glam_poisson)
cat("Grid:", paste(dim(glam_poisson$Y), collapse = " × "),
    " (age × year)\n")

bb <- glam_grid_bases(
  list(age = glam_poisson$age, year = glam_poisson$year),
  k = 10
)
fit <- glam_fit_poisson(
  glam_poisson$Y, bb$B,
  lambda = c(10, 1),
  offset = log(glam_poisson$exposure),
  trace = TRUE
)

rmse_log <- sqrt(mean((log(as.numeric(fit$mu)) -
                         log(as.numeric(glam_poisson$mu)))^2))
cat(sprintf(
  "\nCurrie–Durbán–Eilers GLAM-Poisson: npar=%d  PIRLS=%d  deviance=%.4g  RMSE(log mu)=%.4g\n",
  fit$npar, fit$n_pirls, fit$deviance, rmse_log
))

# Same grid as scattered data → TTPsplines (shared exposure offset)
X <- as.matrix(expand.grid(age = glam_poisson$age, year = glam_poisson$year))
y <- as.numeric(glam_poisson$Y)
off <- log(as.numeric(glam_poisson$exposure))
fit_tt <- ttpspline(
  y, X, family = poisson(), rank = 3, k = 10, lambda = c(10, 1),
  offset = off,
  control = tt_control(
    pirls_maxit = 20, max_sweeps = 8, backend = "R", compute_edf = FALSE
  )
)
rmse_tt <- sqrt(mean((log(fitted(fit_tt)) - log(as.numeric(glam_poisson$mu)))^2))
cat(sprintf(
  "TTPsplines rank-3 (+offset): npar=%d (dense %d, CR %.1fx)  RMSE(log mu)=%.4g  vs GLAM=%.4g\n",
  fit_tt$npar_tt, fit_tt$npar_dense, fit_tt$compression_ratio, rmse_tt,
  sqrt(mean((log(fitted(fit_tt)) - log(as.numeric(fit$mu)))^2))
))
