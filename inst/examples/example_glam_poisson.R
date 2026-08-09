# Currie / GLAM Poisson on an age × year grid
#   Rscript inst/examples/example_glam_poisson.R
# Run from the package root (or with TTPsplines installed).

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
  "\nGLAM-Poisson: npar=%d  PIRLS=%d  deviance=%.4g  RMSE(log mu)=%.4g\n",
  fit$npar, fit$n_pirls, fit$deviance, rmse_log
))
