# Classic test surfaces: Ishigami / Sobol-g / Friedman
#   Rscript inst/examples/example_test_functions.R

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required for this example")
}
root <- (function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), "../..")))
  }
  normalizePath(".")
})()
devtools::load_all(root, quiet = TRUE)

.fit_demo <- function(dat, rank = 2L, k = 6L, lambda = 1) {
  fam <- switch(dat$family,
                gaussian = gaussian(),
                poisson = poisson(),
                binomial = binomial())
  ttps(
    dat$y, dat$X, family = fam, rank = rank, k = k, lambda = lambda,
    control = tt_control(max_sweeps = 10L, pirls_maxit = 15L,
                         lbfgs_maxit = 150L, backend = "R",
                         compute_edf = FALSE)
  )
}

cat("\n=== Ishigami (d=3) ===\n")
dat_i <- simulate_ishigami(n = 500, sigma = 0.15, seed = 1)
print(dat_i)
fit_i <- .fit_demo(dat_i, rank = 2L, k = 8L)
rmse <- sqrt(mean((fit_i$linear.predictors - dat_i$f)^2))
cat(sprintf("  optimizer=%s | RMSE(eta,f)=%.3f | CR=%.1fx\n",
            fit_i$optimizer_used, rmse, fit_i$compression_ratio))

cat("\n=== Sobol g-function (d=4, a=(0,.5,3,9)) ===\n")
dat_s <- simulate_sobol_g(n = 500, d = 4, seed = 2)
print(dat_s)
fit_s <- .fit_demo(dat_s, rank = 2L, k = 6L)
rmse <- sqrt(mean((fit_s$linear.predictors - dat_s$f)^2))
cat(sprintf("  optimizer=%s | RMSE(eta,f)=%.3f | CR=%.1fx\n",
            fit_s$optimizer_used, rmse, fit_s$compression_ratio))

cat("\n=== Friedman #1 (d=5) ===\n")
dat_f <- simulate_friedman(n = 600, sigma = 1, seed = 3)
print(dat_f)
fit_f <- .fit_demo(dat_f, rank = 3L, k = 5L)
rmse <- sqrt(mean((fit_f$linear.predictors - dat_f$f)^2))
cat(sprintf("  optimizer=%s | RMSE(eta,f)=%.3f | CR=%.1fx\n",
            fit_f$optimizer_used, rmse, fit_f$compression_ratio))

cat("\n=== Packaged data(ishigami) ===\n")
utils::data("ishigami", package = "TTPsplines", envir = environment())
X <- as.matrix(ishigami[, c("x1", "x2", "x3")])
fit <- ttps(
  ishigami$y, X, rank = 2, k = 8, lambda = 1,
  control = tt_control(max_sweeps = 8, backend = "R", compute_edf = FALSE)
)
print(summary(fit))
cat("\nDone.\n")
