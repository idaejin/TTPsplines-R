# Quick example: three families (not a full benchmark).
#   Rscript inst/examples/example_three_families.R

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

set.seed(1)
n <- 600
X <- matrix(runif(n * 3), n, 3)
f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.4 * X[, 3]

cat("\n--- Gaussian ---\n")
yg <- f + rnorm(n, 0, 0.3)
fit_g <- ttpspline(yg, X, family = gaussian(), rank = 3, k = 8, lambda = "cGCV",
                   control = tt_control(max_sweeps = 10, backend = "auto"))
print(tt_complexity(fit_g))
print(fit_g)

cat("\n--- Poisson ---\n")
yp <- rpois(n, exp(f - mean(f) + log(3)))
fit_p <- ttpspline(yp, X, family = poisson(), rank = 3, k = 8, lambda = 1,
                   control = tt_control(pirls_maxit = 12, backend = "auto"))
print(fit_p)

cat("\n--- Bernoulli ---\n")
yb <- rbinom(n, 1, plogis(1.5 * (f - mean(f))))
fit_b <- ttpspline(yb, X, family = binomial(), rank = 3, k = 8, lambda = 5,
                   control = tt_control(pirls_maxit = 12, backend = "auto"))
print(fit_b)

cat("\nDone.\n")
