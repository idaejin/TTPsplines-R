# Quick example: three families (not a full benchmark).
#   Rscript inst/examples/example_three_families.R
#
# Demonstrates family-aware optimizer = "auto":
#   Gaussian  -> ALS
#   Poisson   -> PIRLS-ALS
#   binomial  -> LBFGS

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

.show_opt <- function(fit) {
  cat(sprintf(
    "  optimizer: requested=%s | used=%s | reason=%s\n",
    fit$optimizer_requested, fit$optimizer_used, fit$optimizer_reason
  ))
}

cat("\n--- Gaussian (auto -> ALS) ---\n")
yg <- f + rnorm(n, 0, 0.3)
fit_g <- ttpspline(
  yg, X, family = gaussian(), rank = 2, k = 8, lambda = "cGCV",
  control = tt_control(max_sweeps = 10, backend = "auto", compute_edf = FALSE)
)
.stopifnot_opt <- function(fit, used) {
  stopifnot(identical(fit$optimizer_requested, "auto"))
  stopifnot(identical(fit$optimizer_used, used))
}
.stopifnot_opt(fit_g, "ALS")
.show_opt(fit_g)
print(tt_complexity(fit_g))
print(summary(fit_g))

cat("\n--- Poisson (auto -> PIRLS-ALS) ---\n")
yp <- rpois(n, exp(f - mean(f) + log(3)))
fit_p <- ttpspline(
  yp, X, family = poisson(), rank = 2, k = 8, lambda = 1,
  control = tt_control(pirls_maxit = 20, als_sweeps_per_pirls = 3,
                       backend = "auto", compute_edf = FALSE)
)
.stopifnot_opt(fit_p, "PIRLS-ALS")
.show_opt(fit_p)
print(summary(fit_p))

cat("\n--- Bernoulli (auto -> LBFGS) ---\n")
yb <- rbinom(n, 1, plogis(1.5 * (f - mean(f))))
fit_b <- ttpspline(
  yb, X, family = binomial(), rank = 2, k = 8, lambda = 5,
  control = tt_control(lbfgs_maxit = 200, backend = "auto", compute_edf = FALSE)
)
.stopifnot_opt(fit_b, "LBFGS")
.show_opt(fit_b)
print(summary(fit_b))

cat("\nDone.\n")
