# Poisson example (auto -> PIRLS-ALS)
#   Rscript inst/examples/example_poisson.R

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

set.seed(2)
n <- 800
X <- matrix(runif(n * 3), n, 3)
f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.4 * X[, 3]
eta <- f - mean(f) + log(3)
y <- rpois(n, exp(eta))

fit <- ttpspline(
  y, X,
  family = poisson(),
  rank = 2,
  k = 8,
  lambda = 1,
  control = tt_control(
    pirls_maxit = 25, als_sweeps_per_pirls = 4,
    backend = "auto", compute_edf = FALSE
  )
)

stopifnot(identical(fit$optimizer_requested, "auto"))
stopifnot(identical(fit$optimizer_used, "PIRLS-ALS"))
stopifnot(identical(fit$optimizer_reason, "poisson family default"))

cat("Poisson TT-P-spline\n")
print(summary(fit))
cat("\npredict(type='response') [1:5]:\n")
print(predict(fit, X[1:5, ], type = "response"))
cat("\nDone.\n")
