# Bernoulli example (auto -> LBFGS)
#   Rscript inst/examples/example_bernoulli.R

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

set.seed(3)
n <- 800
X <- matrix(runif(n * 3), n, 3)
f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.4 * X[, 3]
eta <- 1.5 * (f - mean(f))
y <- rbinom(n, 1, plogis(eta))

fit <- ttps(
  y, X,
  family = binomial(),
  rank = 2,
  k = 8,
  lambda = 5,
  control = tt_control(
    lbfgs_maxit = 300, backend = "auto", compute_edf = FALSE
  )
)

stopifnot(identical(fit$optimizer_requested, "auto"))
stopifnot(identical(fit$optimizer_used, "LBFGS"))
stopifnot(identical(fit$optimizer_reason, "binomial family default"))
p <- predict(fit, type = "response")
stopifnot(all(p > 0 & p < 1))

cat("Bernoulli TT-P-spline\n")
print(summary(fit))
cat("\npredict(type='response') [1:5]:\n")
print(predict(fit, X[1:5, ], type = "response"))
cat("\nDone.\n")
