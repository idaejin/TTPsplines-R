test_that("auto resolves family-aware defaults transparently", {
  set.seed(41)
  n <- 120
  X <- matrix(runif(n * 2), n, 2)
  ctrl <- tt_control(backend = "R", max_sweeps = 4L, pirls_maxit = 8L,
                     als_sweeps_per_pirls = 2L, lbfgs_maxit = 40L,
                     compute_edf = FALSE)

  y_g <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  fit_g <- ttps(y_g, X, family = gaussian(), rank = 2, k = 5,
                     lambda = 1, optimizer = "auto", control = ctrl)
  expect_equal(fit_g$optimizer_requested, "auto")
  expect_equal(fit_g$optimizer_used, "ALS")
  expect_equal(fit_g$optimizer, "ALS")
  expect_equal(fit_g$optimizer_reason, "gaussian family default")

  y_p <- rpois(n, exp(0.5 * sin(2 * pi * X[, 1])))
  fit_p <- ttps(y_p, X, family = poisson(), rank = 2, k = 5,
                     lambda = 1, optimizer = "auto", control = ctrl)
  expect_equal(fit_p$optimizer_requested, "auto")
  expect_equal(fit_p$optimizer_used, "PIRLS-ALS")
  expect_equal(fit_p$optimizer, "PIRLS-ALS")
  expect_equal(fit_p$optimizer_reason, "poisson family default")

  y_b <- rbinom(n, 1, plogis(sin(2 * pi * X[, 1])))
  fit_b <- ttps(y_b, X, family = binomial(), rank = 2, k = 5,
                     lambda = 1, optimizer = "auto", control = ctrl)
  expect_equal(fit_b$optimizer_requested, "auto")
  expect_equal(fit_b$optimizer_used, "LBFGS")
  expect_equal(fit_b$optimizer, "LBFGS")
  expect_equal(fit_b$optimizer_reason, "binomial family default")

  out <- capture.output(summary(fit_b))
  expect_true(any(grepl("Requested optimizer:\\s+auto", out)))
  expect_true(any(grepl("Selected optimizer:\\s+LBFGS", out)))
  expect_true(any(grepl("Reason:\\s+binomial family default", out)))
})

test_that("manual optimizer override is recorded as user-specified", {
  set.seed(42)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  y <- rbinom(n, 1, 0.4)
  fit <- ttps(
    y, X, family = binomial(), rank = 2, k = 4, lambda = 2,
    optimizer = "PIRLS-ALS",
    control = tt_control(backend = "R", pirls_maxit = 10L,
                         als_sweeps_per_pirls = 2L, compute_edf = FALSE)
  )
  expect_equal(fit$optimizer_requested, "PIRLS-ALS")
  expect_equal(fit$optimizer_used, "PIRLS-ALS")
  expect_equal(fit$optimizer_reason, "user-specified")
})
