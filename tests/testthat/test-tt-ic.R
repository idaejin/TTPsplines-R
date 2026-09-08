test_that("tt_ic Gaussian AIC/BIC from EDF", {
  set.seed(11)
  n <- 120
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2] + rnorm(n, 0, 0.25)
  fit <- ttps(
    y, X, rank = 1, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 4, backend = "R", compute_edf = TRUE)
  )
  expect_true(is.finite(fit$edf))
  aic <- tt_ic(fit, "AIC")
  bic <- tt_ic(fit, "BIC")
  expect_true(is.finite(aic) && is.finite(bic))
  # BIC penalty >= AIC for n >= e^2
  expect_gte(bic, aic)
  # manual Gaussian check
  df <- fit$edf + 1
  expect_equal(aic, n * log(fit$deviance / n) + 2 * df, tolerance = 1e-10)
  expect_equal(bic, n * log(fit$deviance / n) + log(n) * df, tolerance = 1e-10)
})

test_that("tt_ic errors without edf", {
  set.seed(12)
  X <- matrix(runif(40 * 2), 40, 2)
  y <- rnorm(40)
  fit <- ttps(
    y, X, rank = 1, k = 4, lambda = 1,
    control = tt_control(max_sweeps = 2, backend = "R", compute_edf = FALSE)
  )
  fit$edf <- NA_real_
  expect_error(tt_ic(fit, "AIC"), "edf")
})
