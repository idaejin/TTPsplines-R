test_that("tt_rank builds chains", {
  expect_equal(tt_rank(3, d = 4), c(1L, 3L, 3L, 3L, 1L))
  expect_equal(tt_rank(c(2, 4, 3), d = 4), c(1L, 2L, 4L, 3L, 1L))
})

test_that("gaussian fixed-lambda fits and predicts", {
  set.seed(2)
  n <- 300
  X <- matrix(runif(n * 3), n, 3)
  f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + X[, 3]
  y <- f + rnorm(n, 0, 0.3)
  fit <- ttpspline(
    y, X, family = gaussian(), rank = 2, k = 6, lambda = 1,
    control = tt_control(max_sweeps = 6, backend = "R", seed = 1)
  )
  expect_s3_class(fit, "ttpspline")
  expect_equal(length(fitted(fit)), n)
  eta <- predict(fit, X, type = "link")
  expect_equal(eta, fit$linear.predictors, tolerance = 1e-8)
  cx <- tt_complexity(fit)
  expect_equal(cx$npar_dense, 6^3)
  expect_true(cx$npar_tt < cx$npar_dense)
})

test_that("gaussian cGCV runs", {
  set.seed(3)
  n <- 250
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.4)
  fit <- ttpspline(
    y, X, family = gaussian(), rank = 2, k = 6, lambda = "cGCV",
    control = tt_control(max_sweeps = 5, backend = "R",
                         lambda_bounds = c(1e-2, 1e2), seed = 2)
  )
  expect_equal(fit$lambda_method, "cGCV")
  expect_true(all(is.finite(fit$lambda)))
})

test_that("poisson PIRLS finite", {
  set.seed(4)
  n <- 300
  X <- matrix(runif(n * 3), n, 3)
  eta <- 0.6 * sin(2 * pi * X[, 1]) + 0.4 * cos(2 * pi * X[, 2]) + log(2)
  y <- rpois(n, exp(eta))
  fit <- ttpspline(
    y, X, family = poisson(), rank = 2, k = 6, lambda = 1,
    control = tt_control(pirls_maxit = 10, als_sweeps_per_pirls = 2,
                         backend = "R", seed = 3)
  )
  expect_true(all(is.finite(fit$fitted.values)))
  expect_true(all(fit$fitted.values > 0))
  expect_true(is.finite(fit$deviance))
})

test_that("bernoulli probabilities in (0,1)", {
  set.seed(5)
  n <- 350
  X <- matrix(runif(n * 3), n, 3)
  eta <- 1.2 * sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2])
  y <- rbinom(n, 1, plogis(eta))
  fit <- ttpspline(
    y, X, family = binomial(), rank = 2, k = 6, lambda = 5,
    control = tt_control(pirls_maxit = 12, als_sweeps_per_pirls = 2,
                         backend = "R", seed = 4, damping = TRUE)
  )
  p <- predict(fit, type = "response")
  expect_true(all(p > 0 & p < 1))
  expect_true(is.finite(fit$deviance))
})

test_that("Paper 2 lambda hooks are reserved", {
  set.seed(1)
  X <- matrix(runif(50 * 2), 50, 2)
  y <- rnorm(50)
  expect_error(
    ttpspline(y, X, rank = 1, k = 5, lambda = "cFS",
              control = tt_control(max_sweeps = 2, backend = "R")),
    "Paper 2"
  )
})
