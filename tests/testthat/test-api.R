test_that("lambda scalar expands and anisotropic validates", {
  expect_equal(parse_lambda_spec(1, d = 3)$values, c(1, 1, 1))
  expect_equal(parse_lambda_spec(c(0.1, 1, 10), d = 3)$values, c(0.1, 1, 10))
  expect_equal(parse_lambda_spec(1, d = 3)$method, "fixed")
  expect_true(parse_lambda_spec("cGCV", d = 2)$automatic)
  expect_error(parse_lambda_spec(c(1, 2), d = 3), "length")
  expect_error(parse_lambda_spec(c(1, -1, 2), d = 3), "positive")
  expect_error(parse_lambda_spec(c(1, NA, 2), d = 3), "finite")
})

test_that("tt_initialize is reproducible and shared", {
  init1 <- tt_initialize(d = 3, rank = 2, k = 6, seed = 99)
  init2 <- tt_initialize(d = 3, rank = 2, k = 6, seed = 99)
  expect_equal(init1, init2)
  expect_equal(dim(init1[[1]]), c(1L, 6L, 2L))
  expect_equal(attr(init1, "ranks"), c(1L, 2L, 2L, 1L))
})

test_that("same init yields ALS and LBFGS eta close (fixed lambda)", {
  set.seed(11)
  n <- 220
  X <- matrix(runif(n * 3), n, 3)
  f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.5 * X[, 3]
  y <- f + rnorm(n, 0, 0.25)
  init <- tt_initialize(X, rank = 2, k = 5, seed = 7, sd = 0.1)
  ctrl <- tt_control(max_sweeps = 25, lbfgs_maxit = 200, backend = "R",
                     seed = 7, init_sd = 0.1, tol = 1e-8)
  fit_als <- ttps(y, X, rank = 2, k = 5, lambda = 1,
                       optimizer = "ALS", init = init, control = ctrl)
  fit_lbfgs <- ttps(y, X, rank = 2, k = 5, lambda = 1,
                         optimizer = "LBFGS", init = init, control = ctrl)
  # Compare fitted surfaces, not cores (gauge)
  rmse <- sqrt(mean((fit_als$linear.predictors - fit_lbfgs$linear.predictors)^2))
  expect_true(rmse < 0.35)
  expect_equal(fit_als$optimizer, "ALS")
  expect_equal(fit_lbfgs$optimizer, "LBFGS")
})

test_that("Adam fails with clear optional-backend message", {
  set.seed(1)
  X <- matrix(runif(40 * 2), 40, 2)
  y <- rnorm(40)
  expect_error(
    ttps(y, X, rank = 1, k = 4, lambda = 1, optimizer = "Adam",
              control = tt_control(backend = "R", max_sweeps = 2)),
    "Adam"
  )
  st <- tt_keras_status()
  expect_true(is.list(st))
  expect_false(isTRUE(tt_has_keras()) && FALSE) # smoke: callable
})

test_that("joint EDF is finite for small Gaussian fits", {
  set.seed(12)
  n <- 200
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.3)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 8, backend = "R", compute_edf = TRUE)
  )
  expect_true(is.finite(fit$edf))
  expect_true(fit$edf > 1 && fit$edf < fit$n + 1)
  expect_true(fit$edf <= fit$npar_tt + 1)
})

test_that("anisotropic fixed lambda accepted", {
  set.seed(8)
  n <- 180
  X <- matrix(runif(n * 3), n, 3)
  y <- rnorm(n)
  fit <- ttps(
    y, X, rank = 1, k = 5, lambda = c(0.5, 2, 1),
    control = tt_control(max_sweeps = 4, backend = "R")
  )
  expect_equal(fit$lambda, c(0.5, 2, 1), tolerance = 1e-12)
})
