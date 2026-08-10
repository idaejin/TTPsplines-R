test_that("ttps_multistart returns best fit and start table", {
  set.seed(7)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + 0.3 * cos(2 * pi * X[, 2]) + rnorm(n, 0, 0.2)
  ms <- ttps_multistart(
    y, X,
    rank = 2, k = 5, lambda = 1,
    n_starts = 3L, seed = 7L,
    control = tt_control(max_sweeps = 4L, compute_edf = FALSE, backend = "R")
  )
  expect_s3_class(ms, "ttps_multistart")
  expect_s3_class(ms$best, "ttpspline")
  expect_equal(nrow(ms$starts), 3L)
  expect_equal(ms$n_starts, 3L)
  expect_true(ms$selected_start %in% 1:3)
  expect_true(is.finite(ms$best$deviance))
  # best objective is minimal among finite starts (prefer_converged)
  objs <- ms$starts$objective
  ok <- is.finite(objs) & ms$starts$converged
  if (any(ok)) {
    expect_equal(objs[ms$selected_start], min(objs[ok]))
  }
})

test_that("ttps_multistart cGCV records boundary fractions", {
  set.seed(8)
  n <- 60
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.25)
  ms <- ttps_multistart(
    y, X,
    rank = 2, k = 5, lambda = "cGCV",
    n_starts = 2L, seed = 8L,
    control = tt_control(
      max_sweeps = 4L, compute_edf = FALSE, backend = "R",
      lambda_bounds = c(1e-4, 1e4), warn_lambda_boundary = FALSE
    )
  )
  expect_length(ms$frac_boundary, 2L)
  expect_true(all(names(ms$frac_boundary) == c("margin1", "margin2")))
  expect_true(all(is.na(ms$frac_boundary) | (ms$frac_boundary >= 0 & ms$frac_boundary <= 1)))
})

test_that("ttps_multistart is reproducible for fixed seed", {
  set.seed(9)
  n <- 50
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  ctrl <- tt_control(max_sweeps = 3L, compute_edf = FALSE, backend = "R")
  a <- ttps_multistart(y, X, rank = 2, k = 5, lambda = 1,
                       n_starts = 2L, seed = 21L, control = ctrl)
  b <- ttps_multistart(y, X, rank = 2, k = 5, lambda = 1,
                       n_starts = 2L, seed = 21L, control = ctrl)
  expect_equal(a$starts$objective, b$starts$objective, tolerance = 1e-8)
  expect_equal(a$best$deviance, b$best$deviance, tolerance = 1e-8)
})
