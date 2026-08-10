test_that("lambda boundary classifier detects search limits", {
  b <- c(1e-4, 1e4)
  expect_equal(
    .lambda_boundary_status(c(9922.5, 0.043, 0.000100003), b),
    c("upper", "interior", "lower")
  )
  expect_equal(
    .lambda_boundary_status(c(1, 10, 100), b),
    c("interior", "interior", "interior")
  )
  expect_equal(.lambda_boundary_status(1e-4, b), "lower")
  expect_equal(.lambda_boundary_status(1e4, b), "upper")
})

test_that("cGCV fit stores lambda_boundary and summary labels cumulative iters", {
  set.seed(9)
  n <- 150
  X <- matrix(runif(n * 3), n, 3)
  # Strong directional structure tends to push anisotropic cGCV to bounds
  f <- 3 * sin(2 * pi * X[, 1]) + 0.05 * X[, 2] + 0.01 * X[, 3]
  y <- f + rnorm(n, 0, 0.2)
  fit <- suppressWarnings(ttps(
    y, X, family = gaussian(), rank = 2, k = 6, lambda = "cGCV",
    control = tt_control(
      max_sweeps = 8L, backend = "R", compute_edf = FALSE,
      lambda_bounds = c(1e-4, 1e4), warn_lambda_boundary = FALSE
    )
  ))
  expect_true(!is.null(fit$lambda_boundary))
  expect_equal(length(fit$lambda_boundary), 3L)
  expect_true(all(fit$lambda_boundary %in% c("lower", "upper", "interior")))
  expect_equal(fit$lambda_bounds, c(1e-4, 1e4), tolerance = 1e-15)

  out <- capture.output(summary(fit))
  expect_true(any(grepl("Lambda boundary:", out, fixed = TRUE)))
  expect_true(any(grepl("Lambda search bounds:", out, fixed = TRUE)))
})

test_that("boundary hit emits soft warning when enabled", {
  expect_warning(
    .tt_lambda_boundary_info(
      c(9999, 1, 1e-4), "cGCV",
      list(lambda_bounds = c(1e-4, 1e4), warn_lambda_boundary = TRUE)
    ),
    "near the cGCV search boundaries"
  )
  expect_silent(
    .tt_lambda_boundary_info(
      c(9999, 1, 1e-4), "cGCV",
      list(lambda_bounds = c(1e-4, 1e4), warn_lambda_boundary = FALSE)
    )
  )
  # fixed λ: no warning even at bounds
  expect_silent(
    .tt_lambda_boundary_info(
      c(1e-4, 1, 1e4), "fixed",
      list(lambda_bounds = c(1e-4, 1e4), warn_lambda_boundary = TRUE)
    )
  )
})

test_that("summary marks optimizer iterations as cumulative when present", {
  set.seed(11)
  n <- 100
  X <- matrix(runif(n * 2), n, 2)
  y <- rbinom(n, 1, 0.4)
  fit <- ttps(
    y, X, family = binomial(), rank = 2, k = 4, lambda = 1,
    control = tt_control(lbfgs_maxit = 40L, backend = "R", compute_edf = FALSE)
  )
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Optimizer iterations:.*cumulative", out)))
})

test_that("summary does not claim lambda evaluations under fixed λ", {
  set.seed(3)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  fit <- ttps(
    y, X, family = gaussian(), rank = 2, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 4L, backend = "R", compute_edf = FALSE)
  )
  expect_identical(fit$lambda_method, "fixed")
  expect_equal(fit$n_criterion_evals, 0L)
  out <- capture.output(summary(fit))
  expect_false(any(grepl("^Lambda evaluations:", out)))
  expect_true(any(grepl("Criterion evaluations:.*fixed lambda", out)))
})
