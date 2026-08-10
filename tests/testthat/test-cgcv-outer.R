test_that("damped trust update clips log10 steps and respects bounds", {
  step <- TTPsplines:::.cgcv_damped_trust_update(
    lambda_old = c(1, 1),
    lambda_tilde = c(1e6, 1e-6),
    rho = 1,
    max_log10_step = 1,
    bounds = c(1e-4, 1e4)
  )
  expect_equal(step$lambda_new[1], 10, tolerance = 1e-8)
  expect_equal(step$lambda_new[2], 0.1, tolerance = 1e-8)

  damp <- TTPsplines:::.cgcv_damped_trust_update(
    lambda_old = 1,
    lambda_tilde = 100,
    rho = 0.5,
    max_log10_step = Inf,
    bounds = c(1e-4, 1e4)
  )
  # log λ = 0.5 * log(100) = log(10)
  expect_equal(damp$lambda_new[[1]], 10, tolerance = 1e-8)
})

test_that("scale-anisotropy recenters prod omega = 1", {
  sa <- TTPsplines:::.cgcv_scale_anisotropy_from_lambda(c(2, 8, 0.5))
  expect_equal(prod(sa$omega), 1, tolerance = 1e-10)
  upd <- TTPsplines:::.cgcv_update_scale_anisotropy(
    lambda_old = c(1, 1, 1),
    lambda_tilde = c(10, 100, 0.1),
    rho = 0.5,
    max_log10_step = Inf,
    bounds = c(1e-6, 1e6)
  )
  expect_equal(prod(upd$omega_new), 1, tolerance = 1e-8)
})

test_that("outer simultaneous cGCV runs on small Gaussian", {
  skip_on_cran()
  set.seed(21)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2] + rnorm(n, 0, 0.25)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = "cGCV",
    control = tt_control(
      max_sweeps = 6,
      outer_maxit = 4,
      cgcv_fit_sweeps = 4,
      cgcv_update = "outer_simultaneous",
      cgcv_damping = 0.25,
      cgcv_max_log10_step = 1,
      compute_edf = FALSE,
      warn_lambda_boundary = FALSE,
      seed = 21
    )
  )
  expect_true(identical(fit$cgcv$update, "outer_simultaneous"))
  expect_true(!is.null(fit$cgcv$proposals))
  expect_true(all(is.finite(fit$lambda)))
  expect_true(diff(range(fitted(fit))) > 0.05)
})

test_that("sequential cGCV stores trace with P_other diagnostics", {
  skip_on_cran()
  set.seed(22)
  n <- 60
  X <- matrix(runif(n * 2), n, 2)
  y <- X[, 1] + rnorm(n, 0, 0.2)
  fit <- ttps(
    y, X, rank = 1, k = 5, lambda = "cGCV",
    control = tt_control(
      max_sweeps = 3,
      cgcv_update = "sequential",
      cgcv_trace = TRUE,
      compute_edf = FALSE,
      warn_lambda_boundary = FALSE,
      seed = 22
    )
  )
  tr <- fit$cgcv$trace
  expect_true(is.data.frame(tr))
  expect_true(all(c("P_other_op", "ed", "lambda_new", "boundary") %in% names(tr)))
})

test_that("tt_cgcv_frozen_curves returns per-margin grids", {
  skip_on_cran()
  set.seed(23)
  n <- 70
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 5, compute_edf = FALSE, seed = 23)
  )
  fr <- tt_cgcv_frozen_curves(
    fit, grid = exp(seq(log(1e-3), log(1e3), length.out = 21))
  )
  expect_equal(nrow(fr$proposals), 2L)
  expect_true(all(fr$curves$margin %in% 1:2))
  expect_true(all(is.finite(fr$curves$gcv) | is.infinite(fr$curves$gcv)))
})
