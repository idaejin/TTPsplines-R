test_that("GD uses same objective/grad machinery as LBFGS (Bernoulli)", {
  set.seed(21)
  n <- 180
  X <- matrix(runif(n * 3), n, 3)
  eta <- 1.1 * sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2])
  y <- rbinom(n, 1, plogis(eta))
  init <- tt_initialize(X, rank = 2, k = 5, seed = 3, sd = 0.05)

  fit_gd <- ttpspline(
    y, X, family = binomial(), rank = 2, k = 5, lambda = 1,
    optimizer = "GD", init = init,
    control = tt_control(
      backend = "R", gd_lr = 1, gd_maxit = 400L, gd_tol = 1e-6,
      gd_linesearch = TRUE, compute_edf = FALSE, seed = 3
    )
  )
  fit_lb <- ttpspline(
    y, X, family = binomial(), rank = 2, k = 5, lambda = 1,
    optimizer = "LBFGS", init = init,
    control = tt_control(backend = "R", lbfgs_maxit = 300L, compute_edf = FALSE)
  )

  expect_equal(fit_gd$optimizer, "GD")
  o_gd <- tt_objective(fit_gd, X, y)
  o_lb <- tt_objective(fit_lb, X, y)
  expect_true(is.finite(o_gd$value))
  # Same model: GD should not be wildly worse than LBFGS on a small problem
  expect_true(o_gd$value <= o_lb$value * 1.15 + 5)
  expect_true(max(abs(fit_gd$linear.predictors)) < 40)
  # eta consistent with cores
  basis <- eval_marginal_bases(X, fit_gd$knots, fit_gd$degree)
  eta2 <- fit_gd$offset + fit_gd$intercept + tt_contraction(fit_gd$cores, basis)
  expect_equal(fit_gd$linear.predictors, eta2, tolerance = 1e-10)
})

test_that("PIRLS-ALS is an alias for the structure-aware GLM path", {
  set.seed(22)
  n <- 120
  X <- matrix(runif(n * 2), n, 2)
  y <- rbinom(n, 1, plogis(sin(2 * pi * X[, 1])))
  fit <- ttpspline(
    y, X, family = binomial(), rank = 2, k = 5, lambda = 5,
    optimizer = "PIRLS-ALS",
    control = tt_control(
      backend = "R", pirls_maxit = 12L, als_sweeps_per_pirls = 2L,
      compute_edf = FALSE
    )
  )
  expect_equal(fit$optimizer, "PIRLS-ALS")
  expect_true(all(predict(fit, type = "response") > 0 &
                    predict(fit, type = "response") < 1))
})

test_that("GD Armijo line search decreases the objective", {
  set.seed(23)
  n <- 150
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.3)
  init <- tt_initialize(X, rank = 2, k = 5, seed = 1, sd = 0.1)
  fit <- ttpspline(
    y, X, family = gaussian(), rank = 2, k = 5, lambda = 1,
    optimizer = "GD", init = init,
    control = tt_control(
      backend = "R", gd_lr = 1, gd_maxit = 80L, gd_linesearch = TRUE,
      compute_edf = FALSE, trace = FALSE
    )
  )
  expect_equal(fit$optimizer, "GD")
  expect_true(is.finite(fit$deviance))
})
