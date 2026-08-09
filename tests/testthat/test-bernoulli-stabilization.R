test_that("Bernoulli W and z share floored working variance", {
  ctrl <- tt_control(binomial_mu_eps = 1e-5, binomial_weight_floor = 1e-4)
  eta <- c(-40, -20, 0, 20, 40)
  y <- c(0, 0, 1, 1, 1)
  w <- glm_working(binomial(), y, eta, control = ctrl)
  expect_equal(w$weight, w$var_work)
  expect_equal(w$z, eta + (y - w$mu) / w$var_work)
  # flooring active at extremes
  expect_true(any(w$var_raw < ctrl$binomial_weight_floor))
  expect_true(all(w$var_work >= ctrl$binomial_weight_floor - 1e-15))
})

test_that("Bernoulli PIRLS accepted objectives are non-increasing", {
  set.seed(11)
  n <- 200
  X <- matrix(runif(n * 3), n, 3)
  eta <- 1.2 * sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2])
  y <- rbinom(n, 1, plogis(eta))
  init <- tt_initialize(X, rank = 2, k = 5, seed = 2, sd = 0.05)
  fit <- ttpspline(
    y, X, family = binomial(), rank = 2, k = 5, lambda = 1,
    optimizer = "ALS", init = init,
    control = tt_control(
      backend = "R", pirls_maxit = 25L, als_sweeps_per_pirls = 3L,
      pirls_step_halving = TRUE, compute_edf = FALSE, seed = 2
    )
  )
  h <- fit$history
  expect_true(is.data.frame(h) && "objective" %in% names(h))
  if (nrow(h) >= 2L) {
    dL <- diff(h$objective[h$line_ok %||% TRUE])
    expect_true(all(dL <= 1e-8 + 1e-10))
  }
  # eta consistent with cores
  basis <- eval_marginal_bases(X, fit$knots, fit$degree)
  eta2 <- fit$intercept + tt_contraction(fit$cores, basis)
  expect_equal(fit$linear.predictors, eta2, tolerance = 1e-10)
  expect_true(!is.null(fit$convergence$pirls))
})

test_that("Gaussian and Poisson still fit after Bernoulli PIRLS changes", {
  set.seed(12)
  n <- 180
  X <- matrix(runif(n * 3), n, 3)
  f <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2]
  yg <- f + rnorm(n, 0, 0.25)
  yp <- rpois(n, exp(0.4 * f + log(2)))
  fg <- ttpspline(yg, X, family = gaussian(), rank = 2, k = 5, lambda = 1,
                  control = tt_control(backend = "R", max_sweeps = 8, compute_edf = FALSE))
  fp <- ttpspline(yp, X, family = poisson(), rank = 2, k = 5, lambda = 1,
                  control = tt_control(backend = "R", pirls_maxit = 10,
                                       als_sweeps_per_pirls = 2, compute_edf = FALSE))
  expect_true(is.finite(fg$deviance))
  expect_true(all(is.finite(fp$fitted.values)))
  expect_true(all(fp$fitted.values > 0))
})
