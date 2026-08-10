# Conditional inference: Jacobian FD, gauge invariance, API

test_that("Gaussian Bayesian/frequentist vcov is finite, SPD-ish", {
  skip_on_cran()
  set.seed(11)
  n <- 120
  X <- matrix(runif(n * 2), n, 2)
  f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2])
  y <- f + rnorm(n, 0, 0.25)
  fit <- ttps(
    y, X, family = gaussian(), rank = 2, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 8, compute_edf = TRUE, seed = 11)
  )
  Vb <- vcov(fit, type = "bayesian")
  Vf <- vcov(fit, type = "frequentist")
  expect_true(is.matrix(Vb) && is.matrix(Vf))
  expect_equal(dim(Vb), dim(Vf))
  expect_true(all(is.finite(Vb)))
  expect_true(all(is.finite(Vf)))
  expect_lt(max(abs(Vb - t(Vb))), 1e-6)
  expect_lt(max(abs(Vf - t(Vf))), 1e-6)
  expect_true(all(diag(Vb) >= -1e-10))
  expect_true(all(diag(Vf) >= -1e-10))
  expect_error(vcov(fit, unconditional = TRUE), "smoothing-parameter")
})

test_that("predict se.fit and confidence intervals (Gaussian)", {
  skip_on_cran()
  set.seed(12)
  n <- 100
  X <- matrix(runif(n * 2), n, 2)
  y <- X[, 1] + 0.5 * X[, 2] + rnorm(n, 0, 0.2)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = 2,
    control = tt_control(max_sweeps = 6, compute_edf = FALSE, seed = 12)
  )
  pr <- predict(fit, X[1:10, ], se.fit = TRUE)
  expect_named(pr, c("fit", "se.fit"))
  expect_equal(length(pr$fit), 10L)
  expect_true(all(pr$se.fit > 0))
  ci <- predict(fit, X[1:10, ], interval = "confidence", level = 0.9)
  expect_true(all(ci$lower <= ci$fit))
  expect_true(all(ci$fit <= ci$upper))
  # Default predict still returns a vector
  eta <- predict(fit, X[1:5, ])
  expect_true(is.numeric(eta) && length(eta) == 5L)
})

test_that("prediction Jacobian matches finite differences", {
  skip_on_cran()
  set.seed(13)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 6, compute_edf = FALSE, seed = 13)
  )
  fit <- TTPsplines:::tt_prepare_inference(fit)
  cores <- fit$inference$cores_gauge
  x0 <- matrix(c(0.3, 0.7), 1, 2)
  b0 <- eval_marginal_bases(x0, fit$knots, fit$degree)
  J_an <- TTPsplines:::.tt_predict_jacobian(cores, TRUE, b0)
  J_fd <- TTPsplines:::.tt_predict_jacobian_fd(cores, fit$intercept, b0, eps = 1e-5)
  expect_equal(as.numeric(J_an), as.numeric(J_fd), tolerance = 5e-4)
})

test_that("prediction SE is gauge-invariant", {
  skip_on_cran()
  set.seed(14)
  n <- 100
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) * X[, 2] + rnorm(n, 0, 0.2)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 8, compute_edf = FALSE, seed = 14)
  )
  Xte <- matrix(runif(20 * 2), 20, 2)
  se0 <- predict(fit, Xte, se.fit = TRUE)$se.fit

  # Gauge-equivalent cores: transform interface 1
  A <- matrix(c(1.4, 0.3, 0.2, 0.9), 2, 2)
  fit2 <- fit
  fit2$cores <- TTPsplines:::tt_apply_gauge(fit$cores, iface = 1L, A = A)
  fit2$inference <- NULL
  fit2$._inf <- new.env(parent = emptyenv())
  # Same fitted values
  expect_equal(predict(fit, Xte), predict(fit2, Xte), tolerance = 1e-7)
  se1 <- predict(fit2, Xte, se.fit = TRUE)$se.fit
  expect_equal(se0, se1, tolerance = 2e-2)
  # After left-orthogonal gauge fixing, packed vcov may coincide; the
  # scientific requirement is invariance of SE(f̂), already checked above.
})

test_that("cGCV still uses conditional inference", {
  skip_on_cran()
  set.seed(15)
  n <- 120
  X <- matrix(runif(n * 2), n, 2)
  # Mildly noisy additive signal — cGCV should stay interior with default bounds
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2] + rnorm(n, 0, 0.4)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = "cGCV",
    control = tt_control(max_sweeps = 8, compute_edf = FALSE, seed = 15)
  )
  expect_identical(fit$lambda_method, "cGCV")
  pr <- predict(fit, X[1:5, ], se.fit = TRUE)
  expect_true(all(is.finite(pr$se.fit)))
  expect_true(isTRUE(fit$._inf$data$conditional_on_lambda))
  expect_false(isTRUE(fit$._inf$data$smoothing_uncertainty))
})

test_that("Poisson and Bernoulli pointwise intervals (link-scale SE)", {
  skip_on_cran()
  set.seed(16)
  n <- 100
  X <- matrix(runif(n * 2), n, 2)
  eta <- 0.8 * X[, 1] - 0.5 * X[, 2]
  yp <- rpois(n, exp(eta - mean(eta) + log(3)))
  fit_p <- ttps(
    yp, X, family = poisson(), rank = 2, k = 5, lambda = 1,
    control = tt_control(pirls_maxit = 15, compute_edf = FALSE, seed = 16)
  )
  pr <- predict(fit_p, X[1:8, ], type = "response", se.fit = TRUE,
                interval = "confidence")
  expect_true(all(pr$fit > 0))
  expect_true(all(pr$se.fit > 0))
  expect_true(all(pr$lower > 0))
  expect_true(all(pr$lower <= pr$fit & pr$fit <= pr$upper))
  # se.fit remains on link scale even when type=response
  expect_equal(names(pr)[1:2], c("fit", "se.fit"))
  Vb <- vcov(fit_p)
  expect_true(all(is.finite(Vb)))
  expect_equal(fit_p$._inf$data$scale, 1)

  yb <- rbinom(n, 1, plogis(1.2 * (eta - mean(eta))))
  fit_b <- ttps(
    yb, X, family = binomial(), rank = 2, k = 5, lambda = 2,
    control = tt_control(lbfgs_maxit = 80, compute_edf = FALSE, seed = 17)
  )
  cib <- predict(fit_b, X[1:5, ], type = "response", interval = "confidence")
  expect_true(all(cib$lower >= 0 & cib$upper <= 1))
  expect_true(all(cib$lower <= cib$fit & cib$fit <= cib$upper))
})

test_that("Poisson prediction SE is gauge-invariant", {
  skip_on_cran()
  set.seed(18)
  n <- 90
  X <- matrix(runif(n * 2), n, 2)
  eta <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2]
  y <- rpois(n, exp(eta - mean(eta) + log(2)))
  fit <- ttps(
    y, X, family = poisson(), rank = 2, k = 5, lambda = 1,
    control = tt_control(pirls_maxit = 12, compute_edf = FALSE, seed = 18)
  )
  Xte <- matrix(runif(15 * 2), 15, 2)
  se0 <- predict(fit, Xte, se.fit = TRUE)$se.fit
  A <- matrix(c(1.3, 0.2, 0.1, 0.95), 2, 2)
  fit2 <- fit
  fit2$cores <- TTPsplines:::tt_apply_gauge(fit$cores, iface = 1L, A = A)
  fit2$inference <- NULL
  fit2$._inf <- new.env(parent = emptyenv())
  expect_equal(predict(fit, Xte), predict(fit2, Xte), tolerance = 1e-6)
  se1 <- predict(fit2, Xte, se.fit = TRUE)$se.fit
  expect_equal(se0, se1, tolerance = 2e-2)
})
