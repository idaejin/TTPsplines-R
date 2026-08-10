test_that("0/1 weights match subset fit for Poisson", {
  set.seed(42)
  n <- 100
  X <- matrix(runif(n * 2), n, 2)
  f <- sin(2 * pi * X[, 1]) + 0.4 * X[, 2]
  y <- rpois(n, exp(f - mean(f) + log(1.5)))
  w <- as.numeric(runif(n) > 0.35)
  if (sum(w) < 20L) w[seq_len(20L)] <- 1

  ctrl <- tt_control(
    pirls_maxit = 15, max_sweeps = 6, backend = "R",
    compute_edf = FALSE, seed = 7L
  )
  # Shared knots so basis margins match (subset alone would rebuild knots)
  bs <- build_marginal_bases(X, k = 6, degree = 3)
  fit_w <- ttps(
    y, X, family = poisson(), rank = 2, k = 6, lambda = 1,
    weights = w, knots = bs$knots, control = ctrl
  )
  idx <- which(w > 0)
  fit_sub <- ttps(
    y[idx], X[idx, , drop = FALSE], family = poisson(),
    rank = 2, k = 6, lambda = 1, knots = bs$knots, control = ctrl
  )

  expect_equal(fit_w$weights, w)
  expect_equal(fit_w$deviance, fit_sub$deviance, tolerance = 1e-6)
  expect_equal(fit_w$intercept, fit_sub$intercept, tolerance = 1e-6)
  mu_w <- predict(fit_w, newdata = X[idx, , drop = FALSE], type = "response")
  mu_s <- fitted(fit_sub)
  expect_equal(mu_w, mu_s, tolerance = 1e-5)
})

test_that("Gaussian ALS respects observation weights", {
  set.seed(9)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  w <- rep(1, n)
  w[1:20] <- 0
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = 1,
    weights = w,
    control = tt_control(max_sweeps = 6, backend = "R", compute_edf = FALSE,
                         seed = 1L)
  )
  expect_equal(fit$weights, w)
  expect_true(is.finite(fit$deviance))
  # Deviance only from positive-weight rows
  expect_equal(
    fit$deviance,
    sum(w * (y - fitted(fit))^2),
    tolerance = 1e-8
  )
})

test_that("weights validate non-negative and positive sum", {
  set.seed(1)
  n <- 40
  X <- matrix(runif(n * 2), n, 2)
  y <- rnorm(n)
  expect_error(
    ttps(y, X, rank = 1, k = 4, lambda = 1, weights = -1,
              control = tt_control(max_sweeps = 2, compute_edf = FALSE)),
    "non-negative"
  )
  expect_error(
    ttps(y, X, rank = 1, k = 4, lambda = 1, weights = rep(0, n),
              control = tt_control(max_sweeps = 2, compute_edf = FALSE)),
    "sum\\(weights\\)"
  )
})

test_that("Rcpp backend accepts weights (Gaussian + Poisson 0/1)", {
  skip_if_not(exists("tt_fit_d_cpp", mode = "function"))
  set.seed(11)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  f <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2]
  w <- as.numeric(runif(n) > 0.3)
  if (sum(w) < 20L) w[seq_len(20L)] <- 1
  bs <- build_marginal_bases(X, k = 6, degree = 3)

  yg <- f + rnorm(n, 0, 0.25)
  g_r <- ttps(
    yg, X, rank = 2, k = 5, lambda = 1, weights = w, knots = bs$knots,
    control = tt_control(max_sweeps = 8, backend = "R", compute_edf = FALSE,
                         seed = 2L)
  )
  g_c <- ttps(
    yg, X, rank = 2, k = 5, lambda = 1, weights = w, knots = bs$knots,
    control = tt_control(max_sweeps = 8, backend = "Rcpp", compute_edf = FALSE,
                         seed = 2L)
  )
  expect_identical(g_c$backend, "Rcpp")
  expect_equal(g_c$deviance, g_r$deviance, tolerance = 1e-8)

  y <- rpois(n, exp(f - mean(f) + log(1.2)))
  idx <- which(w > 0)
  ctrl <- tt_control(
    pirls_maxit = 15, max_sweeps = 6, backend = "Rcpp",
    compute_edf = FALSE, seed = 7L
  )
  fit_w <- ttps(
    y, X, family = poisson(), rank = 2, k = 6, lambda = 1,
    weights = w, knots = bs$knots, control = ctrl
  )
  fit_sub <- ttps(
    y[idx], X[idx, , drop = FALSE], family = poisson(),
    rank = 2, k = 6, lambda = 1, knots = bs$knots, control = ctrl
  )
  expect_identical(fit_w$backend, "Rcpp")
  expect_equal(fit_w$deviance, fit_sub$deviance, tolerance = 1e-8)
  expect_equal(fit_w$intercept, fit_sub$intercept, tolerance = 1e-8)
})
