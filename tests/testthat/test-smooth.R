test_that("build_smooth_term ps and cc work", {
  set.seed(1)
  n <- 80
  x <- sort(runif(n))
  sm <- build_smooth_term(x, bs = "ps", k = 8, m = 2, name = "t")
  expect_equal(nrow(sm$B), n)
  expect_equal(ncol(sm$B), 8L)
  expect_equal(dim(sm$S), c(8L, 8L))
  expect_equal(sm$m, 2L)

  xc <- (seq_len(n) - 1) / n
  smc <- build_smooth_term(xc, bs = "cc", k = 10, m = 2, period = c(0, 1),
                           name = "cyc")
  expect_equal(ncol(smc$B), 10L)
  expect_equal(smc$bs, "cc")
})

test_that("Gaussian ttps recovers additive smooth of time", {
  set.seed(10)
  n <- 150
  X <- matrix(runif(n * 2), n, 2)
  colnames(X) <- c("x1", "x2")
  t <- seq(0, 1, length.out = n)
  f_tt <- 0.4 * sin(2 * pi * X[, 1])
  f_sm <- sin(2 * pi * t)
  y <- 1 + f_tt + f_sm + rnorm(n, 0, 0.15)
  fit <- ttps(
    y, X,
    rank = 1, k = 6, lambda = 1,
    smooth = list(time = list(x = t, bs = "ps", k = 12, m = 2, lambda = "cGCV")),
    control = tt_control(max_sweeps = 10, compute_edf = FALSE, outer_maxit = 4)
  )
  expect_true(!is.null(fit$smooth))
  expect_true(is.finite(fit$smooth$time$lambda))
  expect_true(is.finite(fit$smooth$time$edf))
  expect_true(fit$smooth$time$edf > 1)
  expect_true(fit$smooth$time$edf < 12)
  s <- summary(fit)
  expect_true(!is.null(s$smooth_table))
  expect_true("s(time)" %in% rownames(s$smooth_table))
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "Smooth terms:")
  # predict needs smooth newdata
  mu <- predict(fit, X, smooth = list(time = t))
  expect_equal(length(mu), n)
  expect_true(cor(mu, y) > 0.7)
  expect_error(predict(fit, X), "smooth")
})

test_that("Poisson ttps with dow linear + s(time) and custom m", {
  set.seed(11)
  n <- 120
  X <- matrix(runif(n * 2), n, 2)
  t <- seq(0, 1, length.out = n)
  dow <- factor(sample(0:6, n, replace = TRUE))
  L <- model.matrix(~ 0 + dow)
  eta <- 0.2 + 0.5 * sin(2 * pi * t) + 0.3 * X[, 1]
  y <- rpois(n, exp(eta))
  fit <- ttps(
    y, X, family = poisson(),
    rank = 1, k = 5, lambda = 1,
    linear = L,
    smooth = list(time = list(x = t, bs = "ps", k = 10, m = 3, lambda = 1)),
    control = tt_control(
      pirls_maxit = 10, outer_maxit = 1, als_sweeps_adaptive = FALSE,
      compute_edf = FALSE
    )
  )
  expect_equal(fit$smooth$time$m, 3L)
  expect_equal(fit$smooth$time$lambda_method, "fixed")
  expect_true(length(fit$beta) == ncol(L))
  s <- summary(fit)
  expect_true("z value" %in% colnames(s$coefficients))
  expect_equal(s$smooth_table["s(time)", "m"], 3)
})

test_that("cyclic smooth bs=cc runs", {
  set.seed(12)
  n <- 100
  X <- matrix(runif(n * 2), n, 2)
  hour <- (seq_len(n) - 1) %% 24
  y <- sin(2 * pi * hour / 24) + rnorm(n, 0, 0.2)
  fit <- ttps(
    y, X, rank = 1, k = 5, lambda = 1,
    smooth = list(hour = list(x = hour, bs = "cc", k = 8, m = 2,
                              period = c(0, 24), lambda = 1)),
    control = tt_control(max_sweeps = 6, compute_edf = FALSE)
  )
  expect_equal(fit$smooth$hour$bs, "cc")
  mu <- predict(fit, X, smooth = list(hour = hour))
  expect_true(cor(mu, y) > 0.5)
})

test_that("target_edf selects lambda with matching edf", {
  set.seed(13)
  n <- 160
  X <- matrix(runif(n * 2), n, 2)
  colnames(X) <- c("x1", "x2")
  t <- seq(0, 1, length.out = n)
  y <- 1 + 0.3 * sin(2 * pi * X[, 1]) + sin(4 * pi * t) + rnorm(n, 0, 0.2)
  fit <- ttps(
    y, X,
    rank = 1, k = 6, lambda = 1,
    smooth = list(time = list(x = t, bs = "ps", k = 40, m = 2, target_edf = 12)),
    control = tt_control(max_sweeps = 8, compute_edf = FALSE, outer_maxit = 3)
  )
  expect_equal(fit$smooth$time$lambda_method, "target_edf")
  expect_equal(fit$smooth$time$target_edf, 12)
  expect_true(abs(fit$smooth$time$edf - 12) < 0.5)
  expect_true(is.finite(fit$smooth$time$lambda) && fit$smooth$time$lambda > 0)
  s <- summary(fit)
  expect_equal(s$smooth_table["s(time)", "target_edf"], 12)
  expect_equal(as.character(s$smooth_table["s(time)", "method"]), "target_edf")
})

test_that("target_edf rejects invalid targets", {
  x <- seq(0, 1, length.out = 50)
  expect_error(
    build_smooth_term(x, bs = "ps", k = 10, m = 2, target_edf = 2),
    "target_edf"
  )
  expect_error(
    build_smooth_term(x, bs = "ps", k = 10, m = 2, target_edf = 11),
    "target_edf"
  )
})

test_that(".lambda_for_target_edf is monotone in target", {
  set.seed(14)
  n <- 100
  x <- seq(0, 1, length.out = n)
  sm <- build_smooth_term(x, bs = "ps", k = 20, m = 2, target_edf = 8)
  Gram <- crossprod(sm$B)
  h_hi <- .lambda_for_target_edf(Gram, sm$S, target = 14, bounds = c(1e-4, 1e6))
  h_lo <- .lambda_for_target_edf(Gram, sm$S, target = 6, bounds = c(1e-4, 1e6))
  expect_true(h_lo$lambda > h_hi$lambda)
  expect_true(abs(h_hi$edf - 14) < 0.3)
  expect_true(abs(h_lo$edf - 6) < 0.3)
})
