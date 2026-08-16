test_that("build_tt_dlnm lag windows and basis dims", {
  set.seed(1)
  n <- 60
  x <- list(temp = rnorm(n), pm10 = runif(n, 0, 50))
  dlnm <- build_tt_dlnm(x, lag = 5L, k = 6L, k_lag = 4L, degree = 3L)
  expect_equal(sum(dlnm$ok), n - 5L)
  expect_equal(length(dlnm$basis_lags), 6L)
  expect_equal(ncol(dlnm$basis_lags[[1]][[1]]), 6L)
  expect_equal(ncol(dlnm$basis_lags[[1]][["lag"]]), 4L)
  expect_equal(nrow(dlnm$basis_lags[[1]][[1]]), sum(dlnm$ok))
})

test_that("ttps_dlnm Gaussian runs and overall predict is centered", {
  set.seed(2)
  n <- 120
  temp <- sin(seq(0, 4 * pi, length.out = n)) + rnorm(n, 0, 0.1)
  pm10 <- runif(n, 10, 40)
  # mild concurrent effect
  y <- 2 + 0.3 * temp + 0.01 * pm10 + rnorm(n, 0, 0.2)
  year <- 2000 + (seq_len(n) - 1) %/% 30
  month <- ((seq_len(n) - 1) %% 12) + 1
  fit <- ttps_dlnm(
    y,
    list(temp = temp, pm10 = pm10),
    lag = 3L,
    family = gaussian(),
    rank = 1, k = 5, k_lag = 4,
    lambda = 1,
    smooth = list(
      year = list(x = year, bs = "ps", k = 6, m = 2, lambda = 10),
      month = list(x = month, bs = "cc", k = 6, m = 2,
                   period = c(0.5, 12.5), lambda = 10)
    ),
    control = tt_control(max_sweeps = 6, compute_edf = FALSE)
  )
  expect_s3_class(fit, "ttps_dlnm")
  expect_true(is.finite(fit$deviance))
  expect_equal(length(fit$lambda), 3L)
  at <- seq(-1, 1, length.out = 21)
  pr <- predict_dlnm(fit, var = "temp", at = at, cen = c(temp = 0, pm10 = 25),
                     type = "overall")
  expect_equal(nrow(pr), length(at))
  # RR at centering value ~ 1
  i0 <- which.min(abs(at - 0))
  expect_true(abs(pr$RR[i0] - 1) < 1e-6)
  sl <- predict_dlnm(fit, var = "temp", at = at, cen = c(temp = 0, pm10 = 25),
                     type = "slice", lag_at = 0L)
  expect_equal(sl$lag[1], 0)
})

test_that("ttps_dlnm Poisson with linear + year/month smooth", {
  set.seed(3)
  n <- 100
  temp <- rnorm(n)
  o3 <- runif(n, 5, 40)
  eta <- 0.5 + 0.2 * temp
  y <- rpois(n, exp(eta))
  dow <- factor(sample(0:6, n, replace = TRUE))
  lin <- model.matrix(~ 0 + dow)
  year <- rep(2001:2004, length.out = n)
  month <- rep(1:12, length.out = n)
  fit <- ttps_dlnm(
    y,
    list(temp = temp, o3 = o3),
    lag = 2L,
    family = poisson(),
    rank = 1, k = 5, k_lag = 4,
    lambda = 1,
    linear = lin,
    smooth = list(
      year = list(x = year, bs = "ps", k = 5, m = 2, target_edf = 3),
      month = list(x = month, bs = "cc", k = 6, m = 2,
                   period = c(0.5, 12.5), target_edf = 3)
    ),
    control = tt_control(
      pirls_maxit = 8, als_sweeps_per_pirls = 1,
      max_sweeps = 4, compute_edf = FALSE, tol = 1e-4
    )
  )
  expect_s3_class(fit, "ttps_dlnm")
  expect_true(fit$deviance > 0)
  expect_true(!is.null(fit$smooth$year))
  expect_true(!is.null(fit$smooth$month))
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Tensor-Train DLNM")
  s <- summary(fit)
  expect_s3_class(s, "summary.ttpspline")
  sout <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(sout, "Tensor-Train DLNM")
  expect_true(!is.null(s$coefficients))
})
