test_that("glam_fit_poisson recovers a smooth age-year surface", {
  dat <- simulate_glam_poisson(n_age = 25L, n_year = 21L, seed = 7L)
  bb <- glam_grid_bases(list(age = dat$age, year = dat$year), k = 8)
  fit <- glam_fit_poisson(
    dat$Y, bb$B,
    lambda = c(5, 1),
    offset = log(dat$exposure),
    pirls_maxit = 30L
  )
  expect_equal(fit$family, "poisson")
  expect_true(fit$n_pirls >= 2L)
  expect_true(is.finite(fit$deviance))
  rmse_log <- sqrt(mean((log(as.numeric(fit$mu)) - log(as.numeric(dat$mu)))^2))
  expect_lt(rmse_log, 0.35)
})

test_that("glam_grid_bases matches grid dimensions", {
  axes <- list(u = seq(0, 1, length.out = 12), v = seq(-1, 1, length.out = 9))
  bb <- glam_grid_bases(axes, k = 6)
  expect_equal(nrow(bb$B[[1]]), 12L)
  expect_equal(nrow(bb$B[[2]]), 9L)
  expect_equal(ncol(bb$B[[1]]), 6L)
})

test_that("glam_fit_gaussian still works", {
  set.seed(1)
  n1 <- 15L; n2 <- 12L
  Y <- outer(seq_len(n1), seq_len(n2), function(i, j) sin(i / 3) + cos(j / 4))
  Y <- Y + rnorm(length(Y), sd = 0.05)
  dim(Y) <- c(n1, n2)
  bb <- glam_grid_bases(list(seq_len(n1), seq_len(n2)), k = 6)
  fit <- glam_fit_gaussian(Y, bb$B, lambda = 1)
  expect_equal(fit$npar, 36L)
  expect_lt(sqrt(mean((fit$mu - Y)^2)), 0.25)
})
