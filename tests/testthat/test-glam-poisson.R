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

test_that("glam_fit_poisson fit_weights = 0 cells do not enter the fit", {
  set.seed(31)
  axes <- list(u = seq(0, 1, length.out = 10), v = seq(0, 1, length.out = 9))
  bb <- glam_grid_bases(axes, k = 5)
  eta <- outer(axes$u, axes$v, function(a, b) log(0.5) + sin(2 * pi * a) + cos(2 * pi * b))
  Y <- array(rpois(90, exp(eta)), c(10, 9))
  W <- array(1, c(10, 9))
  W[1:3, 1:3] <- 0
  fit1 <- glam_fit_poisson(Y, bb$B, lambda = 1, fit_weights = W,
                           pirls_maxit = 100L, tol = 1e-14)
  Y2 <- Y
  Y2[W == 0] <- Y2[W == 0] + 25L       # change only held-out cells
  fit2 <- glam_fit_poisson(Y2, bb$B, lambda = 1, fit_weights = W,
                           pirls_maxit = 100L, tol = 1e-14)
  expect_equal(as.numeric(fit2$mu), as.numeric(fit1$mu), tolerance = 1e-8)
  # weighted penalized MLE by dense PIRLS
  B <- kronecker(bb$B[[2]], bb$B[[1]])
  P <- glam_penalty(c(5L, 5L), c(1, 1))
  y <- as.numeric(Y); wo <- as.numeric(W)
  e <- rep(log(mean(y[wo > 0])), length(y))
  for (it in 1:200) {
    mu <- exp(e)
    th <- solve(crossprod(B, wo * mu * B) + P, crossprod(B, wo * mu * (e + (y - mu) / mu)))
    e <- as.numeric(B %*% th)
  }
  expect_equal(as.numeric(fit1$mu), exp(e), tolerance = 1e-6)
  expect_equal(fit1$deviance, glm_deviance(poisson(), y, exp(e), weights = wo),
               tolerance = 1e-6)
})
