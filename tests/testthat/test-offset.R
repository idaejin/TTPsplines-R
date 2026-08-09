test_that("Poisson offset is used in eta and predict", {
  set.seed(21)
  n <- 120
  X <- matrix(runif(n * 2), n, 2)
  exposure <- runif(n, 50, 200)
  f <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2]
  mu <- exposure * exp(f - mean(f) + log(0.05))
  y <- rpois(n, mu)
  off <- log(exposure)

  fit <- ttpspline(
    y, X, family = poisson(), rank = 2, k = 6, lambda = 1,
    offset = off,
    control = tt_control(pirls_maxit = 20, backend = "R", compute_edf = FALSE)
  )
  expect_equal(fit$offset, off)
  basis <- eval_marginal_bases(X, fit$knots, fit$degree)
  expect_equal(
    fit$linear.predictors,
    off + fit$intercept + tt_contraction(fit$cores, basis),
    tolerance = 1e-8
  )
  # newdata without offset => smooth only + intercept
  eta0 <- predict(fit, X, type = "link", offset = 0)
  expect_equal(eta0, fit$intercept + tt_contraction(fit$cores, basis), tolerance = 1e-8)
  # newdata with training offset recovers linear predictors
  expect_equal(predict(fit, X, type = "link", offset = off),
               fit$linear.predictors, tolerance = 1e-8)
})

test_that("GLAM and TT with same offset recover similar means on grid", {
  dat <- simulate_glam_poisson(n_age = 21L, n_year = 17L, seed = 5L)
  bb <- glam_grid_bases(list(age = dat$age, year = dat$year), k = 8)
  glam <- glam_fit_poisson(dat$Y, bb$B, lambda = c(5, 1),
                           offset = log(dat$exposure))
  X <- as.matrix(expand.grid(age = dat$age, year = dat$year))
  y <- as.numeric(dat$Y)
  off <- log(as.numeric(dat$exposure))
  fit <- ttpspline(
    y, X, family = poisson(), rank = 3, k = 8, lambda = c(5, 1),
    offset = off,
    control = tt_control(pirls_maxit = 20, max_sweeps = 8, backend = "R",
                         compute_edf = FALSE)
  )
  rmse_log <- function(a, b) sqrt(mean((log(as.numeric(a)) - log(as.numeric(b)))^2))
  expect_lt(rmse_log(fitted(fit), glam$mu), 0.15)
})
