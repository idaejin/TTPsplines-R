test_that("lambda = 'CV' first-sweep search then freeze (Gaussian)", {
  set.seed(21)
  n <- 160
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + 0.4 * X[, 2] + rnorm(n, 0, 0.25)
  ctrl <- tt_control(
    max_sweeps = 6, backend = "R", compute_edf = FALSE,
    cv_folds = 3, cv_ngrid = 7, seed = 3, warn_lambda_boundary = FALSE
  )
  ctrl1 <- tt_control(
    max_sweeps = 1, backend = "R", compute_edf = FALSE,
    cv_folds = 3, cv_ngrid = 7, seed = 3, warn_lambda_boundary = FALSE
  )
  fit1 <- ttps(y, X, rank = 1, k = 5, lambda = "CV",
               optimizer = "ALS", control = ctrl1)
  fit <- ttps(y, X, rank = 1, k = 5, lambda = "CV",
              optimizer = "ALS", control = ctrl)
  expect_equal(fit$lambda_method, "CV")
  expect_equal(fit$lambda, fit1$lambda, tolerance = 1e-10)
  expect_true(all(is.finite(fit$lambda) & fit$lambda > 0))
  expect_equal(nrow(fit$cv$trace), 2L)
  expect_equal(fit$n_criterion_evals, 2L * (3L * 7L + 1L))
  expect_true(is.finite(fit$deviance))
})

test_that("lambda = 'CV' Poisson tunes only first PIRLS sweep", {
  set.seed(22)
  n <- 140
  X <- matrix(runif(n * 2), n, 2)
  eta <- 0.4 * sin(2 * pi * X[, 1]) - 0.2 * X[, 2]
  y <- rpois(n, exp(eta - mean(eta) + log(3)))
  ctrl <- tt_control(
    pirls_maxit = 8, max_sweeps = 4, als_sweeps_per_pirls = 1,
    als_sweeps_adaptive = FALSE, backend = "R", compute_edf = FALSE,
    cv_folds = 3, cv_ngrid = 7, seed = 4, warn_lambda_boundary = FALSE
  )
  fit <- ttps(y, X, family = poisson(), rank = 1, k = 5, lambda = "CV",
              optimizer = "PIRLS-ALS", control = ctrl)
  expect_equal(fit$lambda_method, "CV")
  expect_true(all(is.finite(fit$lambda) & fit$lambda > 0))
  expect_equal(nrow(fit$cv$trace), 2L)
  expect_equal(fit$n_criterion_evals, 2L * (3L * 7L + 1L))
  expect_true(is.finite(fit$deviance))
})

test_that("lambda = 'CV' rejects LBFGS", {
  set.seed(1)
  X <- matrix(runif(40 * 2), 40, 2)
  y <- rnorm(40)
  expect_error(
    ttps(y, X, rank = 1, k = 4, lambda = "CV", optimizer = "LBFGS",
         control = tt_control(max_sweeps = 2, compute_edf = FALSE)),
    "ALS"
  )
})

test_that("lambda = 'CV' cv_sweeps retunes after sweep 1", {
  set.seed(23)
  n <- 120
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.3)
  ctrl <- tt_control(
    max_sweeps = 3, backend = "R", compute_edf = FALSE,
    cv_folds = 3, cv_ngrid = 5, cv_sweeps = 2, seed = 5,
    warn_lambda_boundary = FALSE
  )
  fit <- ttps(y, X, rank = 1, k = 5, lambda = "CV",
              optimizer = "ALS", control = ctrl)
  expect_equal(nrow(fit$cv$trace), 4L)
  expect_equal(fit$n_criterion_evals, 2L * 2L * (3L * 5L + 1L))
  expect_equal(fit$cv$sweeps, 2L)
  expect_equal(fit$cv$rule, "min")
})

test_that("tt_control cv_sweeps Inf means retune every sweep", {
  expect_equal(tt_control(cv_sweeps = Inf)$cv_sweeps, .Machine$integer.max)
  expect_equal(tt_control(cv_sweeps = 1)$cv_sweeps, 1L)
  expect_equal(tt_control()$cv_rule, "min")
})

test_that("CV 1se picks the largest lambda on a flat score", {
  set.seed(1)
  n <- 40
  Xw <- cbind(1, rnorm(n))
  yw <- rnorm(n)
  P <- matrix(0, 2, 2)
  S <- crossprod(Xw)
  b <- as.numeric(crossprod(Xw, yw))
  ws <- list(Xw = Xw, yw = yw, P = P, P0 = NULL, lambda0 = 1, S = S, b = b)
  folds <- rep(seq_len(4), length.out = n)
  grid <- 10^seq(-3, 2, length.out = 6)
  fit_min <- update_lambda_cv(ws, folds = folds, grid = grid, rule = "min")
  fit_1se <- update_lambda_cv(ws, folds = folds, grid = grid, rule = "1se")
  expect_equal(fit_min$lambda, grid[[1L]])
  expect_equal(fit_1se$lambda, max(grid))
  expect_equal(fit_1se$cv_rule, "1se")
})
