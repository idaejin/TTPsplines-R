test_that("normalize_linear and tt_update_intercept_beta work", {
  n <- 50
  L <- cbind(z1 = rnorm(n), z2 = runif(n))
  expect_null(normalize_linear(NULL, n))
  Ln <- normalize_linear(L, n)
  expect_equal(dim(Ln), c(n, 2L))
  expect_warning(normalize_linear(cbind(1, L), n), "intercept")

  y <- 2 + 3 * L[, 1] - L[, 2] + rnorm(n, 0, 0.05)
  ab <- tt_update_intercept_beta(y, offset = 0, f = 0, linear = Ln, weights = 1)
  expect_equal(ab$intercept, 2, tolerance = 0.2)
  expect_equal(unname(ab$beta[1]), 3, tolerance = 0.2)
})

test_that("Gaussian ALS recovers linear coef with TT noise surface", {
  set.seed(1)
  n <- 120
  X <- matrix(runif(n * 2), n, 2)
  colnames(X) <- c("x1", "x2")
  z <- rnorm(n)
  f <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2])
  y <- 1.5 + 2 * z + f + rnorm(n, 0, 0.15)
  L <- cbind(z = z)
  fit <- ttps(
    y, X, linear = L, rank = 2, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 8, compute_edf = FALSE)
  )
  expect_true(length(fit$beta) == 1L)
  expect_equal(unname(fit$beta), 2, tolerance = 0.35)
  mu <- predict(fit, X, linear = L)
  expect_equal(length(mu), n)
  expect_true(cor(mu, y) > 0.85)
  expect_error(predict(fit, X), "linear")
  expect_error(
    ttps(y, X, linear = L, optimizer = "LBFGS",
         control = tt_control(compute_edf = FALSE)),
    "linear="
  )
})

test_that("Poisson PIRLS accepts linear= and improves deviance vs TT-only", {
  set.seed(2)
  n <- 100
  X <- matrix(runif(n * 2), n, 2)
  z <- rnorm(n)
  eta <- 0.3 + 0.8 * z + 0.5 * X[, 1]
  y <- rpois(n, exp(eta))
  L <- cbind(z = z)
  ctrl <- tt_control(
    pirls_maxit = 12, outer_maxit = 4, als_sweeps_adaptive = FALSE,
    compute_edf = FALSE, tol = 1e-5
  )
  fit0 <- ttps(y, X, family = poisson(), rank = 2, k = 5, lambda = 1,
               control = ctrl)
  fit1 <- ttps(y, X, family = poisson(), rank = 2, k = 5, lambda = 1,
               linear = L, control = ctrl)
  expect_true(length(fit1$beta) == 1L)
  expect_true(fit1$deviance < fit0$deviance + 1e-6)
  expect_true(sign(fit1$beta) == sign(0.8))
})

test_that("summary() prints glm-style parametric coefficient table", {
  set.seed(4)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  z <- rnorm(n)
  y <- 1 + 1.5 * z + 0.4 * X[, 1] + rnorm(n, 0, 0.2)
  L <- cbind(z = z)
  fit <- ttps(
    y, X, linear = L, rank = 1, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 6, compute_edf = TRUE)
  )
  s <- summary(fit)
  expect_true(is.matrix(s$coefficients))
  expect_equal(colnames(s$coefficients),
               c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  expect_true("(Intercept)" %in% rownames(s$coefficients))
  expect_true("z" %in% rownames(s$coefficients))
  expect_true(is.finite(s$coefficients["z", "Std. Error"]))
  expect_equal(unname(s$coefficients["z", "Estimate"]), unname(fit$beta))
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "Parametric coefficients:")
  expect_match(out, "Std\\. Error")
})

test_that("Poisson summary uses z table; aliased DOW get NA SE", {
  set.seed(5)
  n <- 100
  X <- matrix(runif(n * 2), n, 2)
  dow <- factor(sample(0:6, n, replace = TRUE))
  L <- model.matrix(~ 0 + dow)
  eta <- 0.2 + 0.3 * X[, 1]
  y <- rpois(n, exp(eta))
  fit <- ttps(
    y, X, family = poisson(), linear = L, rank = 1, k = 4, lambda = 1,
    control = tt_control(pirls_maxit = 8, outer_maxit = 1,
                         als_sweeps_adaptive = FALSE, compute_edf = FALSE)
  )
  s <- summary(fit)
  expect_equal(colnames(s$coefficients)[3], "z value")
  # full DOW dummies + intercept => at least one aliased column
  expect_true(any(s$aliased))
  expect_true(anyNA(s$coefficients[, "Std. Error"]))
})

test_that("tt_rank_select subsets linear by fold", {
  set.seed(3)
  n <- 60
  X <- matrix(runif(n * 2), n, 2)
  z <- rnorm(n)
  y <- rnorm(n, 0.5 * z)
  L <- cbind(z = z)
  ctrl <- tt_control(max_sweeps = 3, compute_edf = FALSE)
  sel <- tt_rank_select(
    y, X, ranks = 1:2, k = 4, lambda = 1, folds = 3, seed = 1,
    linear = L, control = ctrl
  )
  expect_equal(nrow(sel$linear), n)
  fit <- tt_rank_refit(sel)
  expect_true(length(fit$beta) == 1L)
})
