test_that("null_space='profiled' rejects cyclic margins", {
  set.seed(78)
  n <- 40
  X <- cbind(runif(n), runif(n))
  y <- rnorm(n)
  expect_error(
    ttps(y, X, rank = 2L, k = 5L, lambda = 1, null_space = "profiled",
         cyclic = c(TRUE, FALSE),
         control = tt_control(max_sweeps = 2L, compute_edf = FALSE, backend = "R")),
    "cyclic"
  )
})

test_that("null_space='profiled' recovers null truth; Q0 identity holds", {
  set.seed(81)
  n <- 120
  X <- cbind(runif(n), runif(n))
  colnames(X) <- c("x1", "x2")
  f <- 1.5 + 2 * X[, 1] - X[, 2]
  y <- f + rnorm(n, sd = 0.05)
  ctrl <- tt_control(max_sweeps = 10L, compute_edf = TRUE, backend = "R", seed = 81)

  fit_p <- ttps(y, X, rank = 2L, k = 6L, lambda = 1e6, penalty_order = 2L,
                null_space = "profiled", control = ctrl)
  expect_identical(fit_p$null_space, "profiled")
  expect_identical(fit_p$null_space_info$method, "profiled")
  expect_true(sqrt(mean((fitted(fit_p) - f)^2)) < 0.08)

  ns <- fit_p$null_space_info
  expect_true(max(abs(ns$ortho_X0_mu_tt_perp)) < 1e-6)
  expect_true(sqrt(mean(ns$mu_tt_perp^2)) < 1e-4)
  expect_true(ns$identity_max_abs < 1e-8)
  p0 <- ns$design_rank
  expect_true(is.finite(fit_p$edf))
  expect_true(fit_p$edf + 1e-8 >= p0)

  pred <- predict(fit_p, X[1:5, ])
  expect_length(pred, 5L)
  expect_true(all(is.finite(pred)))
})

test_that("profiled vs joint: Q0 identity and finite fit (mixed)", {
  set.seed(82)
  n <- 150
  X <- cbind(runif(n), runif(n))
  colnames(X) <- c("x1", "x2")
  f <- 1 + X[, 1] - 0.5 * X[, 2] +
    0.4 * sin(4 * pi * X[, 1]) * sin(4 * pi * X[, 2])
  y <- f + rnorm(n, sd = 0.1)
  ctrl <- tt_control(max_sweeps = 12L, compute_edf = FALSE, backend = "R", seed = 82)
  lam <- 10
  fit_j <- ttps(y, X, rank = 3L, k = 8L, lambda = lam, null_space = "joint",
                control = ctrl)
  fit_p <- ttps(y, X, rank = 3L, k = 8L, lambda = lam, null_space = "profiled",
                control = ctrl)

  ns <- fit_p$null_space_info
  expect_true(max(abs(ns$ortho_X0_mu_tt_perp)) < 1e-6)
  expect_true(ns$identity_max_abs < 1e-8)
  expect_true(sqrt(mean((fitted(fit_p) - f)^2)) < 0.25)
  expect_true(sqrt(mean((fitted(fit_j) - f)^2)) < 0.25)
})
