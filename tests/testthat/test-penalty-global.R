# Classical multidimensional P-spline penalty on Θ (TT representation only).

test_that("own_margin / separable modes are rejected", {
  expect_identical(normalize_penalty_mode("global"), "global")
  expect_identical(normalize_penalty_mode("full"), "global")
  expect_error(normalize_penalty_mode("own_margin"), "removed")
  expect_error(normalize_penalty_mode("separable"), "removed")
  expect_identical(tt_control()$penalty_mode, "global")
})

test_that("P_k^full matches A_k^T S_lambda A_k on small d (every core)", {
  set.seed(11)
  for (r in c(1L, 2L)) {
    for (p in c(4L, 5L)) {
      ranks <- c(1L, rep(r, 2L), 1L)
      d <- 3L
      cores <- initialize_tt_cores(p, ranks, seed = 11L + r + p, sd = 0.25)
      lambda <- c(0.7, 1.1, 0.4)
      p_vec <- rep(p, d)
      S <- glam_penalty(p_vec, lambda, penalty_order = 2L)
      for (k in seq_len(d)) {
        mk <- length(cores[[k]])
        A <- matrix(0, prod(p_vec), mk)
        for (i in seq_len(mk)) {
          g <- numeric(mk)
          g[i] <- 1
          cc <- cores
          cc[[k]] <- array(g, dim(cores[[k]]))
          A[, i] <- as.numeric(tt_full_theta(cc))
        }
        P_ref <- crossprod(A, S %*% A)
        P_ref <- 0.5 * (P_ref + t(P_ref))
        pen <- tt_core_penalty_full(
          cores, k, lambda, penalty_order = 2L, max_dense = 0L
        )
        expect_lt(max(abs(pen$P_full - P_ref)), 1e-8)
      }
    }
  }
})

test_that("rank-1 d=2 closed form for P_a and J", {
  set.seed(2)
  p <- 5L
  a <- rnorm(p)
  b <- rnorm(p)
  cores <- list(
    array(a, c(1L, p, 1L)),
    array(b, c(1L, p, 1L))
  )
  D <- diff(diag(p), differences = 2L)
  DtD <- crossprod(D)
  lam <- c(1.3, 0.6)
  J_analytic <- lam[1] * sum((D %*% a)^2) * sum(b^2) +
    lam[2] * sum(a^2) * sum((D %*% b)^2)
  expect_equal(
    tt_global_penalty_value(cores, lam, penalty_order = 2L),
    0.5 * J_analytic,
    tolerance = 1e-10
  )
  P_a <- lam[1] * sum(b^2) * DtD + lam[2] * sum((D %*% b)^2) * diag(p)
  pen_a <- tt_core_penalty_full(cores, 1L, lam, penalty_order = 2L, max_dense = 0L)
  expect_lt(max(abs(pen_a$P_full - P_a)), 1e-9)
  P_b <- lam[2] * sum(a^2) * DtD + lam[1] * sum((D %*% a)^2) * diag(p)
  pen_b <- tt_core_penalty_full(cores, 2L, lam, penalty_order = 2L, max_dense = 0L)
  expect_lt(max(abs(pen_b$P_full - P_b)), 1e-9)
})

test_that("rank-1 d=3 closed form for J", {
  set.seed(3)
  p <- 4L
  a <- rnorm(p)
  b <- rnorm(p)
  c <- rnorm(p)
  cores <- list(
    array(a, c(1L, p, 1L)),
    array(b, c(1L, p, 1L)),
    array(c, c(1L, p, 1L))
  )
  D <- diff(diag(p), differences = 2L)
  lam <- c(0.5, 1.2, 0.8)
  J_analytic <-
    lam[1] * sum((D %*% a)^2) * sum(b^2) * sum(c^2) +
    lam[2] * sum(a^2) * sum((D %*% b)^2) * sum(c^2) +
    lam[3] * sum(a^2) * sum(b^2) * sum((D %*% c)^2)
  expect_equal(
    tt_global_penalty_value(cores, lam, penalty_order = 2L),
    0.5 * J_analytic,
    tolerance = 1e-10
  )
})

test_that("global J is gauge-invariant", {
  set.seed(8)
  ranks <- c(1L, 2L, 2L, 1L)
  p <- 5L
  cores <- initialize_tt_cores(p, ranks, seed = 8L, sd = 0.3)
  lam <- c(0.4, 1.0, 0.7)
  J0 <- tt_global_penalty_value(cores, lam, penalty_order = 2L)
  A <- matrix(c(1.2, 0.3, -0.1, 0.9), 2, 2)
  cores2 <- tt_apply_gauge(cores, iface = 1L, A = A)
  B <- matrix(c(0.8, -0.2, 0.4, 1.1), 2, 2)
  cores3 <- tt_apply_gauge(cores2, iface = 2L, A = B)
  expect_equal(
    tt_global_penalty_value(cores3, lam, penalty_order = 2L),
    J0,
    tolerance = 1e-10
  )
  expect_equal(tt_full_theta(cores3), tt_full_theta(cores), tolerance = 1e-10)
})

test_that("fixed-λ ALS never increases Q after a core update (P6)", {
  set.seed(13)
  n <- 80
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.3 * X[, 3] +
    rnorm(n, 0, 0.15)
  fit <- ttps(
    y, X, rank = 2, k = 6, lambda = c(0.5, 1.2, 0.8),
    control = tt_control(
      max_sweeps = 6L, backend = "R", compute_edf = FALSE,
      seed = 13L, objective_tol = 1e-9
    )
  )
  expect_identical(fit$penalty_mode, "global")
  expect_true(isTRUE(fit$q_descent$checked))
  expect_equal(fit$q_descent$violations, 0L)
  expect_lte(fit$q_descent$max_increase, 1e-8)
  objs <- vapply(fit$history, `[[`, numeric(1), "objective")
  expect_true(all(diff(objs) <= 1e-8 * pmax(1, abs(objs[-length(objs)]))))
})

test_that("ALS stores global penalty_mode", {
  set.seed(3)
  n <- 80
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.15)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = 1,
    control = tt_control(
      max_sweeps = 4L, backend = "R", compute_edf = FALSE, seed = 3L
    )
  )
  expect_identical(fit$penalty_mode, "global")
  expect_true(is.finite(fit$deviance))
  obj <- tt_objective(fit, X, y)
  expect_true(is.finite(obj$value) && obj$penalty >= 0)
})

test_that("Rcpp P_k^full matches R TT path when exported", {
  skip_if_not(exists("tt_conditional_penalty_full_cpp", mode = "function"))
  set.seed(42)
  ranks <- c(1L, 2L, 2L, 1L)
  p <- 4L
  cores <- initialize_tt_cores(p, ranks, seed = 42, sd = 0.3)
  lambda <- c(0.5, 1.0, 0.75)
  k <- 2L
  DtD_list <- tt_DtD_list(cores, 2L)
  pen_r <- tt_conditional_penalty_full_tt(cores, k, lambda, DtD_list)
  pen_c <- tt_conditional_penalty_full_cpp(cores, as.integer(k), lambda, DtD_list)
  expect_equal(pen_c$P_full, pen_r$P_full, tolerance = 1e-8)
  expect_equal(
    tt_global_penalty_value_cpp(cores, lambda, DtD_list),
    tt_global_penalty_value(cores, lambda, penalty_order = 2L),
    tolerance = 1e-10
  )
})

test_that("LBFGS and GD run with finite global objective", {
  set.seed(5)
  n <- 60
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + 0.5 * cos(2 * pi * X[, 2]) + rnorm(n, 0, 0.2)
  ctrl <- tt_control(
    max_sweeps = 2L, backend = "R", compute_edf = FALSE, seed = 5L,
    lbfgs_maxit = 40L, gd_maxit = 80L, outer_maxit = 1L
  )
  fit_l <- ttps(y, X, rank = 2, k = 5, lambda = 1, optimizer = "LBFGS",
                control = ctrl)
  fit_g <- ttps(y, X, rank = 2, k = 5, lambda = 1, optimizer = "GD",
                control = ctrl)
  expect_identical(fit_l$penalty_mode, "global")
  expect_identical(fit_g$penalty_mode, "global")
  expect_true(is.finite(fit_l$deviance) && is.finite(fit_g$deviance))
})

test_that("cGCV uses P_other + lambda_k P_own", {
  set.seed(17)
  n <- 70
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = "cGCV",
    control = tt_control(
      max_sweeps = 5L, backend = "R", compute_edf = FALSE, seed = 17L
    )
  )
  expect_identical(fit$penalty_mode, "global")
  expect_true(all(is.finite(fit$lambda)))
  expect_length(fit$lambda, 2L)
})
