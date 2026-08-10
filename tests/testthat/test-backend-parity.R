test_that("Rcpp P_k^full helpers agree with R TT path (gauge-free)", {
  skip_if_not(exists("tt_conditional_penalty_full_cpp", mode = "function"),
              "Rcpp backend not compiled")
  set.seed(21)
  ranks <- c(1L, 2L, 2L, 1L)
  p <- 5L
  cores <- initialize_tt_cores(p, ranks, seed = 21L, sd = 0.2)
  lambda <- c(0.9, 1.1, 0.7)
  DtD_list <- tt_DtD_list(cores, 2L)
  for (k in 1:3) {
    pen_r <- tt_conditional_penalty_full_tt(cores, k, lambda, DtD_list)
    pen_c <- tt_conditional_penalty_full_cpp(cores, as.integer(k), lambda, DtD_list)
    expect_lt(max(abs(pen_r$P_full - pen_c$P_full)), 1e-8)
  }
  expect_equal(
    tt_global_penalty_value_cpp(cores, lambda, DtD_list),
    tt_global_penalty_value(cores, lambda, penalty_order = 2L),
    tolerance = 1e-10
  )
})

test_that("requesting backend=Rcpp for ALS falls back to R global path", {
  skip_if_not(exists("tt_fit_d_cpp", mode = "function"),
              "Rcpp backend not compiled")
  set.seed(22)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  expect_warning(
    fit <- ttps(
      y, X, rank = 2, k = 5, lambda = 1,
      control = tt_control(
        max_sweeps = 3L, backend = "Rcpp", compute_edf = FALSE, seed = 22L
      )
    ),
    "global penalty"
  )
  expect_identical(fit$backend, "R")
  expect_identical(fit$penalty_mode, "global")
  expect_true(is.finite(fit$deviance))
})
