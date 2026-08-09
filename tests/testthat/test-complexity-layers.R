test_that("tt_complexity separates storage, intrinsic dim, and EDF", {
  cx <- tt_complexity(d = 8, p = 10, rank = 3)
  expect_equal(cx$n_full, 1e8)
  expect_equal(cx$n_tt_stored, 600)
  expect_equal(cx$compression_ratio, 1e8 / 600)
  # gauge = (d-1)*r^2 = 7*9 = 63
  expect_equal(cx$n_gauge, 63)
  expect_equal(cx$n_tt_intrinsic, 600 - 63)
  expect_true(is.na(cx$edf))
})

test_that("tt_complexity on a fit attaches EDF without confusing it with N_TT", {
  set.seed(21)
  n <- 120
  X <- matrix(runif(n * 3), n, 3)
  y <- rnorm(n)
  fit <- ttpspline(
    y, X, rank = 2, k = 5, lambda = 2,
    control = tt_control(max_sweeps = 5, backend = "R", compute_edf = TRUE)
  )
  cx <- tt_complexity(fit)
  expect_equal(cx$n_tt_stored, fit$npar_tt)
  expect_equal(cx$n_tt_intrinsic, fit$npar_tt_intrinsic)
  expect_true(is.finite(cx$edf))
  expect_false(isTRUE(all.equal(cx$edf, cx$n_tt_stored)))
})
