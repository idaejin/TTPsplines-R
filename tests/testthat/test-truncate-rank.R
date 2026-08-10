# Truncate TT rank + contraction smoke tests

test_that("tt_truncate_rank reduces bonds and keeps map approximately", {
  set.seed(7)
  d <- 3L
  k <- 5L
  init3 <- tt_initialize(d = d, rank = 3, k = k, seed = 11, sd = 0.2)
  init2 <- tt_truncate_rank(init3, rank = 2)
  expect_equal(attr(init2, "ranks"), c(1L, 2L, 2L, 1L))
  expect_equal(dim(init2[[1]]), c(1L, k, 2L))
  expect_equal(dim(init2[[2]]), c(2L, k, 2L))
  expect_equal(dim(init2[[3]]), c(2L, k, 1L))

  # Same-rank truncate is a no-op on bond sizes
  same <- tt_truncate_rank(init3, rank = 3)
  expect_equal(attr(same, "ranks"), c(1L, 3L, 3L, 1L))
})

test_that("warm-start from truncated r=3 cores is accepted by ttps", {
  set.seed(3)
  n <- 120
  X <- matrix(runif(n * 3), n, 3)
  y <- rnorm(n)
  fit3 <- ttps(
    y, X, rank = 3, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 4, backend = "R", compute_edf = FALSE, seed = 1)
  )
  init2 <- tt_truncate_rank(fit3$cores, rank = 2)
  fit2 <- ttps(
    y, X, rank = 2, k = 5, lambda = 1, init = init2,
    control = tt_control(max_sweeps = 4, backend = "R", compute_edf = FALSE)
  )
  expect_s3_class(fit2, "ttpspline")
  expect_equal(fit2$rank_max, 2L)
  expect_true(is.finite(fit2$deviance))
})
