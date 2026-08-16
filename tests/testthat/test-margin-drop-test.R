test_that("tt_margin_drop_test nested keeps signal margin x1", {
  set.seed(41)
  n <- 120
  d <- 5
  X <- matrix(runif(n * d), n, d)
  colnames(X) <- paste0("x", seq_len(d))
  y <- sin(2 * pi * X[, 1]) + sin(2 * pi * X[, 2]) + stats::rnorm(n, sd = 0.25)
  ctrl <- tt_control(max_sweeps = 4L, compute_edf = TRUE, backend = "R")
  tst <- tt_margin_drop_test(
    y, X, rank = 2L, k = 5L, lambda = 1, method = "nested",
    alpha = 0.05, control = ctrl, verbose = FALSE
  )
  expect_s3_class(tst, "tt_margin_drop_test")
  expect_equal(nrow(tst$results), d)
  expect_true(all(is.finite(tst$results$p_value)))
  # Strong signal x1 should not be a drop candidate
  expect_false(isTRUE(tst$results$drop_candidate[tst$results$margin == "x1"]))
  expect_true("x1" %in% tst$keep_names)
})
