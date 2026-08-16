test_that("ttps stores joint block-margin EDF that sums; cond is separate", {
  set.seed(31)
  n <- 80
  X <- matrix(runif(n * 3), n, 3)
  colnames(X) <- c("x1", "x2", "x3")
  y <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + stats::rnorm(n, sd = 0.25)
  fit <- ttps(
    y, X, rank = 2, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 5, compute_edf = TRUE, backend = "R")
  )
  expect_true(is.finite(fit$edf))
  expect_equal(length(fit$edf_margin), 3L)
  expect_true(all(is.finite(fit$edf_margin)))
  expect_equal(names(fit$edf_margin), colnames(X))
  expect_equal(sum(fit$edf_margin), fit$edf, tolerance = 1e-6)
  expect_equal(length(fit$edf_margin_cond), 3L)
  expect_true(all(is.finite(fit$edf_margin_cond)))

  ed <- tt_edf(fit)
  expect_s3_class(ed, "tt_edf")
  expect_equal(ed$joint, fit$edf)
  expect_equal(ed$margin, fit$edf_margin)
  expect_equal(ed$margin_cond, fit$edf_margin_cond)
  expect_true(isTRUE(ed$additive))
})

test_that("tt_edf recomputes when compute_edf was FALSE", {
  set.seed(32)
  n <- 60
  X <- matrix(runif(n * 3), n, 3)
  colnames(X) <- paste0("x", 1:3)
  y <- rowSums(sin(2 * pi * X[, 1:2])) + stats::rnorm(n, sd = 0.3)
  fit <- ttps(
    y, X, rank = 2, k = 4, lambda = 1,
    control = tt_control(max_sweeps = 4, compute_edf = FALSE, backend = "R")
  )
  expect_true(is.na(fit$edf))
  expect_null(fit$edf_margin)
  ed <- tt_edf(fit)
  expect_true(is.finite(ed$joint))
  expect_equal(length(ed$margin), 3L)
  expect_true(all(is.finite(ed$margin)))
  expect_equal(sum(ed$margin), ed$joint, tolerance = 1e-6)
  expect_equal(length(ed$margin_cond), 3L)
})
