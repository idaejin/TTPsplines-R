test_that("monitor enables trace and prints progress", {
  set.seed(7)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  y <- rnorm(n)
  ctrl <- tt_control(monitor = TRUE, max_sweeps = 3L, backend = "R",
                     compute_edf = FALSE)
  expect_true(isTRUE(ctrl$trace))
  expect_true(isTRUE(ctrl$monitor))

  out <- capture.output({
    fit <- ttpspline(
      y, X, family = gaussian(), rank = 2, k = 4, lambda = 1,
      monitor = TRUE,
      control = tt_control(max_sweeps = 3L, backend = "R", compute_edf = FALSE)
    )
    invisible(fit)
  })
  expect_true(any(grepl("^TTPsplines \\|", out)))
  expect_true(any(grepl("ALS sweep", out)))
  expect_equal(fit$optimizer_used, "ALS")
})

test_that("tt_control monitor aliases trace", {
  expect_false(tt_control()$trace)
  expect_true(tt_control(monitor = TRUE)$trace)
  expect_true(tt_control(trace = TRUE)$monitor)
  expect_true(tt_control(trace = FALSE, monitor = TRUE)$trace)
})
