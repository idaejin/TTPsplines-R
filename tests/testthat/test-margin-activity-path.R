test_that("tt_margin_activity_path scores and 1se select recover signal margins", {
  set.seed(21)
  n <- 100
  d <- 5
  X <- matrix(runif(n * d), n, d)
  colnames(X) <- paste0("x", seq_len(d))
  y <- sin(2 * pi * X[, 1]) + sin(2 * pi * X[, 2]) + stats::rnorm(n, sd = 0.3)
  ctrl <- tt_control(
    max_sweeps = 4L, compute_edf = FALSE, outer_maxit = 4L, backend = "R"
  )
  path <- tt_margin_activity_path(
    y, X,
    rank = 2L, k = 5L,
    lambda_path = 10^c(1, 0, -1),
    select = "1se", folds = 3L, seed = 21,
    control = ctrl, verbose = FALSE
  )
  expect_s3_class(path, "tt_margin_activity_path")
  expect_equal(length(path$scores), d)
  expect_true(all(is.finite(path$scores)))
  expect_true(all(path$selected %in% seq_len(d)))
  # Top scores should include the two signal margins under this DGP
  top2 <- names(sort(path$scores, decreasing = TRUE))[1:2]
  expect_true(all(c("x1", "x2") %in% top2) || all(c("x1", "x2") %in% path$selected_names))
  expect_true(!is.null(path$cv))
  expect_true(path$m_hat >= 1L)
})

test_that("tt_margin_activity_path select=none returns ranking only", {
  set.seed(22)
  n <- 60
  X <- matrix(runif(n * 4), n, 4)
  colnames(X) <- paste0("x", 1:4)
  y <- X[, 1] + stats::rnorm(n, sd = 0.2)
  ctrl <- tt_control(max_sweeps = 3L, compute_edf = FALSE, backend = "R")
  path <- tt_margin_activity_path(
    y, X, rank = 2L, k = 4L,
    lambda_path = c(1, 0.1),
    select = "none", control = ctrl, verbose = FALSE
  )
  expect_null(path$cv)
  expect_true(is.na(path$m_hat))
  expect_equal(length(path$selected), 0L)
  expect_null(path$fit)
})
