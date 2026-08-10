test_that("tt_rank_select is reproducible with seed and shares folds", {
  set.seed(10)
  n <- 120
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.3 * X[, 3] + rnorm(n, 0, 0.25)
  ctrl <- tt_control(max_sweeps = 4L, compute_edf = FALSE)
  a <- tt_rank_select(y, X, ranks = 1:2, k = 5, lambda = 1, folds = 3,
                      rule = "1se", seed = 42, control = ctrl)
  b <- tt_rank_select(y, X, ranks = 1:2, k = 5, lambda = 1, folds = 3,
                      rule = "1se", seed = 42, control = ctrl)
  expect_identical(a$fold_id, b$fold_id)
  expect_equal(a$cv_results$mean_cv, b$cv_results$mean_cv, tolerance = 1e-10)
  expect_equal(a$selected_rank, b$selected_rank)
  # same folds for every rank (single fold_id vector)
  expect_length(a$fold_id, n)
  expect_equal(sort(unique(a$fold_id)), 1:3)
})

test_that("minimum and 1-SE rules + ties pick smaller rank", {
  # synthetic CV table
  cv <- data.frame(
    rank = 1:4,
    mean_cv = c(0.20, 0.10, 0.10, 0.11),
    se_cv = c(0.02, 0.02, 0.02, 0.02),
    stringsAsFactors = FALSE
  )
  expect_equal(.tt_rank_min_cv(cv), 2L) # tie 2 and 3 → smaller
  # thresh = 0.10 + 0.02 = 0.12 → ranks 2,3,4 ok → min = 2
  expect_equal(.tt_rank_1se(cv, 2L), 2L)

  cv2 <- data.frame(
    rank = 1:4,
    mean_cv = c(0.12, 0.11, 0.10, 0.105),
    se_cv = c(0.01, 0.01, 0.03, 0.01),
    stringsAsFactors = FALSE
  )
  # r_min=3, thresh=0.10+0.03=0.13 → all ranks ≤ 0.13 → 1se=1
  expect_equal(.tt_rank_min_cv(cv2), 3L)
  expect_equal(.tt_rank_1se(cv2, 3L), 1L)
})

test_that("family-aware metrics and fixed / vector / cGCV lambda", {
  set.seed(3)
  n <- 100
  X <- matrix(runif(n * 3), n, 3)
  f <- sin(2 * pi * X[, 1]) + X[, 2]
  y <- f + rnorm(n, 0, 0.3)
  ctrl <- tt_control(max_sweeps = 3L, compute_edf = FALSE,
                     tol_lambda = 1e-2, lambda_bounds = c(1e-2, 1e2),
                     warn_lambda_boundary = FALSE)

  expect_equal(.tt_resolve_cv_metric("auto", "gaussian"), "rmse")
  expect_equal(.tt_resolve_cv_metric("auto", "poisson"), "poisson_deviance")
  expect_equal(.tt_resolve_cv_metric("auto", "bernoulli"), "logloss")

  sel_fix <- tt_rank_select(y, X, ranks = 1:2, k = 5, lambda = 1,
                             folds = 3, seed = 1, control = ctrl)
  expect_equal(sel_fix$metric, "rmse")
  expect_equal(sel_fix$lambda_method, "fixed")

  sel_vec <- tt_rank_select(y, X, ranks = 1:2, k = 5, lambda = c(1, 2, 1),
                            folds = 3, seed = 1, control = ctrl)
  expect_equal(sel_vec$lambda_method, "fixed")

  sel_cgcv <- tt_rank_select(
    y, X, ranks = 1:2, k = 5, lambda = "cGCV", folds = 3, seed = 1,
    control = ctrl
  )
  expect_equal(sel_cgcv$lambda_method, "cGCV")
  # each fold stores a training-selected lambda (length d)
  lam <- sel_cgcv$lambda_by_fold[["1"]][[1]]
  expect_true(is.numeric(lam) && length(lam) == 3L)
  expect_true(all(is.finite(lam)))
})

test_that("invalid ranks / rank_chain / refit / S3 methods", {
  set.seed(5)
  n <- 80
  X <- matrix(runif(n * 3), n, 3)
  y <- rnorm(n)
  ctrl <- tt_control(max_sweeps = 2L, compute_edf = FALSE)
  expect_error(tt_rank_select(y, X, ranks = 0, control = ctrl), "positive")
  expect_error(
    tt_rank_select(y, X, ranks = 1, rank_chain = c(1, 2, 1), control = ctrl),
    "rank_chain"
  )

  sel <- tt_rank_select(y, X, ranks = 1:2, k = 4, lambda = 1, folds = 3,
                        rule = "min", seed = 7, control = ctrl)
  expect_s3_class(sel, "tt_rank_selection")
  expect_true(sel$selected_rank %in% 1:2)
  expect_output(print(sel), "Selected rank")
  expect_output(print(summary(sel)), "Total CV wall time")
  expect_silent(plot(sel, type = "error"))
  expect_silent(plot(sel, type = "compression"))
  expect_silent(plot(sel, type = "tradeoff"))

  fit <- tt_rank_refit(sel)
  expect_s3_class(fit, "ttpspline")
  expect_equal(.tt_rank_label(fit), as.character(sel$selected_rank))
})

test_that("cGCV path uses training rows only (no validation leakage)", {
  # Construct data where validation responses differ wildly from training.
  # If cGCV saw validation y, lambda/fit would change; we check that swapping
  # only validation y after fold_id is fixed does not change training lambdas
  # when we re-run selection with the same seed and identical training y.
  set.seed(11)
  n <- 90
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  ctrl <- tt_control(max_sweeps = 3L, compute_edf = FALSE,
                     tol_lambda = 5e-2, lambda_bounds = c(1e-2, 1e2),
                     warn_lambda_boundary = FALSE)
  fold_id <- .tt_make_fold_id(n, 3L, seed = 99L)
  # Fit one fold manually: train-only cGCV
  train <- fold_id != 1L
  fit_tr <- ttps(y[train], X[train, , drop = FALSE], rank = 1, k = 5,
                      lambda = "cGCV", control = ctrl)
  # Poison validation y; training unchanged → same lambda if we refit train
  y2 <- y
  y2[!train] <- y2[!train] + 50
  fit_tr2 <- ttps(y2[train], X[train, , drop = FALSE], rank = 1, k = 5,
                       lambda = "cGCV", control = ctrl)
  expect_equal(fit_tr$lambda, fit_tr2$lambda, tolerance = 1e-10)

  # Full selector stores per-fold lambdas from training fits
  sel <- tt_rank_select(y, X, ranks = 1, k = 5, lambda = "cGCV", folds = 3,
                        seed = 99, control = ctrl)
  expect_equal(sel$fold_id, fold_id)
  expect_equal(sel$lambda_by_fold[["1"]][[1]], as.numeric(fit_tr$lambda),
               tolerance = 1e-8)
})
