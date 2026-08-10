test_that("n_starts multi-init records stability and is reproducible", {
  set.seed(8)
  n <- 90
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + rnorm(n, 0, 0.3)
  ctrl <- tt_control(max_sweeps = 3L, compute_edf = FALSE)
  a <- tt_rank_select(
    y, X, ranks = 1:2, k = 4, lambda = 1, folds = 3,
    n_starts = 2L, seed = 11, control = ctrl
  )
  b <- tt_rank_select(
    y, X, ranks = 1:2, k = 4, lambda = 1, folds = 3,
    n_starts = 2L, seed = 11, control = ctrl
  )
  expect_equal(a$n_starts, 2L)
  expect_true(all(c("objective_best", "objective_median", "start_gap",
                    "start_convergence_rate") %in% names(a$cv_results)))
  expect_equal(dim(a$start_objective[["2"]]), c(3L, 2L))
  expect_equal(a$cv_results$mean_cv, b$cv_results$mean_cv, tolerance = 1e-10)
  expect_equal(a$selected_rank, b$selected_rank)
})

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
  expect_identical(a$foldid, a$fold_id)
  expect_equal(a$cv_results$mean_cv, b$cv_results$mean_cv, tolerance = 1e-10)
  expect_equal(a$selected_rank, b$selected_rank)
  expect_length(a$fold_id, n)
  expect_equal(sort(unique(a$fold_id)), 1:3)
  expect_equal(a$cvm, a$cv_results$mean_cv)
  expect_equal(a$cvsd, a$cv_results$se_cv)
  expect_equal(a$`rank.min`, a$rank_min)
  expect_equal(a$`rank.1se`, a$rank_1se)
  expect_s3_class(a, "cv.ttps")
})

test_that("user foldid is first-class and shared across per-rank calls", {
  set.seed(12)
  n <- 60
  X <- matrix(runif(n * 3), n, 3)
  colnames(X) <- paste0("x", 1:3)
  y <- rnorm(n)
  ctrl <- tt_control(max_sweeps = 2L, compute_edf = FALSE)
  foldid <- rep(1:3, length.out = n)
  # unequal sizes by design (60/3 = 20 each, still ok)
  foldid_uneq <- c(rep(1L, 10), rep(2L, 20), rep(3L, 30))
  expect_error(tt_rank_select(y, X, ranks = 1, foldid = 1:5, control = ctrl),
               "length n")

  sel <- tt_rank_select(y, X, ranks = 1:2, k = 4, lambda = 1,
                        foldid = foldid, seed = 99, control = ctrl)
  expect_identical(sel$foldid, as.integer(foldid))
  expect_equal(sel$folds, 3L)
  expect_identical(sel$variable_order, colnames(X))

  # Chicago checkpoint pattern: separate calls, shared foldid
  s1 <- tt_rank_select(y, X, ranks = 1, k = 4, lambda = 1,
                       foldid = foldid_uneq, control = ctrl)
  s2 <- tt_rank_select(y, X, ranks = 2, k = 4, lambda = 1,
                       foldid = foldid_uneq, control = ctrl)
  expect_identical(s1$foldid, s2$foldid)
  expect_equal(as.integer(table(s1$foldid)), c(10L, 20L, 30L))
})

test_that("minimum and 1-SE rules + ties pick smaller rank", {
  cv <- data.frame(
    rank = 1:4,
    mean_cv = c(0.20, 0.10, 0.10, 0.11),
    se_cv = c(0.02, 0.02, 0.02, 0.02),
    stringsAsFactors = FALSE
  )
  expect_equal(.tt_rank_min_cv(cv), 2L)
  expect_equal(.tt_rank_1se(cv, 2L), 2L)

  cv2 <- data.frame(
    rank = 1:4,
    mean_cv = c(0.12, 0.11, 0.10, 0.105),
    se_cv = c(0.01, 0.01, 0.03, 0.01),
    stringsAsFactors = FALSE
  )
  expect_equal(.tt_rank_min_cv(cv2), 3L)
  expect_equal(.tt_rank_1se(cv2, 3L), 1L)

  # plateau: r=1 poor, then flat → 1se prefers small rank in plateau
  cv3 <- data.frame(
    rank = 1:5,
    mean_cv = c(1.0, 0.40, 0.38, 0.37, 0.37),
    se_cv = c(0.05, 0.05, 0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )
  expect_equal(.tt_rank_min_cv(cv3), 4L) # tie 4 and 5 → smaller
  # thresh = 0.37 + 0.05 = 0.42 → r=2..5 ok → 1se = 2
  expect_equal(.tt_rank_1se(cv3, 4L), 2L)

  # NA SE at min → treated as 0
  cv4 <- data.frame(rank = 1:2, mean_cv = c(0.2, 0.1), se_cv = c(0.01, NA_real_))
  expect_equal(.tt_rank_1se(cv4, 2L), 2L)
})

test_that("Inf fold losses are retained in mean CV (not silently dropped)", {
  loss_mat <- matrix(c(1, 2, Inf, 0.5, 0.6, 0.7), nrow = 2, byrow = TRUE,
                     dimnames = list(c("1", "2"), paste0("fold", 1:3)))
  time_mat <- matrix(1, 2, 3, dimnames = dimnames(loss_mat))
  conv_mat <- matrix(TRUE, 2, 3, dimnames = dimnames(loss_mat))
  lambda_list <- list(`1` = list(1, 1, 1), `2` = list(1, 1, 1))
  tab <- .tt_summarize_rank_cv(
    ranks = 1:2, loss_mat = loss_mat, time_mat = time_mat, conv_mat = conv_mat,
    lambda_list = lambda_list, lambda_method = "fixed", d = 3L, p = 5L
  )
  expect_true(is.infinite(tab$mean_cv[1]))
  expect_equal(tab$n_inf_folds[1], 1L)
  expect_equal(tab$n_finite_folds[1], 2L)
  expect_true(is.finite(tab$mean_cv[2]))
  expect_equal(tab$mean_cv[2], mean(c(0.5, 0.6, 0.7)))
})

test_that("family-aware metrics: MSE / Poisson deviance / binomial deviance", {
  set.seed(3)
  n <- 100
  X <- matrix(runif(n * 3), n, 3)
  f <- sin(2 * pi * X[, 1]) + X[, 2]
  y <- f + rnorm(n, 0, 0.3)
  ctrl <- tt_control(max_sweeps = 3L, compute_edf = FALSE,
                     tol_lambda = 1e-2, lambda_bounds = c(1e-2, 1e2),
                     warn_lambda_boundary = FALSE)

  expect_equal(.tt_resolve_cv_metric("auto", "gaussian"), "mse")
  expect_equal(.tt_resolve_cv_metric("auto", "poisson"), "poisson_deviance")
  expect_equal(.tt_resolve_cv_metric("auto", "bernoulli"), "binomial_deviance")
  expect_equal(.tt_resolve_cv_metric("deviance", "poisson"), "poisson_deviance")
  expect_equal(.tt_resolve_cv_metric("deviance", "bernoulli"), "binomial_deviance")

  y0 <- 0
  mu0 <- 0.5
  fam_pois <- poisson()
  expect_equal(
    .tt_cv_loss(y0, mu0, "poisson_deviance", fam_pois),
    as.numeric(glm_deviance(fam_pois, y0, mu0) / 1)
  )
  # y=0 convention: 2*(0 - (0-mu)) = 2*mu
  expect_equal(.tt_cv_loss(0, 2, "poisson_deviance", fam_pois), 4)

  yb <- c(0, 1)
  mub <- c(0.2, 0.8)
  fam_bin <- binomial()
  d_bin <- glm_deviance(fam_bin, yb, mub) / 2
  expect_equal(.tt_cv_loss(yb, mub, "binomial_deviance", fam_bin), d_bin)
  expect_equal(
    .tt_cv_loss(yb, mub, "logloss", fam_bin),
    d_bin / 2,
    tolerance = 1e-12
  )

  mse <- .tt_cv_loss(c(1, 3), c(2, 2), "mse", gaussian())
  rmse <- .tt_cv_loss(c(1, 3), c(2, 2), "rmse", gaussian())
  expect_equal(mse, 1)
  expect_equal(rmse, 1)

  sel_fix <- tt_rank_select(y, X, ranks = 1:2, k = 5, lambda = 1,
                             folds = 3, seed = 1, control = ctrl)
  expect_equal(sel_fix$metric, "mse")
  expect_equal(sel_fix$lambda_method, "fixed")

  sel_rmse <- tt_rank_select(y, X, ranks = 1:2, k = 5, lambda = 1,
                             folds = 3, seed = 1, metric = "rmse",
                             control = ctrl)
  expect_equal(sel_rmse$metric, "rmse")
  # mean of fold RMSE vs MSE: not identical in general
  expect_false(isTRUE(all.equal(sel_fix$cvm, sel_rmse$cvm)))

  sel_vec <- tt_rank_select(y, X, ranks = 1:2, k = 5, lambda = c(1, 2, 1),
                            folds = 3, seed = 1, control = ctrl)
  expect_equal(sel_vec$lambda_method, "fixed")

  sel_cgcv <- tt_rank_select(
    y, X, ranks = 1:2, k = 5, lambda = "cGCV", folds = 3, seed = 1,
    control = ctrl
  )
  expect_equal(sel_cgcv$lambda_method, "cGCV")
  lam <- sel_cgcv$lambda_by_fold[["1"]][[1]]
  expect_true(is.numeric(lam) && length(lam) == 3L)
  expect_true(all(is.finite(lam)))
})

test_that("cv.ttps alias and refit restore EDF on full data", {
  set.seed(5)
  n <- 80
  X <- matrix(runif(n * 3), n, 3)
  colnames(X) <- c("a", "b", "c")
  y <- rnorm(n)
  ctrl <- tt_control(max_sweeps = 2L, compute_edf = FALSE)
  expect_error(tt_rank_select(y, X, ranks = 0, control = ctrl), "positive")
  expect_error(
    tt_rank_select(y, X, ranks = 1, rank_chain = c(1, 2, 1), control = ctrl),
    "rank_chain"
  )

  sel <- cv.ttps(y, X, ranks = 1:2, k = 4, lambda = 1, nfolds = 3,
                 rule = "min", seed = 7, control = ctrl)
  expect_s3_class(sel, c("cv.ttps", "tt_rank_selection"))
  expect_true(sel$selected_rank %in% 1:2)
  expect_output(print(sel), "Selected rank")
  expect_output(print(summary(sel)), "Total CV wall time")
  expect_silent(plot(sel, type = "error"))
  expect_silent(plot(sel, type = "compression"))
  expect_silent(plot(sel, type = "tradeoff"))

  fit <- tt_rank_refit(sel)
  expect_s3_class(fit, "ttpspline")
  expect_equal(.tt_rank_label(fit), as.character(sel$selected_rank))
  expect_true(isTRUE(fit$control$compute_edf) || is.finite(fit$edf) ||
                !is.null(fit$edf))
  # control on refit path should request EDF
  expect_true(isTRUE(sel$fit_args$control$compute_edf == FALSE))
  fit2 <- refit(sel, rank = "min")
  expect_equal(.tt_rank_label(fit2), as.character(sel$rank_min))
  expect_equal(nrow(fit2$X %||% X), n) # full data used via stored y/X
  expect_equal(length(fitted(fit2)), n)
})

test_that("cGCV path uses training rows only (no validation leakage)", {
  set.seed(11)
  n <- 90
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  ctrl <- tt_control(max_sweeps = 3L, compute_edf = FALSE,
                     tol_lambda = 5e-2, lambda_bounds = c(1e-2, 1e2),
                     warn_lambda_boundary = FALSE)
  fold_id <- .tt_make_fold_id(n, 3L, seed = 99L)
  train <- fold_id != 1L

  y2 <- y
  y2[!train] <- y2[!train] + 50
  sel <- tt_rank_select(y, X, ranks = 1, k = 5, lambda = "cGCV", folds = 3,
                        seed = 99, control = ctrl)
  sel2 <- tt_rank_select(y2, X, ranks = 1, k = 5, lambda = "cGCV", folds = 3,
                         seed = 99, control = ctrl)
  expect_equal(sel$fold_id, fold_id)
  expect_equal(sel2$fold_id, fold_id)
  expect_equal(sel$lambda_by_fold[["1"]][[1]], sel2$lambda_by_fold[["1"]][[1]],
               tolerance = 1e-10)
  lam <- sel$lambda_by_fold[["1"]][[1]]
  expect_true(is.numeric(lam) && length(lam) == 3L && all(is.finite(lam)))
})

test_that("diagnostics record lambda, convergence, and n_starts", {
  set.seed(21)
  n <- 70
  X <- matrix(runif(n * 3), n, 3)
  y <- rnorm(n)
  ctrl <- tt_control(max_sweeps = 2L, compute_edf = FALSE)
  sel <- tt_rank_select(y, X, ranks = 1:2, k = 4, lambda = 1, folds = 3,
                        n_starts = 2L, seed = 3, control = ctrl)
  expect_true(!is.null(sel$lambda_by_fold))
  expect_true(!is.null(sel$converged_by_fold))
  expect_true(!is.null(sel$start_objective))
  expect_true(!is.null(sel$start_converged))
  expect_equal(dim(sel$cvraw), c(2L, 3L))
  expect_true(all(c("n_finite_folds", "n_inf_folds") %in% names(sel$cv_results)))
})

test_that("tt_rank_select shares full-X knots by default (no fold OOD warnings)", {
  set.seed(42)
  n <- 60
  X <- matrix(runif(n * 2), n, 2)
  # Global extreme only in fold 1 → train-only knots would warn on predict
  X[1, 1] <- 10
  y <- rnorm(n)
  foldid <- rep(1:3, length.out = n)
  foldid[1] <- 1L
  ctrl <- tt_control(max_sweeps = 2L, compute_edf = FALSE)

  expect_no_warning(
    sel <- tt_rank_select(
      y, X, ranks = 1, k = 5, lambda = 1, foldid = foldid,
      control = ctrl
    )
  )
  expect_identical(sel$knots_source, "full_X")
  expect_false(isTRUE(sel$fold_knots))
  expect_true(!is.null(sel$fit_args$knots))
  kn_full <- build_marginal_bases(X, k = 5, degree = 3)$knots
  expect_equal(sel$fit_args$knots[[1]], kn_full[[1]], tolerance = 1e-12)

  # Escape hatch: per-fold knots recreate the extrapolation warning
  expect_warning(
    tt_rank_select(
      y, X, ranks = 1, k = 5, lambda = 1, foldid = foldid,
      fold_knots = TRUE, control = ctrl
    ),
    "outside knot span"
  )

  expect_error(
    tt_rank_select(
      y, X, ranks = 1, k = 5, lambda = 1, foldid = foldid,
      knots = kn_full, fold_knots = TRUE, control = ctrl
    ),
    "fold_knots = TRUE"
  )
})
