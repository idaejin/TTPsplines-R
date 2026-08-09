test_that("compare separates fit_time and rank_selection_time", {
  cmp <- compare_glam_tt_gaussian(
    d = 3, n_grid = c(8, 7, 6), k = 5, ranks = 1:3,
    lambda = 1, folds = 3, max_sweeps = 3, seed = 2,
    detail_ranks = TRUE, n_test = 200
  )
  expect_s3_class(cmp, "glam_tt_compare")
  expect_true(all(c("fit_time", "rank_selection_time",
                    "total_procedure_time", "oracle_probe_time") %in% names(cmp)))

  glam <- cmp[cmp$method == "GLAM", ][1, ]
  expect_equal(glam$rank_selection_time, 0)
  expect_equal(glam$total_procedure_time, glam$fit_time)

  se1 <- cmp[cmp$method == "TT-1SE", ][1, ]
  expect_true(se1$rank_selection_time > 0)
  expect_equal(se1$total_procedure_time,
               se1$fit_time + se1$rank_selection_time, tolerance = 1e-10)

  ora <- cmp[cmp$method == "TT-oracle", ][1, ]
  expect_equal(ora$rank_selection_time, 0)
  # oracle probe cost recorded, but not added into practical total
  expect_true(is.finite(ora$oracle_probe_time))
  expect_equal(ora$total_procedure_time, ora$fit_time, tolerance = 1e-10)

  sel <- attr(cmp, "selection")
  expect_equal(se1$rank, sel$rank_1se)
  expect_equal(cmp$rank[cmp$method == "TT-minCV"][1], sel$rank_min)

  # same test set used for all RMSE
  expect_true(!is.null(attr(cmp, "test")))
  expect_true(all(is.finite(cmp$rmse_truth[cmp$status == "ok"])))

  expect_silent(plot(cmp, type = "rmse_rank"))
  expect_silent(plot(cmp, type = "tradeoff"))
})

test_that("GLAM infeasible is not_run with NA metrics", {
  sc <- compare_glam_tt_gaussian(
    d = 5, n_grid = 5, k = 6, ranks = 1:2, folds = 3,
    max_sweeps = 2, max_glam_npar = 5000, seed = 4, n_test = 150
  )
  glam <- sc[grepl("^GLAM", sc$method), ][1, ]
  expect_equal(glam$status, "not_run")
  expect_true(is.na(glam$rmse_truth))
  expect_true(is.na(glam$fit_time))
  expect_true(any(sc$method == "TT-1SE"))
})

test_that("minCV/1SE ranks match tt_rank_select; no truth in CV", {
  set.seed(9)
  cmp <- compare_glam_tt_gaussian(
    d = 3, n_grid = 7, k = 5, ranks = 1:3, folds = 3,
    max_sweeps = 3, seed = 9, n_test = 100
  )
  sel <- attr(cmp, "selection")
  expect_equal(cmp$rank[cmp$method == "TT-minCV"], sel$rank_min)
  expect_equal(cmp$rank[cmp$method == "TT-1SE"], sel$rank_1se)
  # selection object has no truth field
  expect_null(sel$truth)
  expect_null(sel$f_true)
})
