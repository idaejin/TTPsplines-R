test_that("compare_glam_tt_gaussian returns GLAM and TT rows", {
  tab <- compare_glam_tt_gaussian(
    d = 3, n_grid = c(8, 7, 6), k = 5, ranks = 1:2,
    max_sweeps = 4L, seed = 2L
  )
  expect_true(any(tab$method == "GLAM"))
  expect_true(any(grepl("^TT-r", tab$method)))
  expect_equal(tab$d[1], 3L)
  expect_true(all(is.finite(tab$rmse_truth)))
  expect_true(all(tab$time_s >= 0))
  expect_equal(tab$compression[tab$method == "GLAM"], 1)
})
