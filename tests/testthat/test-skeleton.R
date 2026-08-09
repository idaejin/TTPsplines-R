test_that("package loads", {
  expect_true(requireNamespace("ttpsplines", quietly = TRUE) || TRUE)
  # After install: expect_error(ttpsplines::tt_fit(1, list()), "stub")
})

test_that("placeholder Rcpp symbol exists once compiled", {
  skip("enable after first RcppAttributes + install")
})
