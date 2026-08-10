test_that("unit-interval covariates get open knots on [0, 1]", {
  set.seed(1)
  x <- runif(50)
  kn <- make_knots(x, k = 6, degree = 3)
  expect_equal(min(kn), 0)
  expect_equal(max(kn), 1)
  # Prediction on seq(0,1) stays interior to the span
  B <- bspline_basis(seq(0, 1, length.out = 21), kn, degree = 3)
  expect_true(all(abs(rowSums(B) - 1) < 1e-10))
})

test_that("non-unit data keep empirical knot range", {
  x <- c(2.1, 3.4, 5.0)
  kn <- make_knots(x, k = 5, degree = 3)
  expect_equal(min(kn), 2.1)
  expect_equal(max(kn), 5.0)
})

test_that("eval_marginal_bases warns outside knot span", {
  kn <- list(structure(make_knots(c(0.2, 0.8), k = 5, degree = 3,
                                  xl = 0.2, xr = 0.8),
                       cyclic = FALSE))
  expect_warning(
    eval_marginal_bases(matrix(c(0, 0.5, 1), 3, 1), kn, degree = 3),
    "outside knot span"
  )
})
