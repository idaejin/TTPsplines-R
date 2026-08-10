test_that("R and Rcpp Gaussian ALS agree under shared init (eta/mu/deviance)", {
  skip_if_not(exists("tt_fit_d_cpp", mode = "function"),
              "Rcpp backend not compiled")
  set.seed(21)
  n <- 180
  X <- matrix(runif(n * 3), n, 3)
  y <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + rnorm(n, 0, 0.2)
  bs <- build_marginal_bases(X, k = 6, degree = 3)
  init <- tt_initialize(X, rank = 2, k = 6, seed = 123L, sd = 0.15)

  ctrl <- function(backend) {
    tt_control(
      max_sweeps = 8L, backend = backend, compute_edf = FALSE,
      seed = 123L, tol = 1e-10
    )
  }
  fit_r <- ttps(
    y, X, rank = 2, k = 6, lambda = 1, init = init, knots = bs$knots,
    control = ctrl("R")
  )
  fit_c <- ttps(
    y, X, rank = 2, k = 6, lambda = 1, init = init, knots = bs$knots,
    control = ctrl("Rcpp")
  )
  expect_identical(fit_c$backend, "Rcpp")

  eta_r <- as.numeric(predict(fit_r, type = "link"))
  eta_c <- as.numeric(predict(fit_c, type = "link"))
  mu_r <- as.numeric(predict(fit_r, type = "response"))
  mu_c <- as.numeric(predict(fit_c, type = "response"))

  # Same algorithm + ridge policy ⇒ near machine agreement on gauge-free
  # quantities (do not compare TT cores).
  expect_lt(max(abs(eta_r - eta_c)), 1e-6)
  expect_lt(max(abs(mu_r - mu_c)), 1e-6)
  expect_equal(fit_r$deviance, fit_c$deviance, tolerance = 1e-8)
  obj_r <- tt_objective(fit_r, X, y)$value
  obj_c <- tt_objective(fit_c, X, y)$value
  expect_equal(obj_r, obj_c, tolerance = 1e-8)
})

test_that("uniform TT storage identity matches tt_complexity", {
  # N_TT = 2 k r + (d-2) k r^2 for r_0 = r_d = 1, uniform interior r, p_j = k
  k <- 10L
  r <- 3L
  for (d in c(3L, 5L, 7L, 9L, 10L)) {
    cx <- tt_complexity(d = d, p = k, rank = r)
    exact <- 2 * k * r + (d - 2) * k * r * r
    expect_equal(cx$n_tt_stored, as.numeric(exact))
    expect_equal(cx$n_full, as.numeric(k)^d)
  }
})
