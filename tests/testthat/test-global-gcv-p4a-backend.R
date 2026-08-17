# P4A: global-GCV lab wired to P3 fixed-λ backends

test_that("P4A dispatcher: R vs Rcpp_fixed agree on RSS/GDF/GCV", {
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"),
              "P3 C++ fitter not compiled")

  set.seed(101)
  n <- 60L
  d <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2] + rnorm(n, 0, 0.15)
  rank <- 2L
  k <- 6L
  lam <- c(1.0, 1.5)
  probes <- TTPsplines:::.tt_lab_rademacher_probes(n, M = 3L, probe_seed = 7L)
  ctrl <- tt_control(max_sweeps = 12L, tol = 1e-8, compute_edf = FALSE, seed = 101L)
  init <- tt_initialize(d = d, rank = rank, k = k, seed = 101L)

  g_r <- TTPsplines:::tt_global_gcv(
    lambda = lam, y = y, X = X, rank = rank, probes = probes,
    M = 3L, control = ctrl, init = init, fit_backend = "R",
    k = k, on_nonconverged = "na", warm_start = TRUE
  )
  g_c <- TTPsplines:::tt_global_gcv(
    lambda = lam, y = y, X = X, rank = rank, probes = probes,
    M = 3L, control = ctrl, init = TTPsplines:::.tt_clone_cores(init),
    fit_backend = "Rcpp_fixed", k = k, on_nonconverged = "na",
    warm_start = TRUE
  )

  expect_true(isTRUE(g_r$fit$fit_backend == "R") ||
                identical(g_r$fit_backend, "R"))
  expect_identical(g_c$fit_backend, "Rcpp_fixed")
  expect_identical(g_c$fit$fit_backend, "Rcpp_fixed")

  # Approximate parity (ttps may use fused Gram; P3 uses dense) — loose but real
  expect_lt(abs(g_r$rss - g_c$rss) / (1 + abs(g_r$rss)), 5e-2)
  expect_true(is.finite(g_r$gdf) && is.finite(g_c$gdf))
  expect_lt(abs(g_r$gdf - g_c$gdf) / (1 + abs(g_r$gdf)), 0.25)
  expect_true(is.finite(g_r$global_gcv) && is.finite(g_c$global_gcv))
  expect_lt(
    abs(g_r$global_gcv - g_c$global_gcv) / (1 + abs(g_r$global_gcv)),
    0.25
  )
})

test_that("P4A: .tt_lab_fit_fixed does not mutate init cores", {
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"))
  set.seed(102)
  n <- 40L
  d <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- rnorm(n)
  init <- tt_initialize(d = d, rank = 2L, k = 5L, seed = 102L)
  snap <- lapply(init, function(C) array(as.numeric(C), dim(C)))
  ctrl <- tt_control(max_sweeps = 5L, compute_edf = FALSE, seed = 102L)
  invisible(TTPsplines:::.tt_lab_fit_fixed(
    y, X, lambda = c(1, 1), rank = 2L, control = ctrl,
    fit_backend = "Rcpp_fixed", init = init, k = 5L
  ))
  for (j in seq_along(init)) {
    expect_equal(as.numeric(init[[j]]), as.numeric(snap[[j]]))
  }
})

test_that("P4A: GDF refits inherit fit_backend from base", {
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"))
  set.seed(103)
  n <- 45L
  d <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- rnorm(n)
  ctrl <- tt_control(max_sweeps = 8L, compute_edf = FALSE, seed = 103L)
  fit <- TTPsplines:::.tt_lab_fit_fixed(
    y, X, lambda = c(0.8, 1.2), rank = 2L, control = ctrl,
    fit_backend = "Rcpp_fixed", k = 5L
  )
  expect_identical(fit$fit_backend, "Rcpp_fixed")
  fit2 <- TTPsplines:::.tt_lab_refit_from_base(
    y + 0.01, fit, warm_start = TRUE, control = ctrl
  )
  expect_identical(fit2$fit_backend, "Rcpp_fixed")
})
