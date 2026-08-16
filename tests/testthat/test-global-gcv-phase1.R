# Phase 1 lab tests (experimental). Heavy cases gated by NOT_CRAN.

test_that("Phase1 cache key and flatness helpers", {
  k1 <- TTPsplines:::.tt_lab_cache_key("d1", 3, c(0, 1), "cold_common", "bank0", M = 5L)
  k2 <- TTPsplines:::.tt_lab_cache_key("d1", 3, c(0, 1), "cold_common", "bank0", M = 5L)
  k3 <- TTPsplines:::.tt_lab_cache_key("d1", 3, c(0, 1.0000001), "cold_common", "bank0", M = 5L)
  expect_identical(k1, k2)
  expect_false(identical(k1, k3))
  flat <- TTPsplines:::.tt_lab_flatness_mask(c(1, 1.005, 1.02, NA, Inf), 0.01)
  expect_equal(flat, c(TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("Phase1 RNG restoration around design + probes", {
  set.seed(123)
  junk <- runif(3)
  before <- .Random.seed
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    "smooth_smooth", n = 40L, n_test = 20L, k = 5L, seed = 7L
  )
  probes <- TTPsplines:::.tt_lab_rademacher_probes(des$n, 4L, probe_seed = 9L)
  after <- .Random.seed
  expect_identical(before, after)
  expect_equal(nrow(probes), des$n)
})

test_that("Phase1 cache hits and core non-mutation", {
  skip_on_cran()
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    "smooth_smooth", n = 50L, n_test = 30L, k = 5L, seed = 11L
  )
  probes <- TTPsplines:::.tt_lab_rademacher_probes(des$n, 6L, probe_seed = 3L)
  init <- TTPsplines:::.tt_with_preserved_seed({
    tt_initialize(d = 2L, rank = 3L, k = des$k, seed = 11L)
  })
  init_copy <- TTPsplines:::.tt_clone_cores(init)
  cache <- TTPsplines:::.tt_lab_new_cache()
  ctrl <- TTPsplines:::.tt_lab_phase1_control(seed = 11L, max_sweeps = 15L)
  e1 <- TTPsplines:::.tt_lab_eval_ggcv(
    lambda = c(1, 1), design = des, rank = 3L, probes = probes, M = 2L,
    init_cores = init, init_policy = "cold_common", mc_bank_id = "b",
    cache = cache, control = ctrl
  )
  e2 <- TTPsplines:::.tt_lab_eval_ggcv(
    lambda = c(1, 1), design = des, rank = 3L, probes = probes, M = 2L,
    init_cores = init, init_policy = "cold_common", mc_bank_id = "b",
    cache = cache, control = ctrl
  )
  expect_false(isTRUE(e1$cache_hit))
  expect_true(isTRUE(e2$cache_hit))
  expect_equal(e1$global_gcv, e2$global_gcv)
  # init cores unchanged
  expect_equal(as.numeric(init[[1]]), as.numeric(init_copy[[1]]))
})

test_that("Phase1 cold grid invariant to traverse order (approx)", {
  skip_on_cran()
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    "smooth_smooth", n = 40L, n_test = 20L, k = 5L, seed = 21L
  )
  probes <- TTPsplines:::.tt_lab_rademacher_probes(des$n, 4L, probe_seed = 5L)
  init <- TTPsplines:::.tt_with_preserved_seed({
    tt_initialize(d = 2L, rank = 3L, k = des$k, seed = 21L)
  })
  ctrl <- TTPsplines:::.tt_lab_phase1_control(seed = 21L, max_sweeps = 12L)
  grid <- c(-1, 1)
  a <- TTPsplines:::.tt_lab_eval_lambda_grid(
    des, 3L, grid, probes, M = 2L, init_policy = "cold_common",
    traverse = "row_asc", common_init = init, control = ctrl
  )
  b <- TTPsplines:::.tt_lab_eval_lambda_grid(
    des, 3L, grid, probes, M = 2L, init_policy = "cold_common",
    traverse = "row_desc", common_init = init, control = ctrl
  )
  a2 <- a[order(a$log10_l1, a$log10_l2), ]
  b2 <- b[order(b$log10_l1, b$log10_l2), ]
  expect_equal(a2$global_gcv_tt, b2$global_gcv_tt, tolerance = 1e-8)
})

test_that("Phase1 sufficient rank matches full TP fitted closely", {
  skip_on_cran()
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    "smooth_smooth", n = 60L, n_test = 30L, k = 5L, seed = 31L
  )
  lam <- c(1, 2)
  full <- TTPsplines:::.tt_lab_eval_full_gcv(des, lam)
  ctrl <- TTPsplines:::.tt_lab_phase1_control(seed = 31L, max_sweeps = 40L)
  fit <- TTPsplines:::.tt_with_preserved_seed({
    ttps(
      des$y, des$X, family = gaussian(), rank = des$k, k = des$k,
      lambda = lam, optimizer = "ALS", backend = "R", control = ctrl
    )
  })
  expect_lt(abs(fit$deviance - full$rss) / max(1, full$rss), 1e-5)
  expect_lt(sqrt(mean((fitted(fit) - full$fitted)^2)), 1e-4)
})

test_that("Phase1 nonconverged handling still records NA under on_nonconverged='na'", {
  skip_on_cran()
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    "weak_signal", n = 40L, n_test = 20L, k = 4L, seed = 41L
  )
  probes <- TTPsplines:::.tt_lab_rademacher_probes(des$n, 2L, probe_seed = 1L)
  ctrl <- TTPsplines:::.tt_lab_phase1_control(seed = 41L, max_sweeps = 8L)
  fit <- TTPsplines:::.tt_with_preserved_seed({
    ttps(
      des$y, des$X, family = gaussian(), rank = 2, k = des$k, lambda = 1,
      optimizer = "ALS", backend = "R", control = ctrl
    )
  })
  bad <- function(y_new, fit_base, ...) {
    out <- fit_base
    out$converged <- FALSE
    out$fitted.values <- fitted(fit_base)
    out$lambda <- fit_base$lambda
    out$rank <- fit_base$rank
    out
  }
  g <- TTPsplines:::tt_global_gdf_mc(
    fit, y = des$y, M = 2L, probes = probes, scheme = "forward",
    refit_fun = bad, on_nonconverged = "na"
  )
  expect_true(all(is.na(g$gdf_contributions)))
})

test_that("Phase1 does not change public ttps formals", {
  f <- formals(ttps)
  expect_true("lambda" %in% names(f))
  expect_true("optimizer" %in% names(f))
  expect_false("global_gcv" %in% names(f))
  expect_false("tt_global_gcv" %in% getNamespaceExports("TTPsplines"))
})
