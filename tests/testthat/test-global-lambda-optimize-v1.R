# Tests for adaptive global-λ optimizer v1 (not exported).

test_that("near_optimal_region uses 1% relative band", {
  theta <- rbind(c(0, 0), c(1, 0), c(0, 3), c(10, 10))
  q <- c(1.0, 1.005, 1.02, 5)
  reg <- TTPsplines:::.tt_lab_near_optimal_region(theta, q, alpha = 0.01)
  expect_equal(reg$n_in, 2L) # 1.0 and 1.005
  expect_equal(unname(reg$intervals[1, ]), c(0, 1))
})

test_that("expand_box widens only touching margins", {
  exp <- TTPsplines:::.tt_lab_expand_box(
    lower = c(-5, -5, -5), upper = c(5, 5, 5),
    theta_best = c(0, -4.8, 1),
    region_intervals = rbind(c(-1, 1), c(-5, -4.5), c(0, 2)),
    boundary_tol = 0.5, soft_tol = 1.0, expand_step = 3
  )
  expect_true(2L %in% exp$expanded_dimensions)
  expect_false(1L %in% exp$expanded_dimensions)
  expect_lt(exp$lower[2], -5)
  expect_equal(exp$lower[1], -5)
})

test_that("expand_box respects force_lo from edge probes", {
  # Interior best on θ2, but edge probe says lower face is competitive
  exp <- TTPsplines:::.tt_lab_expand_box(
    lower = c(-5, -5, -5), upper = c(5, 5, 5),
    theta_best = c(0, -3.3, -2),
    region_intervals = rbind(c(-0.5, 0.5), c(-3.6, -3.2), c(-2.2, -1.8)),
    boundary_tol = 0.5, soft_tol = 1.0, expand_step = 3,
    force_lo = c(FALSE, TRUE, FALSE),
    force_hi = c(FALSE, FALSE, FALSE)
  )
  expect_equal(exp$expanded_dimensions, 2L)
  expect_lt(exp$lower[2], -5)
  expect_equal(exp$lower[1], -5)
  expect_equal(exp$lower[3], -5)
})

test_that("gcv_tied respects SE and relative fallback", {
  expect_true(TTPsplines:::.tt_lab_gcv_tied(1.0, 1.005, NA, NA, 0.01))
  expect_false(TTPsplines:::.tt_lab_gcv_tied(1.0, 1.05, NA, NA, 0.01))
  expect_true(TTPsplines:::.tt_lab_gcv_tied(1.0, 1.02, 0.02, 0.02, 0.01))
})

test_that("classify flags flat edge valley as weakly/unpenalized", {
  intervals <- rbind(c(-1, -0.5), c(-8, -4.5), c(-2, -1.5))
  rownames(intervals) <- paste0("theta", 1:3)
  colnames(intervals) <- c("lo", "hi")
  # flat profile on margin 2 over the low-λ plateau; rises only for large λ
  prof2 <- data.frame(
    margin = 2,
    theta_j = c(seq(-8, -4.5, length.out = 6), seq(-3, 5, length.out = 5)),
    q = c(rep(1.0, 6), 1.02, 1.05, 1.1, 1.15, 1.2),
    gdf = 10
  )
  profiles <- list(
    data.frame(margin = 1, theta_j = seq(-5, 5, length.out = 5),
               q = c(2, 1.2, 1, 1.2, 2), gdf = 1),
    prof2,
    data.frame(margin = 3, theta_j = seq(-5, 5, length.out = 5),
               q = c(2, 1.2, 1, 1.2, 2), gdf = 1)
  )
  id <- TTPsplines:::.tt_lab_classify_lambda_id(
    intervals, lower = c(-5, -8, -5), upper = c(5, 5, 5),
    profile_list = profiles, boundary_tol = 0.5
  )
  expect_true(id$labels[["lambda2"]] %in%
                c("weakly_identified", "effectively_unpenalized"))
})

test_that("tt_global_lambda_optimize_v1 smoke d=2", {
  skip_on_cran()
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    "smooth_smooth", n = 50L, n_test = 30L, k = 5L, seed = 2L
  )
  out <- TTPsplines:::tt_global_lambda_optimize_v1(
    y = des$y, X = des$X, rank = 5L, k = 5L,
    theta_lower = -4, theta_upper = 4,
    sobol_batches = c(8L, 8L),
    n_refine = 1L, n_diverse = 3L,
    M_search = 2L, M_final = 3L, core_starts_final = 1L,
    adaptive_box = FALSE,
    profile_lambda = TRUE,
    classify_identifiability = TRUE,
    seed = 2L,
    control = tt_control(max_sweeps = 15L, tol = 1e-6, compute_edf = FALSE,
                         seed = 2L),
    verbose = FALSE
  )
  expect_equal(length(out$lambda), 2L)
  expect_true(all(is.finite(out$lambda)))
  expect_true(!is.null(out$near_optimal_region$intervals))
  expect_true(!is.null(out$lambda_identifiability))
  expect_equal(out$diagnostics$version, "v1")
})

test_that("phase1 design_dn builds d=4 strong_aniso", {
  des <- TTPsplines:::.tt_lab_phase1_make_design_dn(
    d = 4L, scenario = "strong_aniso", n = 40L, n_test = 20L, k = 5L, seed = 1L
  )
  expect_equal(des$d, 4L)
  expect_equal(ncol(des$X), 4L)
  expect_equal(des$rough_margin, 2L)
  expect_equal(length(des$y), 40L)
})
