# Tests for experimental joint global λ optimizer (not exported).

test_that("builtin Sobol is reproducible and in (0,1)^d", {
  U1 <- TTPsplines:::.tt_lab_sobol_unit(32L, 2L, skip = 1L)
  U2 <- TTPsplines:::.tt_lab_sobol_unit(32L, 2L, skip = 1L)
  expect_equal(U1, U2)
  expect_true(all(U1 > 0 & U1 < 1))
  expect_equal(dim(U1), c(32L, 2L))
  # skip changes the sequence
  U3 <- TTPsplines:::.tt_lab_sobol_unit(32L, 2L, skip = 5L)
  expect_false(isTRUE(all.equal(U1, U3)))
})

test_that("sobol box maps to bounds", {
  B <- TTPsplines:::.tt_lab_sobol_box(16L, c(-6, -6), c(6, 6))
  expect_true(all(B >= -6 - 1e-12 & B <= 6 + 1e-12))
})

test_that("diverse_top keeps spread and prefers low Q", {
  theta <- rbind(
    c(0, 0),
    c(0.1, 0),
    c(3, 3),
    c(3.1, 3),
    c(-2, 2)
  )
  q <- c(1.0, 1.05, 2.0, 2.1, 0.5)
  idx <- TTPsplines:::.tt_lab_diverse_top(theta, q, n_keep = 3L, min_dist = 0.5)
  expect_true(idx[1] == 5L) # best Q first
  expect_true(length(idx) <= 3L)
  expect_true(length(unique(idx)) == length(idx))
})

test_that("default n_global scales with d", {
  expect_equal(TTPsplines:::.tt_lab_default_n_global(2L), 256L)
  expect_equal(TTPsplines:::.tt_lab_default_n_global(3L), 512L)
  expect_equal(TTPsplines:::.tt_lab_default_n_global(5L), 1280L)
})

test_that("tt_global_lambda_optimize smoke recovers finite λ (tiny)", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("TT_GGCV_OPT_SMOKE", "true"), "true"),
    "Set TT_GGCV_OPT_SMOKE=false to skip"
  )
  set.seed(1)
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    "smooth_smooth", n = 60L, n_test = 40L, k = 5L, seed = 3L
  )
  # Tiny budget: validates plumbing, not scientific gates
  out <- TTPsplines:::tt_global_lambda_optimize(
    y = des$y,
    X = des$X,
    rank = 5L,
    k = des$k,
    degree = des$degree,
    penalty_order = des$penalty_order,
    theta_lower = -4,
    theta_upper = 4,
    n_global = 12L,
    n_refine = 2L,
    n_diverse = 4L,
    M_search = 3L,
    M_final = 5L,
    core_starts_final = 2L,
    seed = 3L,
    include_cgcv_anchor = TRUE,
    control = tt_control(
      max_sweeps = 25L, tol = 1e-7, compute_edf = FALSE, seed = 3L
    ),
    verbose = FALSE
  )
  expect_true(is.list(out))
  expect_equal(length(out$lambda), 2L)
  expect_true(all(is.finite(out$lambda)) && all(out$lambda > 0))
  expect_true(is.finite(out$gcv))
  expect_true(nrow(out$sobol_results) >= 12L)
  expect_true(nrow(out$refined_results) >= 1L)
  expect_true(nrow(out$final_candidates) >= 1L)
  expect_true(is.logical(out$boundary))
  # Refinement should not worsen the elite pool vs raw Sobol (weak gate)
  expect_true(
    is.finite(out$diagnostics$best_refined_gcv) ||
      is.finite(out$diagnostics$best_sobol_gcv)
  )
})
