# P4B: safe probe_warm GDF + adaptive fidelity metadata

test_that("P4B: probe_warm does not mutate fit_base cores", {
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"))

  set.seed(201)
  n <- 50L
  d <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
  ctrl <- tt_control(max_sweeps = 10L, tol = 1e-8, compute_edf = FALSE, seed = 201L)
  fit <- TTPsplines:::.tt_lab_fit_fixed(
    y, X, lambda = c(1, 1.2), rank = 2L, control = ctrl,
    fit_backend = "Rcpp_fixed", k = 5L
  )
  snap <- lapply(fit$cores, function(C) array(as.numeric(C), dim(C)))
  probes <- TTPsplines:::.tt_lab_rademacher_probes(n, M = 3L, probe_seed = 9L)

  gdf <- TTPsplines:::tt_global_gdf_mc(
    fit_base = fit, y = y, probes = probes, M = 3L,
    gdf_init = "probe_warm", on_nonconverged = "na", control = ctrl
  )
  expect_identical(gdf$gdf_init, "probe_warm")
  expect_true(isTRUE(gdf$base_cores_unchanged))
  for (j in seq_along(fit$cores)) {
    expect_equal(as.numeric(fit$cores[[j]]), as.numeric(snap[[j]]))
  }
})

test_that("P4B: GDF restores RNG state", {
  set.seed(202)
  n <- 40L
  d <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- rnorm(n)
  ctrl <- tt_control(max_sweeps = 6L, compute_edf = FALSE, seed = 202L)
  fit <- TTPsplines:::.tt_lab_fit_fixed(
    y, X, lambda = c(1, 1), rank = 2L, control = ctrl,
    fit_backend = "R", k = 5L
  )
  probes <- TTPsplines:::.tt_lab_rademacher_probes(n, M = 2L, probe_seed = 3L)
  set.seed(999)
  seed_before <- .Random.seed
  invisible(TTPsplines:::tt_global_gdf_mc(
    fit_base = fit, y = y, probes = probes, M = 2L,
    gdf_init = "cold", on_nonconverged = "na", control = ctrl
  ))
  expect_identical(.Random.seed, seed_before)
})

test_that("P4B: cold vs probe_warm both finite; metadata present", {
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"))
  set.seed(203)
  n <- 55L
  d <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- sin(pi * X[, 1]) + 0.2 * X[, 2] + rnorm(n, 0, 0.15)
  probes <- TTPsplines:::.tt_lab_rademacher_probes(n, M = 3L, probe_seed = 11L)
  ctrl <- tt_control(max_sweeps = 10L, tol = 1e-8, compute_edf = FALSE, seed = 203L)
  init <- tt_initialize(d = d, rank = 2L, k = 5L, seed = 203L)

  g_w <- TTPsplines:::tt_global_gcv(
    lambda = c(1, 1), y = y, X = X, rank = 2L, probes = probes,
    M = 3L, control = ctrl, init = init, fit_backend = "Rcpp_fixed",
    k = 5L, gdf_init = "probe_warm", fidelity = "sobol",
    on_nonconverged = "na"
  )
  g_c <- TTPsplines:::tt_global_gcv(
    lambda = c(1, 1), y = y, X = X, rank = 2L, probes = probes,
    M = 3L, control = ctrl, init = TTPsplines:::.tt_clone_cores(init),
    fit_backend = "Rcpp_fixed", k = 5L, gdf_init = "cold",
    fidelity = "sobol", on_nonconverged = "na"
  )
  expect_identical(g_w$gdf_init, "probe_warm")
  expect_identical(g_c$gdf_init, "cold")
  expect_identical(g_w$fidelity, "sobol")
  expect_true(is.finite(g_w$gdf) || is.finite(g_c$gdf))
  expect_true(isTRUE(g_w$gdf_result$base_cores_unchanged))
})

test_that("P4B: stage fidelity helper applies budgets", {
  base <- tt_control(max_sweeps = 40L, tol = 1e-8, compute_edf = FALSE)
  s <- TTPsplines:::.tt_lab_control_for_stage(base, "sobol", adaptive = TRUE)
  r <- TTPsplines:::.tt_lab_control_for_stage(base, "refine", adaptive = TRUE)
  f <- TTPsplines:::.tt_lab_control_for_stage(base, "final", adaptive = TRUE)
  expect_equal(s$max_sweeps, 12L)
  expect_equal(r$max_sweeps, 25L)
  expect_equal(f$max_sweeps, 50L)
  off <- TTPsplines:::.tt_lab_control_for_stage(base, "sobol", adaptive = FALSE)
  expect_equal(off$max_sweeps, 40L)
})

test_that("P4B: optimize records fidelity and winner_source", {
  skip_on_cran()
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"))
  set.seed(204)
  n <- 80L
  d <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- sin(2 * pi * X[, 1]) + 0.3 * cos(2 * pi * X[, 2]) + rnorm(n, 0, 0.2)
  out <- TTPsplines:::tt_global_lambda_optimize(
    y = y, X = X, rank = 2L,
    theta_lower = -2, theta_upper = 2,
    n_global = 8L, n_refine = 1L, n_diverse = 3L,
    M_search = 3L, M_final = 5L,
    core_starts_final = 1L,
    include_cgcv_anchor = FALSE,
    k = 5L,
    control = tt_control(max_sweeps = 20L, tol = 1e-8, compute_edf = FALSE, seed = 204L),
    fit_backend = "Rcpp_fixed",
    adaptive_fidelity = TRUE,
    gdf_init = "probe_warm",
    verbose = FALSE
  )
  expect_true(isTRUE(out$diagnostics$adaptive_fidelity))
  expect_identical(out$diagnostics$gdf_init, "probe_warm")
  expect_true(!is.null(out$diagnostics$fidelity$sobol_sweeps))
  expect_true("fidelity" %in% names(out$sobol_results))
  expect_true("winner_source" %in% names(out$sobol_results))
  expect_true("fidelity" %in% names(out$final_candidates))
  expect_true(nzchar(out$winner_source))
  expect_true(all(out$final_candidates$fidelity == "final"))
})
