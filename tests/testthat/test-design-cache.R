library(TTPsplines)

test_that("LTR design cache matches full left/right rebuild (Xk)", {
  set.seed(11)
  n <- 80L; d <- 4L; k <- 6L; r <- 2L
  X <- matrix(runif(n * d), n, d)
  bs <- TTPsplines:::build_marginal_bases(X, k = k, degree = 3)
  basis <- bs$basis
  ranks <- c(1L, rep(r, d - 1L), 1L)
  cores <- TTPsplines:::initialize_tt_cores(k, ranks, seed = 2, sd = 0.1)

  R_all <- TTPsplines:::.tt_design_prepare_right(cores, basis)
  L_cur <- matrix(1, n, 1)
  for (kk in seq_len(d)) {
    X_cached <- TTPsplines:::tt_design_core(L_cur, R_all[[kk]], basis[[kk]])
    L_full <- TTPsplines:::left_interfaces(cores, basis)
    R_full <- TTPsplines:::right_interfaces(cores, basis)
    X_ref <- TTPsplines:::tt_design_core(L_full[[kk]], R_full[[kk]], basis[[kk]])
    expect_lt(max(abs(X_cached - X_ref)), 1e-10)

    cores[[kk]] <- cores[[kk]] + array(rnorm(length(cores[[kk]]), 0, 0.01),
                                       dim = dim(cores[[kk]]))
    if (kk < d) {
      L_cur <- TTPsplines:::.tt_design_left_absorb(L_cur, cores[[kk]], basis[[kk]])
    }
  }
})

test_that("RTL design cache matches full rebuild (Xk)", {
  set.seed(12)
  n <- 60L; d <- 3L; k <- 5L; r <- 2L
  X <- matrix(runif(n * d), n, d)
  bs <- TTPsplines:::build_marginal_bases(X, k = k, degree = 3)
  basis <- bs$basis
  ranks <- c(1L, rep(r, d - 1L), 1L)
  cores <- TTPsplines:::initialize_tt_cores(k, ranks, seed = 4, sd = 0.1)

  L_all <- TTPsplines:::.tt_design_prepare_left(cores, basis)
  R_cur <- matrix(1, n, 1)
  for (kk in d:1) {
    X_cached <- TTPsplines:::tt_design_core(L_all[[kk]], R_cur, basis[[kk]])
    L_full <- TTPsplines:::left_interfaces(cores, basis)
    R_full <- TTPsplines:::right_interfaces(cores, basis)
    X_ref <- TTPsplines:::tt_design_core(L_full[[kk]], R_full[[kk]], basis[[kk]])
    expect_lt(max(abs(X_cached - X_ref)), 1e-10)

    cores[[kk]] <- cores[[kk]] + array(rnorm(length(cores[[kk]]), 0, 0.01),
                                       dim = dim(cores[[kk]]))
    if (kk > 1L) {
      R_cur <- TTPsplines:::.tt_design_right_absorb(R_cur, cores[[kk]], basis[[kk]])
    }
  }
})

test_that("design-cache ALS matches rebuild: eta, objective, fitted", {
  skip_on_cran()
  set.seed(21)
  n <- 180L; d <- 3L; k <- 8L; r <- 2L
  X <- matrix(runif(n * d), n, d)
  f <- sin(2 * pi * X[, 1]) + 0.4 * X[, 2] * X[, 3]
  y <- f + rnorm(n, 0, 0.2)

  base_ctrl <- list(
    max_sweeps = 10L,
    backend = "R",
    compute_edf = FALSE,
    warn_lambda_boundary = FALSE,
    cgcv_trace = FALSE,
    seed = 7L,
    use_spectral_gcv = FALSE
  )
  ctrl_cache <- do.call(tt_control, c(base_ctrl, list(design_interface_cache = TRUE)))
  ctrl_ref <- do.call(tt_control, c(base_ctrl, list(design_interface_cache = FALSE)))

  fit_c <- ttps(y, X, rank = r, k = k, lambda = 1, control = ctrl_cache)
  fit_r <- ttps(y, X, rank = r, k = k, lambda = 1, control = ctrl_ref)

  expect_lt(max(abs(fit_c$linear.predictors - fit_r$linear.predictors)), 1e-10)
  expect_lt(max(abs(fit_c$fitted.values - fit_r$fitted.values)), 1e-10)
  expect_lt(abs(fit_c$deviance - fit_r$deviance), 1e-10)
  # packed cores
  for (kk in seq_along(fit_c$cores)) {
    expect_lt(max(abs(fit_c$cores[[kk]] - fit_r$cores[[kk]])), 1e-10)
  }
  expect_lt(abs(fit_c$intercept - fit_r$intercept), 1e-10)
})

test_that("design-cache preserves XtWX / XtWz vs rebuild at frozen cores", {
  set.seed(31)
  n <- 120L; d <- 3L; k <- 7L; r <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- rnorm(n)
  bs <- TTPsplines:::build_marginal_bases(X, k = k, degree = 3)
  basis <- bs$basis
  ranks <- c(1L, rep(r, d - 1L), 1L)
  cores <- TTPsplines:::initialize_tt_cores(k, ranks, seed = 5, sd = 0.12)
  w <- rep(1, n)
  zc <- y - mean(y)
  lambda <- rep(1, d)
  ctrl <- tt_control(compute_edf = FALSE, use_spectral_gcv = FALSE,
                     gram_method = "legacy")

  L_all <- TTPsplines:::left_interfaces(cores, basis)
  R_all <- TTPsplines:::right_interfaces(cores, basis)
  for (kk in seq_len(d)) {
    built_c <- TTPsplines:::.cgcv_core_workspace(
      cores, kk, lambda, basis, zc, ranks, ctrl,
      weight = w, Left = L_all[[kk]], Right = R_all[[kk]],
      return_design = TRUE
    )
    built_r <- TTPsplines:::.cgcv_core_workspace(
      cores, kk, lambda, basis, zc, ranks, ctrl, weight = w,
      return_design = TRUE
    )
    expect_lt(max(abs(built_c$Xk - built_r$Xk)), 1e-10)
    Sc <- crossprod(built_c$Xk, built_c$Xk * w)
    Sr <- crossprod(built_r$Xk, built_r$Xk * w)
    bc <- as.numeric(crossprod(built_c$Xk, w * zc))
    br <- as.numeric(crossprod(built_r$Xk, w * zc))
    expect_lt(max(abs(Sc - Sr)), 1e-8)
    expect_lt(max(abs(bc - br)), 1e-8)
  }
})

test_that("fused Gram ALS matches legacy design+crossprod (eta gate)", {
  skip_on_cran()
  set.seed(55)
  n <- 160L; d <- 3L; k <- 7L; r <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2] * X[, 3] + rnorm(n, 0, 0.2)
  base <- list(
    max_sweeps = 8L, backend = "R", compute_edf = FALSE,
    warn_lambda_boundary = FALSE, cgcv_trace = FALSE, seed = 11L,
    use_spectral_gcv = FALSE, design_interface_cache = TRUE
  )
  fit_f <- ttps(
    y, X, rank = r, k = k, lambda = 1,
    control = do.call(tt_control, c(base, list(gram_method = "fused_blocked")))
  )
  fit_l <- ttps(
    y, X, rank = r, k = k, lambda = 1,
    control = do.call(tt_control, c(base, list(gram_method = "legacy")))
  )
  expect_lt(max(abs(fit_f$linear.predictors - fit_l$linear.predictors)), 1e-8)
  expect_lt(abs(fit_f$deviance - fit_l$deviance), 1e-8)
})

test_that("cGCV ALS cache vs rebuild stays within eta gate", {
  skip_on_cran()
  set.seed(41)
  n <- 120L; d <- 3L; k <- 6L; r <- 2L
  X <- matrix(runif(n * d), n, d)
  y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.25)
  base <- list(
    max_sweeps = 6L, backend = "R", compute_edf = FALSE,
    warn_lambda_boundary = FALSE, cgcv_trace = FALSE, seed = 9L,
    use_spectral_gcv = TRUE, cgcv_update = "sequential"
  )
  fit_c <- ttps(
    y, X, rank = r, k = k, lambda = "cGCV",
    control = do.call(tt_control, c(base, list(design_interface_cache = TRUE)))
  )
  fit_r <- ttps(
    y, X, rank = r, k = k, lambda = "cGCV",
    control = do.call(tt_control, c(base, list(design_interface_cache = FALSE)))
  )
  expect_lt(max(abs(fit_c$linear.predictors - fit_r$linear.predictors)), 1e-10)
  expect_lt(max(abs(as.numeric(fit_c$lambda) - as.numeric(fit_r$lambda))), 1e-10)
})
