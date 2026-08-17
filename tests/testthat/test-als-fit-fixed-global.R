# P3 gate: multi-sweep fixed-λ ALS fitter, R vs C++

.rel_err <- function(a, b) abs(as.numeric(a) - as.numeric(b)) / (1 + abs(as.numeric(b)))
.max_abs <- function(A, B) max(abs(as.numeric(A) - as.numeric(B)))

.make_basis_cores <- function(n, d, rank, k_basis = 6L, seed = 1L) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  knots <- lapply(seq_len(d), function(j) make_knots(X[, j], k = k_basis, degree = 3L))
  basis <- eval_marginal_bases(X, knots, degree = 3L)
  ranks <- tt_rank(rank, d = d)
  cores <- initialize_tt_cores(ncol(basis[[1]]), ranks, seed = seed, sd = 0.25)
  y <- as.numeric(tt_contraction(cores, basis) + rnorm(n, 0, 0.15))
  list(y = y, basis = basis, cores = cores, X = X)
}

.expect_fit_parity <- function(ref, cpp, tol_obj = 1e-8, tol_f = 1e-8,
                               tol_core = 1e-6) {
  expect_equal(ref$n_sweeps, cpp$n_sweeps)
  expect_equal(ref$converged, cpp$converged)
  expect_lt(.rel_err(ref$objective, cpp$objective), tol_obj)
  expect_lt(.rel_err(ref$rss, cpp$rss), tol_obj)
  expect_lt(.rel_err(ref$penalty, cpp$penalty), tol_obj)
  expect_lt(abs(ref$intercept - cpp$intercept), tol_f)
  expect_lt(.max_abs(ref$f, cpp$f), tol_f)
  expect_lt(.max_abs(ref$eta, cpp$eta), tol_f)
  for (j in seq_along(ref$cores)) {
    expect_lt(.max_abs(ref$cores[[j]], cpp$cores[[j]]), tol_core)
  }
}

test_that("P3: R vs Rcpp multi-sweep fixed-λ fitter", {
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"),
              "Rcpp fitter not compiled")

  cfgs <- list(
    list(d = 2L, rank = 1L, seed = 61L),
    list(d = 2L, rank = 2L, seed = 62L),
    list(d = 3L, rank = 2L, seed = 63L),
    list(d = 5L, rank = 2L, seed = 64L)
  )
  for (cfg in cfgs) {
    dat <- .make_basis_cores(n = 80L, d = cfg$d, rank = cfg$rank, seed = cfg$seed)
    snap <- lapply(dat$cores, function(C) array(as.numeric(C), dim(C)))
    lam <- if (cfg$d == 2L) c(0.7, 1.4) else pmax(0.1, 10^seq(-0.5, 0.5, length.out = cfg$d))
    ref <- tt_als_fit_fixed_global(
      dat$y, dat$cores, dat$basis, lam,
      max_sweeps = 12L, tol = 1e-8, backend = "R"
    )
    cpp <- tt_als_fit_fixed_global(
      dat$y, dat$cores, dat$basis, lam,
      max_sweeps = 12L, tol = 1e-8, backend = "Rcpp"
    )
    for (j in seq_len(cfg$d)) {
      expect_equal(as.numeric(dat$cores[[j]]), as.numeric(snap[[j]]))
    }
    .expect_fit_parity(ref, cpp)
  }
})

test_that("P3: warm start + RTL + weights", {
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"),
              "Rcpp fitter not compiled")
  dat <- .make_basis_cores(n = 60L, d = 3L, rank = 2L, seed = 71L)
  set.seed(71)
  w <- runif(length(dat$y), 0.4, 1.5)
  lam <- c(1, 2, 0.5)
  ref <- tt_als_fit_fixed_global(
    dat$y, dat$cores, dat$basis, lam,
    weights = w, margin_order = 3:1,
    max_sweeps = 10L, tol = 1e-8, backend = "R"
  )
  cpp <- tt_als_fit_fixed_global(
    dat$y, dat$cores, dat$basis, lam,
    weights = w, margin_order = 3:1,
    max_sweeps = 10L, tol = 1e-8, backend = "Rcpp"
  )
  .expect_fit_parity(ref, cpp)
})

test_that("P3: offset + convergence / max_sweeps", {
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"),
              "Rcpp fitter not compiled")
  dat <- .make_basis_cores(n = 50L, d = 2L, rank = 2L, seed = 81L)
  off <- rep(0.25, length(dat$y))
  lam <- c(1, 1)

  fit_tol <- tt_als_fit_fixed_global(
    dat$y, dat$cores, dat$basis, lam,
    offset = off, max_sweeps = 40L, tol = 1e-6, backend = "Rcpp"
  )
  expect_true(fit_tol$n_sweeps >= 3L)
  expect_true(fit_tol$converged || fit_tol$n_sweeps == 40L)

  fit_cap <- tt_als_fit_fixed_global(
    dat$y, dat$cores, dat$basis, lam,
    offset = off, max_sweeps = 2L, tol = 1e-16, backend = "Rcpp"
  )
  expect_equal(fit_cap$n_sweeps, 2L)
  expect_false(fit_cap$converged)
  expect_identical(fit_cap$convergence_reason, "max_sweeps")

  ref <- tt_als_fit_fixed_global(
    dat$y, dat$cores, dat$basis, lam,
    offset = off, max_sweeps = 15L, tol = 1e-8, backend = "R"
  )
  cpp <- tt_als_fit_fixed_global(
    dat$y, dat$cores, dat$basis, lam,
    offset = off, max_sweeps = 15L, tol = 1e-8, backend = "Rcpp"
  )
  .expect_fit_parity(ref, cpp)
})

test_that("P3: agrees with tt_als_fit_sequential (legacy Gram)", {
  skip_if_not(exists("tt_als_fit_fixed_global_cpp", mode = "function"),
              "Rcpp fitter not compiled")
  dat <- .make_basis_cores(n = 70L, d = 3L, rank = 2L, seed = 91L)
  lam <- c(0.9, 1.1, 1.3)
  ranks <- tt_rank(2L, d = 3L)
  ctrl <- tt_control(
    max_sweeps = 15L, tol = 1e-8, compute_edf = FALSE,
    gram_method = "legacy", design_interface_cache = FALSE, seed = 91L
  )
  seq_fit <- tt_als_fit_sequential(
    dat$y, dat$basis, ranks,
    lambda_spec = list(method = "fixed", values = lam, automatic = FALSE),
    control = ctrl, penalty_order = 2L, init_cores = dat$cores
  )
  cpp <- tt_als_fit_fixed_global(
    dat$y, dat$cores, dat$basis, lam,
    max_sweeps = 15L, tol = 1e-8, backend = "Rcpp"
  )
  expect_lt(.rel_err(seq_fit$history[[length(seq_fit$history)]]$objective,
                     cpp$objective), 1e-6)
  expect_lt(.max_abs(seq_fit$eta, cpp$eta) / (1 + max(abs(seq_fit$eta))), 1e-5)
})
