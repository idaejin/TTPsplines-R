# P2 gate: one fixed-λ ALS sweep, R vs C++ (global P_k^full)

.rel_err <- function(a, b) abs(as.numeric(a) - as.numeric(b)) / (1 + abs(as.numeric(b)))
.max_abs <- function(A, B) max(abs(as.numeric(A) - as.numeric(B)))

.make_basis_cores <- function(n, d, rank, k_basis = 6L, seed = 1L) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  knots <- lapply(seq_len(d), function(j) {
    make_knots(X[, j], k = k_basis, degree = 3L)
  })
  basis <- eval_marginal_bases(X, knots, degree = 3L)
  ranks <- tt_rank(rank, d = d)
  p <- ncol(basis[[1]])
  cores <- initialize_tt_cores(p, ranks, seed = seed, sd = 0.25)
  y <- as.numeric(tt_contraction(cores, basis) + rnorm(n, 0, 0.15))
  list(y = y, basis = basis, cores = cores)
}

.expect_sweep_parity <- function(ref, cpp, tol_f = 1e-8, tol_obj = 1e-8,
                                 tol_core = 1e-7) {
  expect_lt(.rel_err(ref$objective, cpp$objective), tol_obj)
  expect_lt(.rel_err(ref$rss, cpp$rss), tol_obj)
  expect_lt(.rel_err(ref$penalty, cpp$penalty), tol_obj)
  expect_lt(.max_abs(ref$f, cpp$f), tol_f)
  expect_lt(.max_abs(ref$eta, cpp$eta), tol_f)
  for (j in seq_along(ref$cores)) {
    expect_lt(.max_abs(ref$cores[[j]], cpp$cores[[j]]), tol_core)
  }
  expect_identical(as.integer(ref$margin_order), as.integer(cpp$margin_order))
}

test_that("P2: one fixed-λ sweep LTR parity (d,r,iso/aniso)", {
  skip_if_not(exists("tt_als_sweep_global_cpp", mode = "function"),
              "Rcpp sweep not compiled")

  cfgs <- list(
    list(d = 2L, rank = 1L, seed = 21L),
    list(d = 2L, rank = 2L, seed = 22L),
    list(d = 3L, rank = 2L, seed = 23L),
    list(d = 3L, rank = 3L, seed = 24L),
    list(d = 5L, rank = 2L, seed = 25L)
  )
  for (cfg in cfgs) {
    dat <- .make_basis_cores(n = 70L, d = cfg$d, rank = cfg$rank, seed = cfg$seed)
    intercept <- mean(dat$y)
    d <- cfg$d
    for (lam in list(rep(1.0, d), pmax(0.05, 10^seq(-1, 1, length.out = d)))) {
      snap <- lapply(dat$cores, function(C) array(as.numeric(C), dim(C)))
      ref <- tt_als_sweep_global(
        dat$y, dat$cores, intercept, dat$basis, lam, backend = "R"
      )
      cpp <- tt_als_sweep_global(
        dat$y, dat$cores, intercept, dat$basis, lam, backend = "Rcpp"
      )
      for (j in seq_len(d)) {
        expect_equal(as.numeric(dat$cores[[j]]), as.numeric(snap[[j]]))
      }
      .expect_sweep_parity(ref, cpp)
    }
  }
})

test_that("P2: RTL margin order parity", {
  skip_if_not(exists("tt_als_sweep_global_cpp", mode = "function"),
              "Rcpp sweep not compiled")
  dat <- .make_basis_cores(n = 60L, d = 3L, rank = 2L, seed = 31L)
  intercept <- mean(dat$y)
  lam <- c(0.5, 1.2, 2.0)
  ord <- 3:1
  ref <- tt_als_sweep_global(
    dat$y, dat$cores, intercept, dat$basis, lam,
    margin_order = ord, backend = "R"
  )
  cpp <- tt_als_sweep_global(
    dat$y, dat$cores, intercept, dat$basis, lam,
    margin_order = ord, backend = "Rcpp"
  )
  .expect_sweep_parity(ref, cpp)
})

test_that("P2: weights + extreme λ sweep", {
  skip_if_not(exists("tt_als_sweep_global_cpp", mode = "function"),
              "Rcpp sweep not compiled")
  dat <- .make_basis_cores(n = 55L, d = 2L, rank = 2L, seed = 41L)
  set.seed(41)
  w <- runif(length(dat$y), 0.3, 1.8)
  intercept <- mean(dat$y)
  for (lam in list(c(1e-6, 1e4), c(2.0, 0.5))) {
    ref <- tt_als_sweep_global(
      dat$y, dat$cores, intercept, dat$basis, lam,
      weights = w, backend = "R"
    )
    cpp <- tt_als_sweep_global(
      dat$y, dat$cores, intercept, dat$basis, lam,
      weights = w, backend = "Rcpp"
    )
    expect_lt(.rel_err(ref$objective, cpp$objective), 1e-6)
    expect_lt(.max_abs(ref$f, cpp$f) / (1 + max(abs(ref$f))), 1e-5)
  }
})

test_that("P2: sequential P1 core updates match one P2 sweep (R)", {
  dat <- .make_basis_cores(n = 40L, d = 3L, rank = 2L, seed = 51L)
  intercept <- mean(dat$y)
  lam <- c(0.8, 1.1, 1.5)
  sw <- tt_als_sweep_global(
    dat$y, dat$cores, intercept, dat$basis, lam, backend = "R"
  )
  cores <- dat$cores
  for (k in 1:3) {
    step <- tt_als_core_update_global(
      dat$y, cores, intercept, dat$basis, k, lam, backend = "R"
    )
    cores <- step$cores
  }
  expect_lt(.max_abs(sw$cores[[2]], cores[[2]]), 1e-10)
  expect_lt(.rel_err(sw$objective, tt_gaussian_Q(
    dat$y, cores, intercept, dat$basis, lam
  )$value), 1e-12)
})
