# P1 gate: single-core fixed-λ update, R vs C++ (global P_k^full)

.rel_err <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  abs(a - b) / (1 + abs(b))
}

.max_abs <- function(A, B) max(abs(as.numeric(A) - as.numeric(B)))

.expect_core_parity <- function(ref, cpp, tol_mat = 1e-8, tol_obj = 1e-8) {
  expect_lt(.max_abs(ref$S, cpp$S), tol_mat)
  expect_lt(.max_abs(ref$b, cpp$b), tol_mat)
  expect_lt(.max_abs(ref$P_own, cpp$P_own), tol_mat)
  expect_lt(.max_abs(ref$P_other, cpp$P_other), tol_mat)
  expect_lt(.max_abs(ref$P_full, cpp$P_full), tol_mat)
  expect_lt(.max_abs(ref$g, cpp$g), tol_mat)
  expect_lt(.rel_err(ref$rss, cpp$rss), tol_obj)
  expect_lt(.rel_err(ref$penalty, cpp$penalty), tol_obj)
  expect_lt(.rel_err(ref$objective, cpp$objective), tol_obj)
  # Input cores must not be mutated: compare via separate clone check in caller
}

.make_basis_cores <- function(n, d, rank, k_basis = 6L, seed = 1L,
                              penalty_order = 2L) {
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
  list(X = X, y = y, basis = basis, cores = cores, ranks = ranks,
       knots = knots, p = p)
}

test_that("P1: R vs Rcpp single-core update (iso/aniso, edge/interior)", {
  skip_if_not(exists("tt_als_core_update_global_cpp", mode = "function"),
              "Rcpp core-update not compiled")

  cfgs <- list(
    list(d = 2L, rank = 1L, seed = 11L),
    list(d = 2L, rank = 2L, seed = 12L),
    list(d = 3L, rank = 2L, seed = 13L),
    list(d = 3L, rank = 3L, seed = 14L),
    list(d = 5L, rank = 2L, seed = 15L)
  )

  for (cfg in cfgs) {
    dat <- .make_basis_cores(n = 80L, d = cfg$d, rank = cfg$rank, seed = cfg$seed)
    intercept <- mean(dat$y)
    d <- cfg$d
    # isotropic and anisotropic λ
    lam_list <- list(
      rep(1.0, d),
      pmax(0.05, 10^seq(-1, 1, length.out = d))
    )
    # edge + interior margins
    ks <- unique(c(1L, as.integer(ceiling(d / 2)), d))
    for (lam in lam_list) {
      for (k in ks) {
        cores_snap <- lapply(dat$cores, function(C) array(as.numeric(C), dim(C)))
        ref <- tt_als_core_update_global(
          dat$y, dat$cores, intercept, dat$basis, k, lam, backend = "R"
        )
        cpp <- tt_als_core_update_global(
          dat$y, dat$cores, intercept, dat$basis, k, lam, backend = "Rcpp"
        )
        # no mutation of inputs
        for (j in seq_len(d)) {
          expect_equal(as.numeric(dat$cores[[j]]), as.numeric(cores_snap[[j]]))
        }
        .expect_core_parity(ref, cpp)
      }
    }
  }
})

test_that("P1: near-singular / extreme λ still match R", {
  skip_if_not(exists("tt_als_core_update_global_cpp", mode = "function"),
              "Rcpp core-update not compiled")

  dat <- .make_basis_cores(n = 60L, d = 3L, rank = 2L, k_basis = 8L, seed = 42L)
  intercept <- mean(dat$y)
  # very small and very large λ (conditioning stress)
  for (lam in list(rep(1e-8, 3), rep(1e6, 3), c(1e-6, 1e2, 1e-4))) {
    for (k in 1:3) {
      ref <- tt_als_core_update_global(
        dat$y, dat$cores, intercept, dat$basis, k, lam, backend = "R"
      )
      cpp <- tt_als_core_update_global(
        dat$y, dat$cores, intercept, dat$basis, k, lam, backend = "Rcpp"
      )
      # allow slightly looser g on extreme conditioning
      expect_lt(.max_abs(ref$P_full, cpp$P_full), 1e-7)
      expect_lt(.rel_err(ref$objective, cpp$objective), 1e-6)
      expect_lt(.max_abs(ref$g, cpp$g) / (1 + max(abs(ref$g))), 1e-5)
    }
  }
})

test_that("P1: observation weights parity", {
  skip_if_not(exists("tt_als_core_update_global_cpp", mode = "function"),
              "Rcpp core-update not compiled")
  dat <- .make_basis_cores(n = 50L, d = 2L, rank = 2L, seed = 7L)
  set.seed(7)
  w <- runif(length(dat$y), 0.2, 2)
  intercept <- mean(dat$y)
  lam <- c(0.7, 1.3)
  ref <- tt_als_core_update_global(
    dat$y, dat$cores, intercept, dat$basis, k = 1L, lambda = lam,
    weights = w, backend = "R"
  )
  cpp <- tt_als_core_update_global(
    dat$y, dat$cores, intercept, dat$basis, k = 1L, lambda = lam,
    weights = w, backend = "Rcpp"
  )
  .expect_core_parity(ref, cpp)
})
