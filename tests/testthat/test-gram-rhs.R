library(TTPsplines)

test_that("tt_gram_rhs methods match BLAS reference (S, b)", {
  skip_if_not(exists("tt_gram_rhs_cpp", envir = asNamespace("TTPsplines"),
                     inherits = FALSE))
  set.seed(42)
  n <- 80L
  rl <- 2L; p <- 6L; rr <- 2L
  Left <- matrix(rnorm(n * rl), n, rl)
  Right <- matrix(rnorm(n * rr), n, rr)
  Bk <- matrix(runif(n * p), n, p)
  z <- rnorm(n)
  w <- runif(n, 0.2, 1.5)

  ref <- TTPsplines:::tt_gram_rhs(Left, Right, Bk, z, w, method = "blas")
  for (m in c("fused", "fused_blocked", "kron")) {
    got <- TTPsplines:::tt_gram_rhs(Left, Right, Bk, z, w, method = m)
    expect_lt(max(abs(got$S - ref$S)), 1e-9)
    expect_lt(max(abs(got$b - ref$b)), 1e-9)
  }

  ref0 <- TTPsplines:::tt_gram_rhs(Left, Right, Bk, z, NULL, method = "blas")
  got0 <- TTPsplines:::tt_gram_rhs(Left, Right, Bk, z, NULL, method = "fused")
  expect_lt(max(abs(got0$S - ref0$S)), 1e-9)
})

test_that("Kronecker identity matches design row outer product", {
  skip_if_not(exists("tt_gram_rhs_cpp", envir = asNamespace("TTPsplines"),
                     inherits = FALSE))
  set.seed(7)
  rl <- 3L; p <- 4L; rr <- 2L
  L <- rnorm(rl); B <- rnorm(p); R <- rnorm(rr)
  x <- as.numeric(kronecker(R, kronecker(B, L)))
  S1 <- tcrossprod(x)
  S2 <- kronecker(tcrossprod(R), kronecker(tcrossprod(B), tcrossprod(L)))
  expect_lt(max(abs(S1 - S2)), 1e-12)

  Left <- matrix(L, 1, rl)
  Right <- matrix(R, 1, rr)
  Bk <- matrix(B, 1, p)
  X <- TTPsplines:::tt_design_core(Left, Right, Bk)
  expect_lt(max(abs(as.numeric(X) - x)), 1e-12)
})

test_that("OpenMP fused_blocked matches serial (T=1,2,4,8)", {
  skip_if_not(exists("tt_gram_omp_available",
                     envir = asNamespace("TTPsplines"), inherits = FALSE))
  skip_if_not(isTRUE(TTPsplines:::tt_gram_omp_available()))
  set.seed(99)
  n <- 400L; rl <- 2L; p <- 8L; rr <- 2L
  Left <- matrix(rnorm(n * rl), n, rl)
  Right <- matrix(rnorm(n * rr), n, rr)
  Bk <- matrix(runif(n * p), n, p)
  z <- rnorm(n)
  w <- runif(n, 0.5, 1.5)
  ref <- TTPsplines:::tt_gram_rhs(
    Left, Right, Bk, z, w, method = "fused_blocked", n_threads = 1L
  )
  for (T in c(2L, 4L, 8L)) {
    got <- TTPsplines:::tt_gram_rhs(
      Left, Right, Bk, z, w, method = "fused_blocked", n_threads = T
    )
    expect_lt(max(abs(got$S - ref$S)), 1e-8)
    expect_lt(max(abs(got$b - ref$b)), 1e-8)
    expect_true(got$n_threads >= 1L)
  }
})
