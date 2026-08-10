test_that("spectral cGCV matches same-ridge dense path (g, ed, rss, value)", {
  set.seed(101)
  n <- 200
  m <- 25
  X <- matrix(rnorm(n * m), n, m)
  z <- rnorm(n)
  S <- crossprod(X)
  b <- as.numeric(crossprod(X, z))
  P_own <- crossprod(diff(diag(m), differences = 2))
  P_other <- 0.3 * crossprod(diff(diag(m), differences = 1))
  yw <- z
  Xw <- X

  cache <- TTPsplines:::.make_gcv_spectral_cache(
    S, P_own, b, yw = yw, P0 = P_other
  )
  expect_false(is.null(cache))
  expect_equal(length(cache$ed_weights), m)
  expect_equal(length(cache$transformed_b), m)

  # Reference: same fixed ridge baked into A as the spectral factorization.
  same_ridge_dense <- function(lam) {
    M <- S + P_other + cache$ridge * diag(m) + lam * P_own
    g <- as.numeric(TTPsplines:::solve_spd(M, b))
    ed <- sum(diag(TTPsplines:::solve_spd(M, S)))
    rss <- as.numeric(sum(yw^2) - 2 * crossprod(b, g) + crossprod(g, S %*% g))
    rss <- max(rss, 0)
    denom <- (n - ed)^2
    value <- if (!is.finite(denom) || denom < 1e-12) Inf else n * rss / denom
    list(g = g, ed = ed, rss = rss, value = value)
  }

  lams <- exp(seq(log(1e-4), log(1e4), length.out = 12))
  for (lam in lams) {
    dense <- same_ridge_dense(lam)
    spec <- TTPsplines:::.conditional_gcv_spectral(
      yw, Xw, S, P_own, b, lam, cache, P0 = P_other
    )
    expect_equal(spec$g, dense$g, tolerance = 1e-8)
    expect_equal(spec$ed, dense$ed, tolerance = 1e-8)
    expect_equal(spec$rss, dense$rss, tolerance = 1e-8)
    expect_equal(spec$value, dense$value, tolerance = 1e-8)
    # Also close to package .ed_S / .conditional_gcv (adaptive ridge)
    pack <- TTPsplines:::.conditional_gcv(
      yw, Xw, S, P_own, b, lam, P0 = P_other
    )
    expect_equal(spec$ed, pack$ed, tolerance = 5e-4)
    expect_equal(spec$value, pack$value, tolerance = 5e-4)
  }
})

test_that("make_core_workspace attaches spectral once; Brent reuses it", {
  set.seed(102)
  n <- 120
  m <- 18
  X <- matrix(rnorm(n * m), n, m)
  z <- rnorm(n)
  P <- crossprod(diff(diag(m), differences = 2))
  P0 <- 0.2 * diag(m)
  ws <- TTPsplines:::make_core_workspace(
    z, X, P, lambda0 = 1,
    bounds = c(1e-4, 1e4), tol = 1e-6,
    use_spectral = TRUE, P0 = P0
  )
  expect_true(!is.null(ws$spectral))
  cache_ptr <- ws$spectral
  upd <- TTPsplines:::update_lambda_cgcv(ws)
  expect_true(all(is.finite(c(upd$lambda, upd$value, upd$ed))))
  expect_identical(ws$spectral, cache_ptr)

  grid <- exp(seq(log(1e-3), log(1e3), length.out = 41))
  vals <- vapply(
    grid,
    function(lam) TTPsplines:::.cgcv_eval_at(ws, lam)$value,
    numeric(1)
  )
  expect_true(all(is.finite(vals) | is.infinite(vals)))
  expect_identical(ws$spectral, cache_ptr)
})

test_that("spectral EDF uses tr((A+λB)^{-1}S), not tr((A+λB)^{-1}A)", {
  set.seed(103)
  m <- 12
  n <- 80
  X <- matrix(rnorm(n * m), n, m)
  z <- rnorm(n)
  S <- crossprod(X)
  b <- as.numeric(crossprod(X, z))
  P_own <- crossprod(diff(diag(m), differences = 1))
  P_other <- 2 * diag(m)
  cache <- TTPsplines:::.make_gcv_spectral_cache(
    S, P_own, b, yw = z, P0 = P_other
  )
  lam <- 3.5
  spec <- TTPsplines:::.conditional_gcv_spectral(
    z, X, S, P_own, b, lam, cache, P0 = P_other
  )
  wrong <- sum(1 / (1 + lam * cache$values))
  expect_true(abs(spec$ed - wrong) > 1e-3)
  M <- S + P_other + cache$ridge * diag(m) + lam * P_own
  ed_ref <- sum(diag(TTPsplines:::solve_spd(M, S)))
  expect_equal(spec$ed, ed_ref, tolerance = 1e-8)
})

test_that("spectral frozen grid is faster than dense refactorization path", {
  skip_on_cran()
  set.seed(104)
  n <- 2000
  m <- 60
  X <- matrix(rnorm(n * m), n, m)
  z <- rnorm(n)
  P <- crossprod(diff(diag(m), differences = 2))
  P0 <- 0.5 * crossprod(diff(diag(m), differences = 1))
  grid <- exp(seq(log(1e-4), log(1e4), length.out = 80))

  ws_spec <- TTPsplines:::make_core_workspace(
    z, X, P, 1, c(1e-4, 1e4), 1e-8,
    use_spectral = TRUE, P0 = P0
  )
  ws_dense <- TTPsplines:::make_core_workspace(
    z, X, P, 1, c(1e-4, 1e4), 1e-8,
    use_spectral = FALSE, P0 = P0
  )
  expect_true(!is.null(ws_spec$spectral))

  t_spec <- system.time({
    for (lam in grid) TTPsplines:::.cgcv_eval_at(ws_spec, lam)
  })[["elapsed"]]
  t_dense <- system.time({
    for (lam in grid) TTPsplines:::.cgcv_eval_at(ws_dense, lam)
  })[["elapsed"]]

  expect_lt(t_spec, t_dense)
})
