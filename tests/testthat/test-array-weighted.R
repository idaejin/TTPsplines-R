## Weighted array-mode Gram (GLAM row tensors): X_k' W X_k and X_k' W z
## without forming X_k. Parity vs the dense scattered computation, and
## end-to-end Poisson PIRLS parity array vs scattered.

test_that("tt_gram_rhs_array_weighted matches dense weighted Gram", {
  set.seed(11)
  n_grid <- c(5L, 4L, 3L)
  n <- prod(n_grid)
  rl <- 2L; p <- 4L; rr <- 3L
  # Scattered interfaces with the grid's Kronecker structure (dim1 fastest):
  # for k = 2, Left rows repeat over (n2*n3) blocks, Right over (n1*n2).
  k <- 2L
  n_left  <- n_grid[1]
  n_right <- n_grid[3]
  L_uniq <- matrix(rnorm(n_left * rl), n_left, rl)
  R_uniq <- matrix(rnorm(n_right * rr), n_right, rr)
  Bk_marg <- matrix(rnorm(n_grid[k] * p), n_grid[k], p)
  # Expand to scattered layout (row order: i1 fastest, then i2, then i3).
  idx <- expand.grid(i1 = seq_len(n_grid[1]), i2 = seq_len(n_grid[2]),
                     i3 = seq_len(n_grid[3]), KEEP.OUT.ATTRS = FALSE)
  Left  <- L_uniq[idx$i1, , drop = FALSE]
  Right <- R_uniq[idx$i3, , drop = FALSE]
  Bk_sc <- Bk_marg[idx$i2, , drop = FALSE]
  w <- runif(n, 0.1, 2)
  z <- rnorm(n)

  # Dense reference: X = rowwise (Right x Bk x Left), kron column order.
  X <- TTPsplines:::tt_design_core(Left, Right, Bk_sc)
  S_ref <- crossprod(X * sqrt(w))
  b_ref <- as.numeric(crossprod(X, w * z))

  gr <- TTPsplines:::tt_gram_rhs_array_weighted(
    k = k, Left = Left, Right = Right, Bk = Bk_marg,
    w = w, z = z, n_grid = n_grid, marginal_iface = FALSE
  )
  expect_equal(gr$S, S_ref, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(gr$b, b_ref, tolerance = 1e-12)
  # Also via the marginal-interface entry point.
  gr2 <- TTPsplines:::tt_gram_rhs_array_weighted(
    k = k, Left = L_uniq, Right = R_uniq, Bk = Bk_marg,
    w = w, z = z, n_grid = n_grid, marginal_iface = TRUE
  )
  expect_equal(gr2$S, gr$S, tolerance = 1e-13, ignore_attr = TRUE)
  expect_equal(gr2$b, gr$b, tolerance = 1e-13)
})

test_that("array mode Poisson parity vs scattered (fixed lambda)", {
  set.seed(21)
  n_grid <- c(9L, 8L)
  axes <- list(
    x1 = seq(0, 1, length.out = n_grid[1]),
    x2 = seq(0, 1, length.out = n_grid[2])
  )
  eta_true <- outer(axes$x1, axes$x2,
                    function(a, b) 0.5 + sin(pi * a) + 0.5 * cos(pi * b))
  Y <- array(rpois(prod(n_grid), exp(eta_true)), dim = n_grid)
  ctrl <- tt_control(max_sweeps = 3L, compute_edf = FALSE, seed = 1L,
                     trace = FALSE, pirls_maxit = 4L)
  idx <- expand.grid(lapply(n_grid, seq_len), KEEP.OUT.ATTRS = FALSE)
  X_sc <- do.call(cbind,
                  lapply(seq_along(n_grid), function(j) axes[[j]][idx[[j]]]))
  y_sc <- as.numeric(Y)
  fit_sc <- ttps(y_sc, X_sc, rank = 2L, k = 6L, lambda = 1,
                 family = "poisson", optimizer = "ALS", backend = "R",
                 control = ctrl)
  fit_ar <- ttps(Y, axes = axes, rank = 2L, k = 6L, lambda = 1,
                 family = "poisson", optimizer = "ALS", array = TRUE,
                 control = ctrl)
  expect_equal(fit_ar$fitted.values, fit_sc$fitted.values,
               tolerance = 1e-9,
               label = "array vs scattered fitted values (Poisson, fixed λ)")
  expect_identical(fit_ar$backend, "R-array")
})

test_that("array mode Poisson parity vs scattered (cGCV sequential)", {
  set.seed(33)
  n_grid <- c(10L, 9L)
  axes <- list(
    x1 = seq(0, 1, length.out = n_grid[1]),
    x2 = seq(0, 1, length.out = n_grid[2])
  )
  eta_true <- outer(axes$x1, axes$x2,
                    function(a, b) 1 + 0.8 * sin(2 * pi * a) * cos(pi * b))
  Y <- array(rpois(prod(n_grid), exp(eta_true)), dim = n_grid)
  ctrl <- tt_control(max_sweeps = 2L, compute_edf = FALSE, seed = 1L,
                     trace = FALSE, pirls_maxit = 3L,
                     cgcv_update = "sequential")
  idx <- expand.grid(lapply(n_grid, seq_len), KEEP.OUT.ATTRS = FALSE)
  X_sc <- do.call(cbind,
                  lapply(seq_along(n_grid), function(j) axes[[j]][idx[[j]]]))
  fit_sc <- ttps(as.numeric(Y), X_sc, rank = 2L, k = 6L, lambda = "cGCV",
                 family = "poisson", optimizer = "ALS", backend = "R",
                 control = ctrl)
  fit_ar <- ttps(Y, axes = axes, rank = 2L, k = 6L, lambda = "cGCV",
                 family = "poisson", optimizer = "ALS", array = TRUE,
                 control = ctrl)
  expect_equal(fit_ar$fitted.values, fit_sc$fitted.values,
               tolerance = 1e-7,
               label = "array vs scattered fitted values (Poisson, cGCV)")
})
