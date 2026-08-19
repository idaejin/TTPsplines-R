## Parity tests: ttps(array = TRUE) vs ttps(array = FALSE) on a complete grid.
## Criterion: fitted values must be identical (up to machine precision).

test_that("array mode parity 3D, rank=2, fixed lambda", {
  set.seed(42)
  n_grid <- c(8L, 7L, 6L)
  axes <- list(
    x1 = seq(0, 1, length.out = n_grid[1]),
    x2 = seq(0, 1, length.out = n_grid[2]),
    x3 = seq(0, 1, length.out = n_grid[3])
  )
  truth <- outer(
    outer(sin(2 * pi * axes$x1), cos(2 * pi * axes$x2), "+"),
    axes$x3,
    function(a, b) a + 0.5 * b
  )
  set.seed(7)
  Y <- truth + array(rnorm(prod(n_grid), 0, 0.2), dim(truth))
  ctrl <- tt_control(max_sweeps = 15L, compute_edf = FALSE, seed = 1L,
                     trace = FALSE)
  # Scattered (array = FALSE)
  idx <- expand.grid(lapply(n_grid, seq_len), KEEP.OUT.ATTRS = FALSE)
  X_sc <- do.call(cbind,
                  lapply(seq_along(n_grid), function(j) axes[[j]][idx[[j]]]))
  y_sc <- as.numeric(Y)
  fit_sc  <- ttps(y_sc, X_sc, rank = 2L, k = 6L, lambda = 1,
                  optimizer = "ALS", control = ctrl)
  # Array (array = TRUE)
  fit_arr <- ttps(Y, axes = axes, rank = 2L, k = 6L, lambda = 1,
                  optimizer = "ALS", array = TRUE, control = ctrl)
  expect_equal(fit_arr$fitted.values, fit_sc$fitted.values,
               tolerance = 1e-12,
               label = "array vs scattered fitted values (3D, rank=2, fixed λ)")
})

test_that("array mode parity 2D, rank=3, cGCV lambda (sequential)", {
  set.seed(123)
  n_grid <- c(12L, 10L)
  axes <- list(
    x1 = seq(0, 1, length.out = n_grid[1]),
    x2 = seq(0, 1, length.out = n_grid[2])
  )
  Y <- outer(axes$x1, axes$x2, function(a, b) sin(pi * a) * cos(pi * b))
  Y <- Y + matrix(rnorm(prod(n_grid), 0, 0.1), n_grid[1], n_grid[2])
  # Force sequential so both paths use the same update mode
  ctrl <- tt_control(max_sweeps = 10L, compute_edf = FALSE, seed = 1L,
                     trace = FALSE, cgcv_update = "sequential")
  idx <- expand.grid(lapply(n_grid, seq_len), KEEP.OUT.ATTRS = FALSE)
  X_sc <- do.call(cbind,
                  lapply(seq_along(n_grid), function(j) axes[[j]][idx[[j]]]))
  y_sc <- as.numeric(Y)
  fit_sc  <- ttps(y_sc, X_sc, rank = 3L, k = 8L, lambda = "cGCV",
                  optimizer = "ALS", control = ctrl)
  fit_arr <- ttps(Y, axes = axes, rank = 3L, k = 8L, lambda = "cGCV",
                  optimizer = "ALS", array = TRUE, control = ctrl)
  expect_equal(fit_arr$fitted.values, fit_sc$fitted.values,
               tolerance = 1e-12,
               label = "array vs scattered fitted values (2D, rank=3, cGCV)")
})

test_that("array mode: default axes (unit interval) produces valid fit", {
  set.seed(9)
  n_grid <- c(6L, 5L)
  Y <- matrix(rnorm(prod(n_grid)), n_grid[1], n_grid[2])
  ctrl <- tt_control(max_sweeps = 5L, compute_edf = FALSE, seed = 1L,
                     trace = FALSE)
  fit <- ttps(Y, rank = 2L, k = 4L, lambda = 1, optimizer = "ALS",
              array = TRUE, control = ctrl)
  expect_length(fit$fitted.values, prod(n_grid))
  expect_false(anyNA(fit$fitted.values))
})

test_that("array mode: input validation errors", {
  expect_error(ttps(1:10, rank = 2L, k = 4L, lambda = 1, array = TRUE),
               regexp = "d-way array")
  Y <- matrix(1:12, 3, 4)
  expect_error(ttps(Y, axes = list(1:3), rank = 2L, k = 4L, lambda = 1,
                    array = TRUE),
               regexp = "length d")
  expect_error(ttps(Y, axes = list(1:5, 1:4), rank = 2L, k = 4L, lambda = 1,
                    array = TRUE),
               regexp = "axes\\[\\[1\\]\\]")
})
