# gGCV for array data (K-ALS / TT-GLAM): GDF against exact dense GLAM EDF at
# saturating rank (d = 2, r = K), selection against exact dense GLAM GCV.

.ga_setup <- function() {
  ng <- c(9L, 8L)
  axes <- list(x1 = seq(0, 1, length.out = ng[1]), x2 = seq(0, 1, length.out = ng[2]))
  K <- 5L
  bb <- glam_grid_bases(axes, k = K)
  list(ng = ng, axes = axes, K = K, B = kronecker(bb$B[[2]], bb$B[[1]]),
       ctrl = tt_control(max_sweeps = 30L, pirls_maxit = 30L, seed = 1L,
                         compute_edf = FALSE))
}

.ga_dense_poisson <- function(y, B, P, wo) {
  e <- rep(log(mean(y[wo > 0])), length(y))
  for (it in 1:200) {
    mu <- exp(e)
    z <- e + (y - mu) / mu
    W <- wo * mu
    th <- solve(crossprod(B, W * B) + P, crossprod(B, W * z))
    e <- as.numeric(B %*% th)
  }
  exp(e)
}

test_that("Gaussian array GDF equals exact GLAM EDF at saturating rank", {
  s <- .ga_setup()
  set.seed(11)
  f <- outer(s$axes$x1, s$axes$x2, function(a, b) sin(2 * pi * a) + cos(2 * pi * b))
  Y <- array(f + rnorm(prod(s$ng), sd = 0.3), s$ng)
  lam <- c(0.5, 2)
  P <- glam_penalty(c(s$K, s$K), lam)
  edf <- sum(diag(solve(crossprod(s$B) + P, crossprod(s$B))))
  gu <- tt_gdf_array(Y, lambda = lam, axes = s$axes, rank = s$K, k = s$K,
                     probes = "unit", control = s$ctrl)
  expect_equal(gu$gdf, edf, tolerance = 1e-3)
  gh <- tt_gdf_array(Y, lambda = lam, axes = s$axes, rank = s$K, k = s$K,
                     M = 40L, control = s$ctrl)
  expect_lt(abs(gh$gdf - edf), 4 * gh$gdf_se)
  mu_d <- s$B %*% solve(crossprod(s$B) + P, crossprod(s$B, as.numeric(Y)))
  expect_equal(as.numeric(gu$fit$fitted.values), as.numeric(mu_d), tolerance = 1e-6)
})

test_that("weighted zero-heavy Poisson array GDF equals exact EDF (no clipping)", {
  s <- .ga_setup()
  set.seed(12)
  eta <- outer(s$axes$x1, s$axes$x2,
               function(a, b) log(0.2) + sin(2 * pi * a) + cos(2 * pi * b))
  Y <- array(rpois(prod(s$ng), exp(eta)), s$ng)
  expect_gt(mean(Y == 0), 0.5)
  Wt <- array(1, s$ng)
  Wt[1:3, 1:3] <- 0
  lam <- c(1, 1)
  P <- glam_penalty(c(s$K, s$K), lam)
  wo <- as.numeric(Wt)
  mu_d <- .ga_dense_poisson(as.numeric(Y), s$B, P, wo)
  WB <- crossprod(s$B, wo * mu_d * s$B)
  edf <- sum(diag(solve(WB + P, WB)))
  pu <- tt_gdf_array(Y, lambda = lam, axes = s$axes, family = poisson(),
                     rank = s$K, k = s$K, weights = Wt, probes = "unit",
                     control = s$ctrl)
  expect_equal(pu$n_eff, sum(wo > 0))
  expect_equal(as.numeric(pu$fit$fitted.values), mu_d, tolerance = 1e-6)
  expect_equal(pu$gdf, edf, tolerance = 1e-3)
})

test_that("isotropic gGCV-array lands near the exact dense GLAM GCV argmin", {
  s <- .ga_setup()
  set.seed(13)
  f <- outer(s$axes$x1, s$axes$x2, function(a, b) sin(2 * pi * a) + cos(2 * pi * b))
  Y <- array(f + rnorm(prod(s$ng), sd = 0.4), s$ng)
  y <- as.numeric(Y)
  n <- length(y)
  grid <- seq(-3, 3, by = 0.05)
  gcv <- vapply(grid, function(t) {
    P <- glam_penalty(c(s$K, s$K), rep(10^t, 2))
    M <- crossprod(s$B) + P
    mu <- s$B %*% solve(M, crossprod(s$B, y))
    n * sum((y - mu)^2) / (n - sum(diag(solve(M, crossprod(s$B)))))^2
  }, numeric(1))
  sel <- tt_ggcv_array(Y, axes = s$axes, rank = s$K, k = s$K,
                       theta_lower = -3, theta_upper = 3,
                       M_search = 8L, M_final = 24L, control = s$ctrl)
  expect_lt(abs(sel$theta - grid[which.min(gcv)]), 0.3)
  expect_true(is.finite(sel$gcv))
  expect_s3_class(sel$fit, "ttpspline")
})

test_that("ttps(array = TRUE, lambda = 'gGCV') dispatches to the array selector", {
  s <- .ga_setup()
  set.seed(14)
  Y <- array(rpois(prod(s$ng), 2), s$ng)
  cc <- s$ctrl
  cc$ggcv_array_M_search <- 2L
  cc$ggcv_array_M_final <- 4L
  cc$ggcv_array_n_refine <- 0L
  cc$ggcv_array_groups <- c(1L, 2L)
  cc$ggcv_array_n_global <- 4L
  cc$warn_lambda_boundary <- FALSE
  set.seed(99); u1 <- runif(1)
  set.seed(99)
  fit <- ttps(Y, axes = s$axes, array = TRUE, family = poisson(), rank = 2L,
              k = s$K, lambda = "gGCV", control = cc)
  u2 <- runif(1)
  expect_s3_class(fit, "ttpspline")
  expect_identical(fit$lambda_method, "gGCV")
  expect_identical(fit$ggcv$mode, "array")
  expect_length(fit$lambda, 2L)
  expect_true(is.finite(fit$ggcv$gdf))
  expect_identical(u1, u2)
})

test_that("ubre needs a known scale for Gaussian", {
  s <- .ga_setup()
  Y <- array(rnorm(prod(s$ng)), s$ng)
  expect_error(
    tt_ggcv_array(Y, axes = s$axes, rank = 2L, k = s$K, criterion = "ubre",
                  control = s$ctrl),
    "known `scale`"
  )
})
