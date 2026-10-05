# Exact EDF (tt_edf_exact): analytic Hessian, gauge projection, Newton
# polishing. References: dense EDF at saturating rank, central differences of
# the training gradient, and finite-difference traces of the polished map.

.ee_ctrl <- function(...) {
  tt_control(max_sweeps = 2000L, pirls_maxit = 200L, tol = 1e-12, seed = 1L,
             compute_edf = FALSE, warn_lambda_boundary = FALSE, ...)
}

test_that("exact EDF equals the dense EDF at saturating rank", {
  set.seed(11)
  axes <- list(x1 = seq(0, 1, length.out = 9), x2 = seq(0, 1, length.out = 8))
  X <- as.matrix(expand.grid(axes))
  f <- sin(2 * pi * X[, 1]) + cos(2 * pi * X[, 2])
  y <- f + rnorm(nrow(X), sd = 0.3)
  yp <- rpois(nrow(X), exp(log(0.8) + f))
  bb <- glam_grid_bases(axes, k = 5L)
  B <- kronecker(bb$B[[2]], bb$B[[1]])
  fg <- ttps(y, X, rank = 5, k = 5, lambda = c(0.5, 2), control = .ee_ctrl())
  dense_g <- sum(diag(solve(crossprod(B) + glam_penalty(c(5, 5), c(0.5, 2)),
                            crossprod(B))))
  expect_equal(tt_edf_exact(fg, y, X)$edf, dense_g, tolerance = 1e-6)
  fp <- ttps(yp, X, family = poisson(), rank = 5, k = 5, lambda = c(1, 1),
             control = .ee_ctrl())
  ep <- tt_edf_exact(fp, yp, X)
  expect_true(ep$converged)
  WB <- crossprod(B, ep$mu * B)
  dense_p <- sum(diag(solve(WB + glam_penalty(c(5, 5), c(1, 1)), WB)))
  expect_equal(ep$edf, dense_p, tolerance = 1e-6)
})

test_that("analytic Hessian matches central differences of the gradient", {
  ns <- asNamespace("TTPsplines")
  ng <- c(5L, 4L, 4L)
  axes <- lapply(ng, function(m) seq(0, 1, length.out = m))
  names(axes) <- paste0("x", 1:3)
  X <- as.matrix(expand.grid(axes))
  set.seed(3)
  Y <- array(sin(2 * pi * X[, 1]) * cos(pi * X[, 2]) + X[, 3]^2 +
               rnorm(nrow(X), sd = 0.2), ng)
  fit <- ttps(Y, axes = axes, array = TRUE, rank = 2, k = 4,
              lambda = rep(0.05, 3),
              control = tt_control(max_sweeps = 20L, seed = 1L,
                                   compute_edf = FALSE))
  basis <- ns$eval_marginal_bases(X, fit$knots, fit$degree, cyclic = fit$cyclic)
  y <- as.numeric(Y)
  lam <- rep(0.05, 3)
  p <- ns$.tt_exact_parts(fit$cores, fit$intercept, basis, y, "gaussian",
                          rep(1, length(y)), rep(0, length(y)), lam, 2L,
                          attr(basis, "cyclic"))
  th <- c(ns$.tt_pack_cores(fit$cores), fit$intercept)
  np <- length(th)
  gr <- function(t) {
    o <- ns$.tt_gaussian_objective(t[-np], y, t[np], basis, fit$cores, NULL, lam)
    c(o$grad, -sum(o$resid))
  }
  Hf <- vapply(seq_len(np), function(j) {
    e <- numeric(np)
    e[j] <- 1e-3
    (gr(th + e) - gr(th - e)) / 2e-3
  }, numeric(np))
  expect_lt(max(abs(p$grad - gr(th))), 1e-10)
  # the Gaussian gradient is quadratic along each coordinate: central
  # differences are exact up to rounding
  expect_lt(max(abs(p$H - (Hf + t(Hf)) / 2)) / max(abs(p$H)), 1e-9)
})

test_that("low-rank exact EDF is the trace of the polished map", {
  ns <- asNamespace("TTPsplines")
  ng <- c(4L, 4L, 3L)
  axes <- lapply(ng, function(m) seq(0, 1, length.out = m))
  names(axes) <- paste0("x", 1:3)
  X <- as.matrix(expand.grid(axes))
  set.seed(5)
  Y <- array(rpois(nrow(X), exp(0.5 + sin(2 * pi * X[, 1]) + X[, 3])), ng)
  fit <- ttps(Y, axes = axes, array = TRUE, family = poisson(), rank = 2, k = 4,
              lambda = rep(0.1, 3), control = .ee_ctrl())
  e <- tt_edf_exact(fit, Y = Y, axes = axes)
  expect_true(e$converged)
  expect_identical(e$n_negative, 0L)
  basis <- ns$eval_marginal_bases(X, fit$knots, fit$degree, cyclic = fit$cyclic)
  y <- as.numeric(Y)
  h <- 1e-4
  tr <- 0
  for (i in seq_along(y)) {
    mu <- vapply(c(h, -h), function(s) {
      yy <- y
      yy[i] <- yy[i] + s
      ns$.tt_exact_edf_core(e$cores, e$intercept, basis, yy, poisson(),
                            offset = fit$offset, lambda = fit$lambda)$mu[i]
    }, numeric(1))
    tr <- tr + (mu[1] - mu[2]) / (2 * h)
  }
  expect_equal(e$edf, tr, tolerance = 1e-4)
  # without polishing the trace is taken off the stationary point
  e0 <- tt_edf_exact(fit, Y = Y, axes = axes, polish = FALSE)
  expect_true(is.finite(e0$edf))
})
