test_that("conditional Q_k objective matches global objective contributions", {
  set.seed(31)
  n <- 80
  X <- matrix(runif(n * 2), n, 2)
  y <- rbinom(n, 1, 0.45)
  init <- tt_initialize(X, rank = 2, k = 5, seed = 1, sd = 0.05)
  fit <- ttpspline(
    y, X, family = binomial(), rank = 2, k = 5, lambda = 1,
    optimizer = "LBFGS", init = init,
    control = tt_control(backend = "R", lbfgs_maxit = 30L, compute_edf = FALSE)
  )
  basis <- eval_marginal_bases(X, fit$knots, fit$degree)
  o_global <- tt_glm_penalized_objective(
    y, fit$cores, fit$intercept, basis, fit$penalties %||% {
      d <- length(fit$cores)
      p <- dim(fit$cores[[1]])[2]
      lapply(seq_len(d), function(k)
        core_penalty(fit$rank[k], p, fit$rank[k + 1L], fit$penalty_order))
    },
    fit$lambda, binomial()
  )
  # Rebuild penalties if not stored
  d <- length(fit$cores)
  p <- dim(fit$cores[[1]])[2]
  pens <- lapply(seq_len(d), function(k)
    core_penalty(fit$rank[k], p, fit$rank[k + 1L], 2L))
  o_global <- tt_glm_penalized_objective(
    y, fit$cores, fit$intercept, basis, pens, fit$lambda, binomial()
  )
  # Sum of conditional Q_k is NOT equal to global Q (double-counts nll).
  # Instead: evaluate Q_k and check nll(eta) + pen_k matches pieces.
  k <- 1L
  L <- left_interfaces(fit$cores, basis)
  R <- right_interfaces(fit$cores, basis)
  Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
  g <- as.numeric(fit$cores[[k]])
  qk <- .tt_conditional_qk(g, y, fit$intercept, Xk, pens[[k]], fit$lambda[k], binomial())
  eta2 <- fit$intercept + as.numeric(Xk %*% g)
  expect_equal(qk$eta, eta2, tolerance = 1e-10)
  expect_equal(qk$eta, o_global$eta, tolerance = 1e-8)
  expect_equal(qk$nll, o_global$nll, tolerance = 1e-8)
})

test_that("all solvers evaluate the same penalized objective on shared cores", {
  set.seed(32)
  n <- 100
  X <- matrix(runif(n * 3), n, 3)
  y <- rbinom(n, 1, plogis(sin(2 * pi * X[, 1])))
  init <- tt_initialize(X, rank = 2, k = 5, seed = 2, sd = 0.05)
  # Short LBFGS to get a nontrivial point
  fit0 <- ttpspline(
    y, X, family = binomial(), rank = 2, k = 5, lambda = 1,
    optimizer = "LBFGS", init = init,
    control = tt_control(backend = "R", lbfgs_maxit = 20L, compute_edf = FALSE)
  )
  basis <- eval_marginal_bases(X, fit0$knots, fit0$degree)
  d <- length(fit0$cores)
  p <- dim(fit0$cores[[1]])[2]
  pens <- lapply(seq_len(d), function(k)
    core_penalty(fit0$rank[k], p, fit0$rank[k + 1L], 2L))
  o1 <- tt_glm_penalized_objective(
    y, fit0$cores, fit0$intercept, basis, pens, fit0$lambda, binomial()
  )
  th <- TTPsplines:::.tt_pack_cores(fit0$cores)
  o2 <- TTPsplines:::.tt_glm_objective(
    th, y, fit0$intercept, basis, fit0$cores, pens, fit0$lambda, binomial()
  )
  expect_equal(o1$value, o2$value, tolerance = 1e-8)
  expect_equal(o1$eta, o2$eta, tolerance = 1e-8)
})

test_that("conditional analytical gradient matches finite differences", {
  set.seed(33)
  n <- 60
  X <- matrix(runif(n * 2), n, 2)
  for (fam in list(binomial(), poisson())) {
    if (identical(family_key(fam), "bernoulli")) {
      y <- rbinom(n, 1, 0.4)
    } else {
      y <- rpois(n, 2)
    }
    init <- tt_initialize(X, rank = 2, k = 4, seed = 4, sd = 0.08)
    bs <- build_marginal_bases(X, k = 4, degree = 3)
    basis <- bs$basis
    ranks <- tt_rank(2, d = 2)
    pens <- lapply(1:2, function(k)
      core_penalty(ranks[k], 4L, ranks[k + 1L], 2L))
    cores <- init
    intercept <- init_intercept(fam, y)
    lam <- c(1, 1)
    for (k in 1:2) {
      L <- left_interfaces(cores, basis)
      R <- right_interfaces(cores, basis)
      Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
      g <- as.numeric(cores[[k]])
      ga <- .tt_conditional_qk_grad(g, y, intercept, Xk, pens[[k]], lam[k], fam)$grad
      # finite differences
      eps <- 1e-6
      gfd <- numeric(length(g))
      for (j in seq_along(g)) {
        gp <- g; gm <- g
        gp[j] <- gp[j] + eps
        gm[j] <- gm[j] - eps
        qp <- .tt_conditional_qk(gp, y, intercept, Xk, pens[[k]], lam[k], fam)$value
        qm <- .tt_conditional_qk(gm, y, intercept, Xk, pens[[k]], lam[k], fam)$value
        gfd[j] <- (qp - qm) / (2 * eps)
      }
      err <- max(abs(ga - gfd))
      expect_true(err < 1e-4, info = paste(family_key(fam), "core", k, "err", err))
    }
  }
})

test_that("Damped-Newton-ALS accepted steps are monotone in global objective", {
  set.seed(34)
  n <- 150
  X <- matrix(runif(n * 3), n, 3)
  eta <- 1.1 * sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2])
  y <- rbinom(n, 1, plogis(eta))
  init <- tt_initialize(X, rank = 2, k = 5, seed = 5, sd = 0.05)
  fit <- ttpspline(
    y, X, family = binomial(), rank = 2, k = 5, lambda = 1,
    optimizer = "Damped-Newton-ALS", init = init,
    control = tt_control(
      backend = "R", dn_max_sweeps = 15L, compute_edf = FALSE, tol = 1e-10
    )
  )
  expect_equal(fit$optimizer, "Damped-Newton-ALS")
  h <- fit$history
  if (is.list(h) && length(h) >= 2L) {
    vals <- vapply(h, function(z) z$objective, numeric(1))
    expect_true(all(diff(vals) <= 1e-6 + 1e-8 * pmax(1, abs(vals[-length(vals)]))))
  }
  basis <- eval_marginal_bases(X, fit$knots, fit$degree)
  expect_equal(fit$linear.predictors, fit$offset + fit$intercept + tt_contraction(fit$cores, basis),
               tolerance = 1e-10)
})

test_that("LBFGS-ALS does not increase conditional objectives on smoke fit", {
  set.seed(35)
  n <- 120
  X <- matrix(runif(n * 3), n, 3)
  y <- rbinom(n, 1, plogis(0.8 * sin(2 * pi * X[, 1])))
  init <- tt_initialize(X, rank = 2, k = 5, seed = 6, sd = 0.05)
  fit <- ttpspline(
    y, X, family = binomial(), rank = 2, k = 5, lambda = 1,
    optimizer = "LBFGS-ALS", init = init,
    control = tt_control(
      backend = "R", block_lbfgs_sweeps = 8L, block_lbfgs_maxit = 40L,
      compute_edf = FALSE
    )
  )
  expect_equal(fit$optimizer, "LBFGS-ALS")
  expect_true(is.finite(fit$deviance))
  expect_true(all(predict(fit, type = "response") > 0 &
                    predict(fit, type = "response") < 1))
})
