test_that("Bernoulli penalized objective matches LBFGS internal evaluator", {
  set.seed(9)
  n <- 120
  X <- matrix(runif(n * 3), n, 3)
  eta <- 1.1 * sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2])
  y <- rbinom(n, 1, plogis(eta))
  init <- tt_initialize(X, rank = 2, k = 5, seed = 1, sd = 0.05)
  fit <- ttps(
    y, X, family = binomial(), rank = 2, k = 5, lambda = 1,
    optimizer = "ALS", init = init,
    control = tt_control(backend = "R", pirls_maxit = 2L,
                         als_sweeps_per_pirls = 1L, compute_edf = FALSE,
                         damping = FALSE)
  )
  o1 <- tt_objective(fit, X, y)
  basis <- eval_marginal_bases(X, fit$knots, fit$degree)
  d <- 3L
  p <- ncol(basis[[1]])
  pens <- lapply(seq_len(d), function(k) {
    core_penalty(fit$rank[k], p, fit$rank[k + 1L], 2L)
  })
  th <- TTPsplines:::.tt_pack_cores(fit$cores)
  o2 <- TTPsplines:::.tt_glm_objective(
    th, y, fit$intercept, basis, fit$cores, pens, fit$lambda, binomial()
  )
  expect_equal(o1$value, o2$value, tolerance = 1e-10)
  expect_equal(o1$eta, o2$eta, tolerance = 1e-10)
})
