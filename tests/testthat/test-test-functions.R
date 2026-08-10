test_that("test-function truth evaluators match known properties", {
  # Ishigami at origin-ish: sin(0)+a sin(0)^2 + ... = 0
  expect_equal(f_ishigami(c(0, 0, 0)), 0, tolerance = 1e-12)
  # Sobol g with a=0 reduces first factor to |4x-2|
  expect_equal(f_sobol_g(matrix(0.5, 1, 1), a = 0), 0, tolerance = 1e-12)
  # Friedman at (0.5,0.5,0.5,0,0): 10 sin(pi/4) + 0 + 0 + 0
  expect_equal(
    f_friedman(matrix(c(0.5, 0.5, 0.5, 0, 0), 1)),
    10 * sin(pi * 0.25),
    tolerance = 1e-10
  )
})

test_that("simulators return coherent X/y/f and as.data.frame works", {
  dat <- simulate_ishigami(n = 50, sigma = 0.1, seed = 7)
  expect_s3_class(dat, "ttps_sim")
  expect_equal(dim(dat$X), c(50L, 3L))
  expect_equal(length(dat$y), 50L)
  expect_equal(length(dat$f), 50L)
  df <- as.data.frame(dat)
  expect_equal(nrow(df), 50L)
  expect_true(all(c("x1", "x2", "x3", "y", "f") %in% names(df)))
})

test_that("packaged datasets load and fit with ttpspline", {
  rda <- if (file.exists("../../data/ishigami.rda")) {
    "../../data/ishigami.rda"
  } else if (file.exists("data/ishigami.rda")) {
    "data/ishigami.rda"
  } else {
    skip("packaged ishigami.rda not found")
  }
  e <- new.env(parent = emptyenv())
  load(rda, envir = e)
  X <- as.matrix(e$ishigami[, c("x1", "x2", "x3")])
  fit <- ttps(
    e$ishigami$y, X, rank = 2, k = 5, lambda = 1,
    control = tt_control(max_sweeps = 4, backend = "R", compute_edf = FALSE)
  )
  expect_equal(fit$d, 3L)
  expect_true(is.finite(fit$deviance))
})

test_that("sobol_g and friedman simulators accept GLM families", {
  sp <- simulate_sobol_g(n = 40, d = 3, seed = 1, family = "poisson")
  expect_equal(sp$family, "poisson")
  expect_true(all(sp$y >= 0))
  sb <- simulate_friedman(n = 40, seed = 2, family = "binomial")
  expect_equal(sb$family, "binomial")
  expect_true(all(sb$y %in% c(0L, 1L)))
})
