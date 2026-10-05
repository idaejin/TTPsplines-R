# Unified gGCV engine (R/ggcv_engine.R) on its scattered and array maps:
# determinism across workers, parity of scattered rows and the array, GDF
# against the exact dense EDF at saturating rank (d = 2, rank = k, where the
# TT fit is the dense tensor-product P-spline fit), Berman-Turner-like Poisson
# rows with offsets (UBRE, non-negative probes, dense UBRE argmin), the paired
# Monte Carlo SE, parsimony of the grouped search (with GDFs checked against
# the dense smoother), cache reuse, convergence flags, ttps(lambda = "gGCV")
# with its refit, the deprecated working mode (tt_ggcv_poisson(), the array
# routes), the noise-aware final choice (.ggcv_choose() on synthetic
# candidates), the decision-aware budget check (the first winner kept when
# no converged, stable candidate can replace it; the warning texts), the
# adaptive iteration budget (the pilot's objective rule, B and the cap),
# tt_gdf() against tt_gdf_array() and the dense EDF, probe validation (the
# GDF-dependent CV limit) and unstable evaluations, the numerics (the
# minimum-norm global-mode core solve, probe seeds beyond the integer range,
# UBRE on the AIC scale, the Gaussian probe step), Poisson smoothing groups
# in every entry point, and fit errors surfacing with their message.

# Tiny selection budgets leave the fixed-lambda fits unconverged on purpose.
# The convergence warning has its own test; tests that only need a selection
# muffle that warning and nothing else.
.ge_quiet_budget <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("not converged", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# Value of `expr` and the messages of all warnings it raised.
.ge_warnings <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

# 9 x 8 grid on the unit square; rows in expand.grid order = as.numeric(Y).
.ge_grid <- function() {
  axes <- list(x1 = seq(0, 1, length.out = 9L), x2 = seq(0, 1, length.out = 8L))
  list(ng = c(9L, 8L), axes = axes, X = as.matrix(expand.grid(axes)))
}

# Scattered map and selector state as tt_ggcv() builds them (cold probes,
# probe seed 1, adaptive budget unless `budget = "fixed"`).
.ge_state <- function(y, X, family = stats::gaussian(), rank, k, offset = NULL,
                      groups = NULL, control, budget = "adaptive") {
  map <- .ggcv_map_scattered(y, X, family, rank = rank, k = k, degree = 3L,
                             penalty_order = 2L, cyclic = NULL, period = NULL,
                             knots = NULL, offset = offset)
  .ggcv_setup(map, groups = groups, criterion = "auto", scale = NULL,
              probe_init = "cold", probe_budget = NULL, epsilon_rel = 1e-3,
              n_cores = 1L, seed = 1L, control = control, budget = budget)
}

# Dense tensor-product design of scattered rows (x1 index fastest, the
# coefficient order of glam_penalty()).
.ge_row_tensor <- function(B1, B2) {
  K1 <- ncol(B1)
  K2 <- ncol(B2)
  B1[, rep(seq_len(K1), times = K2), drop = FALSE] *
    B2[, rep(seq_len(K2), each = K1), drop = FALSE]
}

# Exact penalized Poisson fit with an offset (dense PIRLS) and its EDF.
.ge_dense_poisson <- function(y, B, P, offset) {
  eta <- offset + log(sum(y) / sum(exp(offset)))
  for (it in 1:200) {
    mu <- exp(eta)
    z <- eta - offset + (y - mu) / mu
    beta <- solve(crossprod(B, mu * B) + P, crossprod(B, mu * z))
    eta_new <- offset + as.numeric(B %*% beta)
    done <- max(abs(eta_new - eta)) < 1e-12
    eta <- eta_new
    if (done) break
  }
  mu <- exp(eta)
  WB <- crossprod(B, mu * B)
  list(mu = mu, edf = sum(diag(solve(WB + P, WB))))
}

# Berman-Turner-like rows as in the lab's point-process runs: the occupied
# cells of an m x m lattice (counts, offset log(cell area)) plus half of the
# empty cells as y = 0 rows whose offset log(area * n_empty / n_sampled) lets
# them stand for all empty cells.
.ge_bt_rows <- function(seed, m = 16L, n_events = 400, amp = 1.5) {
  set.seed(seed)
  cen <- (seq_len(m) - 0.5) / m
  Xc <- as.matrix(expand.grid(x1 = cen, x2 = cen))
  area <- 1 / m^2
  rate <- exp(amp * (sin(2 * pi * Xc[, 1]) + cos(2 * pi * Xc[, 2])))
  yc <- stats::rpois(nrow(Xc), n_events * area * rate / mean(rate))
  occ <- which(yc > 0)
  emp <- which(yc == 0)
  n_q <- ceiling(length(emp) / 2)
  q <- emp[sample.int(length(emp), n_q)]
  list(X = Xc[c(occ, q), , drop = FALSE], y = c(yc[occ], numeric(n_q)),
       offset = c(rep(log(area), length(occ)),
                  rep(log(area * length(emp) / n_q), n_q)))
}

# Gaussian rows on a random design, y = truth(X) + N(0, 0.3^2), with an
# isotropic truth (both margins equally wiggly) and an anisotropic one (x1
# wiggly, x2 linear, i.e. in the null space of the second-order penalty).
.ge_gauss_xy <- function(seed, truth, n = 150) {
  set.seed(seed)
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  list(X = X, y = truth(X) + stats::rnorm(n, sd = 0.3))
}
.ge_truth_iso <- function(X) sin(2 * pi * X[, 1]) + sin(2 * pi * X[, 2])
.ge_truth_aniso <- function(X) sin(2 * pi * X[, 1]) + X[, 2]

# Scattered Gaussian and Poisson responses (with an exposure offset) on the
# same random design.
.ge_scattered <- function() {
  set.seed(1)
  n <- 120
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  y_gauss <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2] + stats::rnorm(n, sd = 0.3)
  y_pois <- stats::rpois(n, exp(0.5 + sin(2 * pi * X[, 1])))
  list(X = X, y_gauss = y_gauss, y_pois = y_pois,
       offset = log(stats::runif(n, 0.5, 1.5)))
}

# Dense tensor-product smoother S(lambda) of scattered rows (d = 2, rank = k,
# where the converged TT fit is the dense P-spline fit).
.ge_dense_smoother <- function(X, K) {
  bs <- build_marginal_bases(X, k = K, degree = 3L)
  B <- .ge_row_tensor(bs$basis[[1L]], bs$basis[[2L]])
  function(lambda) B %*% solve(crossprod(B) + glam_penalty(c(K, K), lambda), t(B))
}

# Seeded probes 1..M of the engine on n rows (probe seed 1).
.ge_probes <- function(n, M) {
  st <- list(seed = 1L, n_rows = n, support = rep(TRUE, n), n_eff = n)
  vapply(seq_len(M), .ggcv_probe, numeric(n), st = st)
}

# Fits of one adaptive evaluation with M probes (tt_gdf() value `ev`): the
# pilot's cold runs (b0, 2 b0, ... up to the cap, until a run reaches the
# iteration where the objective settled), the estimator at B unless the last
# pilot run is it, the reference and the M probes; M + 2 more per
# probe-validation round.
.ge_adaptive_fits <- function(ev, M, b0, cap) {
  b <- min(b0, cap)
  runs <- 1L
  while (b < cap && (!ev$converged || b < ev$n_iter)) {
    b <- min(2L * b, cap)
    runs <- runs + 1L
  }
  as.integer(runs + (b != ev$budget) + 1L + M + ev$reevals * (M + 2L))
}

# Mock of .ggcv_eval_batch() (`eval_batch`, the original) that scales the
# probe contributions of the k-th budget check by f[k] and re-scores it: a
# made-up GDF change at twice the budget.
.ge_scaled_checks <- function(f, eval_batch) {
  force(eval_batch)
  n_budget <- 0L
  function(thetas, M, st, stage, budget_mult = 1L, keep_fit = FALSE) {
    recs <- eval_batch(thetas, M, st, stage, budget_mult, keep_fit)
    if (identical(stage, "budget")) {
      n_budget <<- n_budget + 1L
      if (n_budget <= length(f)) {
        r <- recs[[1L]]
        sc <- .ggcv_score(r$deviance, f[n_budget] * r$contrib, M, st)
        recs[[1L]][names(sc)] <- sc
      }
    }
    recs
  }
}

# Mock of .ggcv_eval_batch() (`eval_batch`, the original): every final
# candidate but the first (the isotropic optimum) is flagged unconverged, and
# the probe contributions of every budget check are multiplied by `f` (`NA`:
# no usable GDF at twice the budget) and re-scored.
.ge_only_first_converged <- function(f, eval_batch) {
  force(eval_batch)
  function(thetas, M, st, stage, budget_mult = 1L, keep_fit = FALSE) {
    recs <- eval_batch(thetas, M, st, stage, budget_mult, keep_fit)
    if (identical(stage, "final")) {
      for (q in seq_along(recs)[-1L]) recs[[q]]$converged <- FALSE
    }
    if (identical(stage, "budget")) {
      r <- recs[[1L]]
      sc <- .ggcv_score(r$deviance, f * r$contrib, M, st)
      recs[[1L]][names(sc)] <- sc
    }
    recs
  }
}

# Mock of .ggcv_eval_batch() that flags the records whose theta satisfies
# `pred` as unstable.
.ge_mark_unstable <- function(pred, eval_batch) {
  force(eval_batch)
  function(thetas, M, st, stage, budget_mult = 1L, keep_fit = FALSE) {
    lapply(eval_batch(thetas, M, st, stage, budget_mult, keep_fit), function(r) {
      if (pred(r$theta)) r$stable <- FALSE
      r
    })
  }
}

test_that("one and two workers give the identical selection", {
  skip_on_cran()
  skip_on_os("windows")
  dat <- .ge_gauss_xy(61L, .ge_truth_aniso)
  run <- function(n_cores) {
    tt_ggcv(dat$y, dat$X, rank = 6L, k = 6L, groups = 1:2, theta_lower = -3,
            theta_upper = 3, n_grid = 5L, n_refine = 1L, M_search = 3L,
            M_final = 6L, n_cores = n_cores,
            control = tt_control(max_sweeps = 20L, seed = 1L, compute_edf = FALSE))
  }
  s1 <- run(1L)
  s2 <- run(2L)
  expect_identical(s2$lambda, s1$lambda)
  expect_identical(s2$score_dev, s1$score_dev)
  expect_identical(s2$gdf, s1$gdf)
  # the whole search path (grouped pattern moves included), wall times aside
  expect_true(any(s1$search$stage == "pattern"))
  keep <- setdiff(names(s1$search), "time_s")
  expect_identical(s2$search[keep], s1$search[keep])
  expect_identical(s2$paired, s1$paired)
})

test_that("scattered rows and the array give the same fixed-lambda map and GDF", {
  g <- .ge_grid()
  ctrl <- tt_control(max_sweeps = 10L, pirls_maxit = 10L, seed = 1L,
                     compute_edf = FALSE)
  lam <- c(0.5, 2)
  eta <- sin(2 * pi * g$X[, 1]) + cos(2 * pi * g$X[, 2])
  set.seed(21)
  n <- nrow(g$X)
  cases <- list(
    list(family = stats::gaussian(), y = eta + stats::rnorm(n, sd = 0.3)),
    list(family = stats::poisson(), y = stats::rpois(n, exp(log(0.2) + eta)))
  )
  expect_gt(mean(cases[[2L]]$y == 0), 0.5)
  for (cs in cases) {
    arr <- tt_gdf_array(array(cs$y, g$ng), lambda = lam, axes = g$axes,
                        family = cs$family, rank = 3L, k = 6L, M = 8L,
                        control = ctrl)
    st <- .ge_state(cs$y, g$X, cs$family, rank = 3L, k = 6L, groups = 1:2,
                    control = ctrl)
    sca <- .ggcv_eval_batch(matrix(log10(lam), nrow = 1L), 8L, st, "gdf",
                            keep_fit = TRUE)[[1L]]
    # same step and probes, so the per-probe contributions must agree too
    expect_equal(st$eps, arr$epsilon)
    expect_equal(as.numeric(sca$fit$fitted.values),
                 as.numeric(arr$fit$fitted.values), tolerance = 1e-6)
    expect_equal(sca$contrib, arr$contrib, tolerance = 1e-6)
    expect_equal(sca$gdf, arr$gdf, tolerance = 1e-6)
    expect_equal(sca$deviance, arr$deviance, tolerance = 1e-6)
  }
})

test_that("scattered Gaussian GDF equals the exact dense EDF at saturating rank", {
  set.seed(31)
  n <- 80
  K <- 5L
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  y <- sin(2 * pi * X[, 1]) + cos(2 * pi * X[, 2]) + stats::rnorm(n, sd = 0.3)
  lam <- c(0.5, 2)
  bs <- build_marginal_bases(X, k = K, degree = 3L)
  B <- .ge_row_tensor(bs$basis[[1L]], bs$basis[[2L]])
  A <- crossprod(B) + glam_penalty(c(K, K), lam)
  edf <- sum(diag(solve(A, crossprod(B))))
  mu <- as.numeric(B %*% solve(A, crossprod(B, y)))
  st <- .ge_state(y, X, rank = K, k = K, groups = 1:2,
                  control = tt_control(max_sweeps = 20L, seed = 1L,
                                       compute_edf = FALSE))
  unit <- .ggcv_eval_unit(log10(lam), st)
  expect_true(unit$converged)
  expect_equal(as.numeric(unit$fit$fitted.values), mu, tolerance = 1e-6)
  expect_equal(unit$gdf, edf, tolerance = 1e-6)
  mc <- .ggcv_eval_batch(matrix(log10(lam), nrow = 1L), 40L, st, "gdf")[[1L]]
  expect_lt(abs(mc$gdf - edf), 4 * mc$gdf_se)
})

test_that("Berman-Turner Poisson rows: UBRE, non-negative probes, dense argmin", {
  skip_if_not_installed("testthat", "3.1.7")  # local_mocked_bindings()
  bt <- .ge_bt_rows(42L)
  n <- length(bt$y)
  K <- 6L
  zero <- bt$y == 0
  expect_gt(sum(zero), 0)
  expect_false(any(bt$offset[zero] %in% bt$offset[!zero]))
  # record every response the selector fits through the scattered map
  seen <- new.env(parent = emptyenv())
  seen$y <- list()
  map_scattered <- .ggcv_map_scattered
  local_mocked_bindings(.ggcv_map_scattered = function(...) {
    map <- map_scattered(...)
    fit <- map$fit
    map$fit <- function(yvec, lambda, ctrl, init = NULL) {
      seen$y[[length(seen$y) + 1L]] <- yvec
      fit(yvec, lambda, ctrl, init)
    }
    map
  })
  sel <- tt_ggcv(bt$y, bt$X, rank = K, k = K, family = stats::poisson(),
                 offset = bt$offset, theta_lower = -3, theta_upper = 3,
                 control = tt_control(pirls_maxit = 15L, seed = 1L,
                                      compute_edf = FALSE))
  expect_identical(sel$criterion, "ubre")
  expect_equal(sel$score_dev, sel$deviance + 2 * sel$gdf)
  expect_equal(sel$ubre, sel$deviance / n + 2 * sel$gdf / n - 1)
  expect_true(sel$fit$ggcv$converged)

  # Probes are y + eps (r_j + 1) / 2 around the shared reference y + eps / 2:
  # no response falls below y >= 0, every row of every probe (zero rows
  # included) sits exactly eps / 2 from the reference, i.e. nothing is
  # clipped, and the steps are the Rademacher probes 1..M_final.
  eps <- 1e-3 * max(sqrt(mean(bt$y^2)), 1)
  dy <- vapply(seen$y, function(v) v - bt$y, numeric(n))
  expect_identical(ncol(dy), sel$n_fits)
  expect_gte(min(dy), 0)
  is_est <- colSums(abs(dy)) == 0
  is_ref <- colSums(abs(dy - eps / 2)) < 1e-9
  dp <- dy[, !is_est & !is_ref, drop = FALSE]
  expect_gt(ncol(dp), 0L)
  expect_lt(max(abs(abs(dp - eps / 2) - eps / 2)), 1e-9)
  expect_gt(mean(dp[zero, ] > 0), 0.4)
  st_probe <- list(seed = 1L, n_rows = n, support = rep(TRUE, n), n_eff = n)
  bank <- vapply(1:16, .ggcv_probe, numeric(n), st = st_probe)
  r <- 2 * dp / eps - 1
  in_bank <- apply(r, 2L, function(rj) any(colSums(abs(bank - rj)) < 1e-6))
  expect_true(all(in_bank))

  # Exact dense UBRE (deviance units) on the same rows. The curve rises by
  # more than one deviance unit within +-0.4 of its argmin, which is what
  # makes the 0.4 tolerance on the selected theta meaningful.
  bs <- build_marginal_bases(bt$X, k = K, degree = 3L)
  B <- .ge_row_tensor(bs$basis[[1L]], bs$basis[[2L]])
  grid <- seq(-3, 3, by = 0.05)
  ubre_dev <- vapply(grid, function(t) {
    d <- .ge_dense_poisson(bt$y, B, glam_penalty(c(K, K), rep(10^t, 2L)), bt$offset)
    glm_deviance(stats::poisson(), bt$y, d$mu) + 2 * d$edf
  }, numeric(1))
  i <- which.min(ubre_dev)
  side <- c(max(1L, i - 8L), min(length(grid), i + 8L))
  expect_gt(min(ubre_dev[side] - ubre_dev[i]), 1)
  expect_lt(abs(sel$theta - grid[i]), 0.4)
  # at saturating rank the selected fit is the dense fit with these offsets
  d_sel <- .ge_dense_poisson(bt$y, B, glam_penalty(c(K, K), sel$lambda), bt$offset)
  expect_equal(as.numeric(sel$fit$fitted.values), d_sel$mu, tolerance = 1e-6)
  # the nearest smoother isotropic evaluation is a final candidate, and here
  # it is the most regular one that the minimum does not beat
  p <- sel$paired
  expect_true("smoother" %in% p$candidate)
  expect_identical(p$candidate[p$winner], "smoother")
  expect_false(p$best[p$winner])
})

test_that("paired SE of two nearby lambdas is far below the unpaired SE", {
  g <- .ge_grid()
  set.seed(21)
  eta <- sin(2 * pi * g$X[, 1]) + cos(2 * pi * g$X[, 2])
  y <- stats::rpois(nrow(g$X), exp(log(0.8) + eta))
  # fixed 20-iteration budget: under the adaptive rule the rank-3 fit at
  # theta = 0 drifts slowly (relative change ~5e-7 per iteration) and stays
  # unconverged at fit_tol; the paired-SE comparison does not depend on it
  st <- .ge_state(y, g$X, stats::poisson(), rank = 3L, k = 6L,
                  control = tt_control(pirls_maxit = 20L, seed = 1L,
                                       compute_edf = FALSE),
                  budget = "fixed")
  ev <- .ggcv_eval_batch(matrix(c(0, 0.25), ncol = 1L), 16L, st, "pair")
  a <- ev[[1L]]
  b <- ev[[2L]]
  expect_true(a$converged && b$converged)
  p <- .ggcv_paired(a, b, tol = 1)
  expect_identical(p$J, 16L)
  expect_equal(p$diff, a$score_dev - b$score_dev)
  expect_equal(p$se, stats::sd(a$g * a$contrib - b$g * b$contrib) / sqrt(16))
  unpaired <- sqrt((a$g * a$gdf_se)^2 + (b$g * b$gdf_se)^2)
  expect_lt(p$se, unpaired)
  # common random numbers, not chance: independent probe banks would give a
  # ratio near 1
  expect_lt(p$se / unpaired, 0.5)
  expect_identical(p$improves, p$diff > max(1, 2 * p$se))
})

test_that("parsimony: isotropic truth stays isotropic, anisotropic is grouped", {
  select <- function(seed, truth, budget = "adaptive", max_sweeps = 400L) {
    dat <- .ge_gauss_xy(seed, truth)
    tt_ggcv(dat$y, dat$X, rank = 6L, k = 6L, groups = 1:2, theta_lower = -3,
            theta_upper = 3, budget = budget,
            control = tt_control(max_sweeps = max_sweeps, seed = 1L,
                                 compute_edf = FALSE))
  }
  # the winner is converged, its GDF is stable and survives the doubled
  # budget at once
  verified <- function(sel) {
    w <- sel$paired$winner
    expect_true(sel$paired$converged[w])
    expect_true(sel$paired$stable[w])
    expect_identical(nrow(sel$budget_checks), 1L)
    expect_true(sel$budget_checks$ok)
  }
  # the nearest stable isotropic evaluation on the smoother side of the
  # isotropic optimum is a final candidate
  has_smoother <- function(sel) {
    s <- sel$search
    t_iso <- sel$paired$theta1[sel$paired$candidate == "isotropic"]
    up <- s$stage != "budget" & s$theta1 == s$theta2 & s$stable &
      s$theta1 > t_iso + 1e-9
    expect_identical(sel$paired$theta1[sel$paired$candidate == "smoother"],
                     min(s$theta1[up]))
  }
  iso <- select(62L, .ge_truth_iso)
  expect_identical(iso$decision, "isotropic")
  expect_identical(iso$theta[1L], iso$theta[2L])
  expect_true(any(iso$search$stage == "pattern"))
  # a grouped final candidate scores lower, but not by max(tol, 2 * SE): the
  # final choice (isotropic before grouped), not the search, keeps the
  # isotropic point
  best <- iso$paired[which.min(iso$paired$score_dev), ]
  expect_false(best$theta1 == best$theta2)
  expect_lt(best$score_dev, iso$score_dev)
  expect_lte(-best$diff_vs_winner, max(1, 2 * best$se_vs_winner))
  verified(iso)
  has_smoother(iso)

  # Anisotropic truth at rank = k, so the converged TT fit is the dense
  # tensor-product fit: the grouped decision must hold with GDFs that are
  # accurate (each final candidate's GDF is the exact mean of r_j' S r_j
  # over its 16 probes) and converged, under the adaptive budget and a fixed
  # one. By the exact criterion the grouped winner beats the isotropic
  # candidate by about 19 units. (With heavy smoothing on x2 the core systems
  # are gauge singular; their solve once made GDFs jump and the decision
  # flip with the budget.)
  dat <- .ge_gauss_xy(61L, .ge_truth_aniso)
  n <- length(dat$y)
  S_of <- .ge_dense_smoother(dat$X, 6L)
  R <- .ge_probes(n, 16L)
  exact <- function(theta) {
    S <- S_of(10^theta)
    dev <- sum((dat$y - S %*% dat$y)^2)
    c(mc = mean(colSums(R * (S %*% R))),
      score = n * log(n * dev / (n - sum(diag(S)))^2))
  }
  for (b in list(list("adaptive", 400L), list("fixed", 40L))) {
    aniso <- select(61L, .ge_truth_aniso, budget = b[[1L]], max_sweeps = b[[2L]])
    expect_identical(aniso$decision, "grouped")
    expect_gt(aniso$theta[2L] - aniso$theta[1L], 1)
    p <- aniso$paired
    iso_row <- p[p$candidate == "isotropic", ]
    expect_gt(iso_row$diff_vs_winner, max(1, 2 * iso_row$se_vs_winner))
    verified(aniso)
    has_smoother(aniso)
    ex <- t(vapply(seq_len(nrow(p)), function(i) {
      exact(c(p$theta1[i], p$theta2[i]))
    }, numeric(2)))
    expect_lt(max(abs(p$gdf - ex[, "mc"])), 0.01)
    expect_gt(ex[p$candidate == "isotropic", "score"] - ex[p$winner, "score"], 10)
  }
})

test_that("a final-stage re-evaluation runs only M_final - M_search + 2 fits", {
  g <- .ge_grid()
  set.seed(21)
  y <- stats::rpois(nrow(g$X), exp(log(0.2) + sin(2 * pi * g$X[, 1]) +
                                     cos(2 * pi * g$X[, 2])))
  ctrl <- tt_control(pirls_maxit = 10L, seed = 1L, compute_edf = FALSE)
  st <- .ge_state(y, g$X, stats::poisson(), rank = 3L, k = 6L, control = ctrl)
  M_search <- 4L
  M_final <- 16L
  # search: estimator, reference and M_search probes; final: estimator,
  # reference and only the probes M_search + 1..M_final
  a <- .ggcv_eval_batch(matrix(-0.5), M_search, st, "search")[[1L]]
  b <- .ggcv_eval_batch(matrix(-0.5), M_final, st, "final", keep_fit = TRUE)[[1L]]
  expect_identical(.ggcv_log_table(st)$n_fits,
                   c(M_search + 2L, M_final - M_search + 2L))
  expect_s3_class(b$fit, "ttpspline")
  # probes are nested, and the extended cache entry equals a fresh evaluation
  expect_identical(b$contrib[seq_len(M_search)], a$contrib)
  fresh <- .ggcv_eval_batch(matrix(-0.5), M_final,
                            .ge_state(y, g$X, stats::poisson(), rank = 3L, k = 6L,
                                      control = ctrl), "fresh")[[1L]]
  expect_identical(b$contrib, fresh$contrib)
  expect_identical(b$gdf, fresh$gdf)
  # a full cache hit runs nothing and logs nothing
  again <- .ggcv_eval_batch(matrix(-0.5), M_final, st, "again")[[1L]]
  expect_identical(again$n_fits, 0L)
  expect_identical(nrow(.ggcv_log_table(st)), 2L)

  # the same accounting in the search log of a whole selection
  dat <- .ge_gauss_xy(62L, .ge_truth_iso)
  sel <- tt_ggcv(dat$y, dat$X, rank = 6L, k = 6L, theta_lower = -3,
                 theta_upper = 3, n_grid = 5L, n_refine = 1L, M_search = 3L,
                 M_final = 7L,
                 control = tt_control(max_sweeps = 20L, seed = 1L,
                                      compute_edf = FALSE))
  fin <- sel$search$n_fits[sel$search$stage == "final"]
  expect_length(fin, 3L)
  expect_true(all(fin == 7L - 3L + 2L))
  expect_identical(sel$search$n_fits[sel$search$stage == "budget"], 7L + 2L)
  expect_identical(sum(sel$search$n_fits), sel$n_fits)
})

test_that("evaluations at a two-iteration PIRLS budget are flagged unconverged", {
  g <- .ge_grid()
  set.seed(21)
  y <- stats::rpois(nrow(g$X), exp(log(0.2) + sin(2 * pi * g$X[, 1]) +
                                     cos(2 * pi * g$X[, 2])))
  res <- .ge_warnings(tt_ggcv(
    y, g$X, rank = 3L, k = 6L, family = stats::poisson(), n_grid = 5L,
    n_refine = 0L, M_search = 2L, M_final = 4L,
    control = tt_control(pirls_maxit = 2L, seed = 1L, compute_edf = FALSE,
                         warn_lambda_boundary = FALSE)
  ))
  sel <- res$value
  # every evaluation at the two-iteration budget (the budget-check row runs
  # four iterations)
  at_budget <- sel$search[sel$search$stage != "budget", ]
  expect_gt(nrow(at_budget), 0L)
  expect_false(any(at_budget$converged))
  expect_true(all(at_budget$last_rel_change > 1e-5))
  expect_false(any(sel$paired$converged))
  expect_false(sel$fit$ggcv$converged)
  # the winner warning fires, and every warning is about convergence
  has <- function(txt) grepl(txt, res$warnings, fixed = TRUE)
  expect_true(any(has("selected lambda is not converged")))
  expect_true(any(has("pirls_maxit")))
  expect_true(all(has("not converged")))
})

test_that("ttps(lambda = 'gGCV') on scattered rows refits with the user control", {
  dat <- .ge_scattered()
  ctl <- tt_control(max_sweeps = 30L, pirls_maxit = 30L, compute_edf = TRUE,
                    seed = 1L, ggcv_max_sweeps = 6L, ggcv_pirls_maxit = 8L,
                    ggcv_n_grid = 5L, ggcv_n_refine = 1L, ggcv_M_search = 3L,
                    ggcv_M_final = 6L, warn_lambda_boundary = FALSE)
  ctl0 <- ctl
  ctl0$ggcv_refit <- FALSE
  cases <- list(
    list(y = dat$y_gauss, family = stats::gaussian(), offset = NULL),
    list(y = dat$y_pois, family = stats::poisson(), offset = dat$offset)
  )
  for (cs in cases) {
    select <- function(control) {
      .ge_quiet_budget(ttps(cs$y, dat$X, family = cs$family, offset = cs$offset,
                            rank = 2L, k = 6L, lambda = "gGCV",
                            control = control))
    }
    set.seed(99)
    u <- stats::runif(1)
    set.seed(99)
    fit <- select(ctl)
    expect_identical(stats::runif(1), u)
    expect_s3_class(fit, "ttpspline")
    expect_identical(fit$lambda_method, "gGCV")
    expect_identical(fit$ggcv$mode, "scattered")
    expect_true(is.finite(fit$ggcv$gdf))
    expect_identical(fit$call$lambda, "gGCV")  # the caller's call, not the refit's
    # the refit runs the caller's control, not the fixed selection budget
    expect_true(fit$ggcv$refit)
    expect_true(is.finite(fit$edf))
    expect_identical(fit$control$tol, ctl$tol)
    expect_identical(fit$control$max_sweeps, 30L)
    if (!is.null(cs$offset)) expect_equal(fit$offset, cs$offset)

    set.seed(99)
    fit0 <- select(ctl0)
    expect_identical(stats::runif(1), u)
    expect_s3_class(fit0, "ttpspline")
    expect_identical(fit0$lambda_method, "gGCV")
    expect_false(fit0$ggcv$refit)
    expect_identical(fit0$lambda, fit$lambda)
    # ggcv_refit = FALSE returns the engine's estimator fit: the cold fit at
    # the winner's budget B (at most the caps) with tol = 0, the map the
    # probes differentiated
    B <- fit0$ggcv$budget$B
    expect_identical(fit0$control$tol, 0)
    if (identical(cs$family$family, "poisson")) {
      expect_identical(fit0$control$pirls_maxit, B)
      expect_lte(B, 8L)
    } else {
      expect_identical(fit0$control$max_sweeps, B)
      expect_lte(B, 6L)
    }
    expect_false(fit0$control$compute_edf)
    expect_true(is.na(fit0$edf))
    again <- ttps(cs$y, dat$X, family = cs$family, offset = cs$offset,
                  rank = 2L, k = 6L, lambda = fit0$lambda, optimizer = "ALS",
                  control = fit0$control)
    expect_identical(fitted(again), fitted(fit0))
  }
})

test_that("tt_ggcv_poisson(mode = 'working') warns and uses the exact map", {
  dat <- .ge_scattered()
  ctl <- tt_control(pirls_maxit = 8L, compute_edf = FALSE, seed = 1L,
                    ggcv_n_grid = 5L, ggcv_n_refine = 1L, ggcv_M_search = 2L,
                    ggcv_M_final = 4L, warn_lambda_boundary = FALSE)
  run <- function(...) {
    .ge_quiet_budget(tt_ggcv_poisson(dat$y_pois, dat$X, rank = 2L, k = 6L,
                                     offset = dat$offset, M_search = 2L,
                                     M_final = 4L, ...))
  }
  expect_warning(sel <- run(mode = "working", control = ctl), "deprecated")
  expect_s3_class(sel$fit, "ttpspline")
  expect_length(sel$lambda, 2L)
  expect_true(all(is.finite(sel$lambda)))
  expect_true(is.finite(sel$gdf))
  expect_identical(sel$criterion, "ubre")
  expect_identical(sel$ggcv_glm_mode, "algorithmic")
  # the working-response proxy is gone: same selection as the default
  ref <- run(control = ctl)
  expect_identical(sel$lambda, ref$lambda)
  expect_identical(sel$score_dev, ref$score_dev)
  # the same through the control knob, on every route, warned exactly once
  ctl_w <- ctl
  ctl_w$ggcv_glm_mode <- "working"
  expect_warning(run(control = ctl_w), "deprecated")
  n_dep <- function(res) sum(grepl("deprecated", res$warnings, fixed = TRUE))
  expect_identical(n_dep(.ge_warnings(run(control = ctl_w))), 1L)
  expect_identical(n_dep(.ge_warnings(run(mode = "working", control = ctl_w))), 1L)
  via_ggcv <- .ge_warnings(.ge_quiet_budget(tt_ggcv(
    dat$y_pois, dat$X, rank = 2L, k = 6L, family = stats::poisson(),
    offset = dat$offset, n_grid = 5L, n_refine = 1L, M_search = 2L,
    M_final = 4L, control = ctl_w
  )))
  expect_identical(n_dep(via_ggcv), 1L)
  expect_identical(via_ggcv$value$lambda, ref$lambda)
  via_ttps <- .ge_warnings(.ge_quiet_budget(ttps(
    dat$y_pois, dat$X, family = stats::poisson(), offset = dat$offset,
    rank = 2L, k = 6L, lambda = "gGCV", control = ctl_w
  )))
  expect_identical(n_dep(via_ttps), 1L)
})

test_that("the array routes warn on the working mode like the scattered ones", {
  # Review case (finding 9): ttps(array = TRUE, family = poisson(),
  # lambda = "gGCV") with ggcv_glm_mode = "working" gave no warning. The
  # array route never used the working proxy, so the selection is the same.
  g <- .ge_grid()
  set.seed(4)
  Y <- array(stats::rpois(72L, exp(0.3 + sin(2 * pi * g$X[, 1]))), g$ng)
  ctl <- tt_control(pirls_maxit = 6L, max_sweeps = 10L, seed = 1L,
                    compute_edf = FALSE, warn_lambda_boundary = FALSE,
                    ggcv_n_grid = 3L, ggcv_n_refine = 0L, ggcv_M_search = 2L,
                    ggcv_M_final = 2L, ggcv_pirls_maxit = 6L,
                    ggcv_max_sweeps = 10L, ggcv_budget_check = FALSE)
  ctl_w <- ctl
  ctl_w$ggcv_glm_mode <- "working"
  dep <- function(res) res$warnings[grepl("deprecated", res$warnings, fixed = TRUE)]
  on_array <- function(Y, family, control) {
    .ge_warnings(ttps(Y, axes = g$axes, array = TRUE, family = family,
                      rank = 2L, k = 5L, lambda = "gGCV", control = control))
  }
  via_ttps <- on_array(Y, stats::poisson(), ctl_w)
  expect_length(dep(via_ttps), 1L)
  expect_match(dep(via_ttps), "ttps(array = TRUE, lambda = \"gGCV\")", fixed = TRUE)
  ref <- on_array(Y, stats::poisson(), ctl)
  expect_length(dep(ref), 0L)
  expect_identical(via_ttps$value$lambda, ref$value$lambda)
  via_fun <- .ge_warnings(tt_ggcv_array(
    Y, axes = g$axes, family = stats::poisson(), rank = 2L, k = 5L,
    n_grid = 3L, n_refine = 0L, M_search = 2L, M_final = 2L,
    budget_check = FALSE, control = ctl_w
  ))
  expect_length(dep(via_fun), 1L)
  expect_match(dep(via_fun), "tt_ggcv_array()", fixed = TRUE)
  # Gaussian data ignore the knob without a warning, as in tt_ggcv()
  expect_length(dep(on_array(Y + 0.5, stats::gaussian(), ctl_w)), 0L)
})

# ---------------------------------------------------------------------------
# noise-aware final choice, budget-verified winner, adaptive budget, tt_gdf()
# ---------------------------------------------------------------------------

# Synthetic final candidate for .ggcv_choose(): 16 probe contributions made of
# a part shared by all candidates (common random numbers, cancels in paired
# differences) plus an own +-`noise` part; g = 2 (UBRE slope). The paired SE
# of two candidates is then 2 * |noise_a - noise_b| * sd(+-1) / 4.
.ge_cand <- function(theta, score_dev, noise = 0, converged = TRUE) {
  J <- 16L
  list(theta = theta, iso = diff(range(theta)) < 1e-9, score_dev = score_dev,
       g = 2, contrib = 20 + 5 * sin(seq_len(J)) + noise * rep(c(1, -1), J / 2),
       converged = converged)
}

test_that(".ggcv_choose() prefers the most regular candidate the minimum does not beat", {
  # (i) large paired SE at the argmin: the smoother isotropic candidate is
  # within 2 SE of the minimum and wins; the far one does not qualify
  cs <- list(.ge_cand(-2, 0, noise = 10), .ge_cand(-1.5, 3), .ge_cand(-1, 30))
  ch <- .ggcv_choose(cs, tol = 1)
  expect_identical(ch$best, 1L)
  expect_identical(ch$winner, 2L)
  expect_identical(ch$qualifies, c(TRUE, TRUE, FALSE))
  expect_equal(ch$diff, c(0, 3, 30))
  expect_gt(2 * ch$se[2L], 3)
  # (ii) small paired SE: the minimum beats the smoother candidate by more
  # than max(tol, 2 * SE) and wins
  cs <- list(.ge_cand(-2, 0, noise = 0.2), .ge_cand(-1.5, 3))
  ch <- .ggcv_choose(cs, tol = 1)
  expect_identical(ch$winner, 1L)
  expect_identical(ch$qualifies, c(TRUE, FALSE))
  expect_lt(2 * ch$se[2L], 1)
  # (iii) a grouped minimum that does not beat the isotropic candidate by
  # max(tol, 2 * SE): isotropic wins (decision "isotropic")
  cs <- list(.ge_cand(c(-1, -1), 0.6), .ge_cand(c(-1.5, -0.5), 0, noise = 0.5))
  ch <- .ggcv_choose(cs, tol = 1)
  expect_identical(ch$best, 2L)
  expect_identical(ch$winner, 1L)
  expect_true(cs[[ch$winner]]$iso)
  # with tol = 0 (pure Monte Carlo rule) 0.6 exceeds 2 * SE = 0.52
  expect_identical(.ggcv_choose(cs, tol = 0)$winner, 2L)
  # (iv) a grouped minimum that clearly beats the isotropic candidate wins
  cs <- list(.ge_cand(c(-1, -1), 10), .ge_cand(c(-1.5, -0.5), 0, noise = 0.5))
  ch <- .ggcv_choose(cs, tol = 1)
  expect_identical(ch$winner, 2L)
  expect_false(cs[[ch$winner]]$iso)
  # (v) unconverged candidates are excluded unless all are unconverged
  cs <- list(.ge_cand(-2, 0, converged = FALSE), .ge_cand(-1.5, 4),
             .ge_cand(-1, 4.5))
  ch <- .ggcv_choose(cs, tol = 1)
  expect_false(ch$none_converged)
  expect_identical(ch$eligible, c(FALSE, TRUE, TRUE))
  expect_identical(ch$best, 2L)
  expect_identical(ch$winner, 3L)
  expect_false(ch$qualifies[1L])
  # budget-check failures are excluded as well
  expect_identical(.ggcv_choose(cs, tol = 1, exclude = 3L)$winner, 2L)
  for (q in seq_along(cs)) cs[[q]]$converged <- FALSE
  ch <- .ggcv_choose(cs, tol = 1)
  expect_true(ch$none_converged)
  expect_identical(ch$best, 1L)
  expect_identical(ch$winner, 1L)
  # nothing eligible
  expect_true(is.na(.ggcv_choose(cs, tol = 1, exclude = 1:3)$winner))
})

test_that("budget check: an unconverged winner gets one diagnostic check, no fallback", {
  skip_if_not_installed("testthat", "3.1.7")  # local_mocked_bindings()
  g <- .ge_grid()
  set.seed(21)
  y <- stats::rpois(nrow(g$X), exp(log(0.2) + sin(2 * pi * g$X[, 1]) +
                                     cos(2 * pi * g$X[, 2])))
  select <- function() {
    .ge_warnings(tt_ggcv(
      y, g$X, rank = 3L, k = 6L, family = stats::poisson(), n_grid = 5L,
      n_refine = 0L, M_search = 2L, M_final = 4L,
      control = tt_control(pirls_maxit = 2L, seed = 1L, compute_edf = FALSE,
                           warn_lambda_boundary = FALSE)
    ))
  }
  res <- select()
  sel <- res$value
  expect_false(any(sel$paired$converged))
  # the first winner's check runs (its GDF moves at twice the budget), but
  # no fallback check is spent on the other unconverged candidates, which
  # could not be verified either
  bc <- sel$budget_checks
  expect_identical(nrow(bc), 1L)
  expect_false(bc$ok)
  expect_gt(abs(bc$rel_change), 0.01)
  expect_false(sel$budget_verified)
  expect_false(sel$fit$ggcv$budget_verified)
  expect_identical(sel$theta, bc$theta1)
  expect_identical(sel$budget$gdf_2x, bc$gdf_2x)
  expect_identical(sum(sel$search$stage == "budget"), 1L)
  # the check runs at twice the candidate's budget: 2 * cap = 4 here
  expect_identical(sel$search$budget[sel$search$stage == "budget"], 4L)
  # one warning, which reports the check
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "not converged", fixed = TRUE)
  expect_match(res$warnings, sprintf("GDF changes by %+.1f%%", 100 * bc$rel_change),
               fixed = TRUE)
  # even a check made to fail by a GDF 50% higher at 2B, which would make
  # the winner lose when rescored, moves nothing: the other candidates are
  # unconverged too and could not be checked
  local_mocked_bindings(.ggcv_eval_batch = .ge_scaled_checks(1.5, .ggcv_eval_batch))
  sel2 <- select()$value
  expect_identical(sel2$theta, sel$theta)
  expect_identical(nrow(sel2$budget_checks), 1L)
  expect_identical(sel2$paired$score_dev_used, sel2$paired$score_dev)
})

test_that("budget check: a failed check changes the choice only when the score can", {
  skip_if_not_installed("testthat", "3.1.7")  # local_mocked_bindings()
  g <- .ge_grid()
  set.seed(21)
  y <- stats::rpois(nrow(g$X), exp(log(0.8) + sin(2 * pi * g$X[, 1]) +
                                     cos(2 * pi * g$X[, 2])))
  select <- function() {
    .ge_warnings(tt_ggcv(
      y, g$X, rank = 6L, k = 6L, family = stats::poisson(), theta_lower = -3,
      theta_upper = 3, n_grid = 5L, n_refine = 1L, M_search = 3L,
      M_final = 6L, control = tt_control(seed = 1L, compute_edf = FALSE,
                                         warn_lambda_boundary = FALSE)
    ))
  }
  # generous budget (saturating rank, default caps): one check, verified
  res <- select()
  sel <- res$value
  expect_length(res$warnings, 0L)
  expect_identical(nrow(sel$budget_checks), 1L)
  expect_true(sel$budget_checks$ok)
  expect_true(sel$budget_verified)
  expect_identical(sel$budget$mode, "adaptive")
  expect_identical(sel$budget$B_2x, 2L * sel$budget$B)
  # diagnostics of the choice
  p <- sel$paired
  expect_identical(sum(p$best), 1L)
  expect_identical(sum(p$winner), 1L)
  expect_true(p$qualifies[p$winner])
  expect_identical(p$score_dev_used, p$score_dev)
  expect_identical(unname(unlist(p[p$best, "theta1"])), sel$best_theta)
  expect_true(is.character(sel$rule) && length(sel$rule) == 1L)
  expect_identical(sel$fit$ggcv$rule, sel$rule)
  s <- sel$search
  ok <- s$M_ok >= 2L
  expect_equal(s$gdf_cv[ok], s$gdf_se[ok] * sqrt(s$M_ok[ok]) / s$gdf[ok])
  expect_identical(s$stable, s$M_ok == s$M & is.finite(s$gdf_cv) &
                     s$gdf_cv <= pmax(1, 3 * sqrt(2 / pmax(s$gdf, 1))))

  # (i) the first check is made to fail with a GDF 50% higher at 2B: with
  # the less favourable of its two scores the winner loses, and the choice
  # is redone by the same rule, whose winner's own check passes. The new
  # winner is the most regular qualifying candidate under the rescored
  # values, not their minimum.
  orig <- .ggcv_eval_batch
  local_mocked_bindings(.ggcv_eval_batch = .ge_scaled_checks(1.5, orig))
  res2 <- select()
  sel2 <- res2$value
  bc <- sel2$budget_checks
  expect_identical(nrow(bc), 2L)
  expect_identical(bc$ok, c(FALSE, TRUE))
  expect_true(sel2$budget_verified)
  expect_identical(bc$theta1[1L], sel$theta)
  expect_identical(sel2$theta, bc$theta1[2L])
  p2 <- sel2$paired
  first <- p2$theta1 == sel$theta
  expect_identical(p2$score_dev_used[first], bc$score_dev_2x[1L])
  expect_gt(p2$score_dev_used[first], p2$score_dev[first])
  expect_true(p2$qualifies[p2$winner])
  expect_identical(p2$theta1[p2$winner], max(p2$theta1[p2$qualifies]))
  expect_false(p2$best[p2$winner])
  expect_identical(p2$budget_ok[first], FALSE)
  expect_identical(p2$budget_ok[p2$winner], TRUE)
  expect_identical(sel2$budget$rel_change, bc$rel_change[2L])
  expect_length(res2$warnings, 1L)
  expect_match(res2$warnings, "moved to the next candidate by the same rule",
               fixed = TRUE)

  # (ii) a GDF 2% higher at 2B fails the check (budget_tol = 1%) but cannot
  # reverse the choice: the winner is kept, unverified, and the warning
  # gives the change and the margin
  local_mocked_bindings(.ggcv_eval_batch = .ge_scaled_checks(1.02, orig))
  res3 <- select()
  sel3 <- res3$value
  expect_identical(sel3$theta, sel$theta)
  expect_identical(nrow(sel3$budget_checks), 1L)
  expect_false(sel3$budget_checks$ok)
  expect_false(sel3$budget_verified)
  expect_length(res3$warnings, 1L)
  expect_match(res3$warnings, "remains the choice", fixed = TRUE)
  expect_match(res3$warnings, sprintf("%+.1f%%", 100 * sel3$budget$rel_change),
               fixed = TRUE)
  # the runner-up it names is a candidate the rule could choose instead
  # (converged, with a stable GDF)
  p3 <- sel3$paired
  runner_up <- sub(".*runner-up ([[:alnum:]]+) at .*", "\\1", res3$warnings)
  expect_true(runner_up %in% p3$candidate[p3$converged & p3$stable & !p3$winner])
})

test_that("budget check: the first winner is kept when no converged, stable candidate can replace it", {
  skip_if_not_installed("testthat", "3.1.7")  # local_mocked_bindings()
  # Review case (finding 11, fix 2): the first winner is the only converged
  # final candidate and its check at 2B gives no usable GDF (score not
  # finite). The redone choice fell back to the unconverged candidates and
  # returned one of them, unchecked and about 4.5 score units worse, with a
  # warning that printed "GDF NA%".
  g <- .ge_grid()
  set.seed(21)
  y <- stats::rpois(nrow(g$X), exp(log(0.8) + sin(2 * pi * g$X[, 1]) +
                                     cos(2 * pi * g$X[, 2])))
  select <- function(n_final = 3L) {
    .ge_warnings(tt_ggcv(
      y, g$X, rank = 6L, k = 6L, family = stats::poisson(), theta_lower = -3,
      theta_upper = 3, n_grid = 5L, n_refine = 1L, n_final = n_final,
      M_search = 3L, M_final = 6L,
      control = tt_control(seed = 1L, compute_edf = FALSE,
                           warn_lambda_boundary = FALSE)
    ))
  }
  orig <- .ggcv_eval_batch
  local_mocked_bindings(.ggcv_eval_batch = .ge_only_first_converged(NA_real_, orig))
  res <- select()
  sel <- res$value
  bc <- sel$budget_checks
  p <- sel$paired
  expect_identical(nrow(bc), 1L)
  expect_false(bc$ok)
  expect_true(is.na(bc$rel_change))
  # the first winner (the only converged candidate) is kept, unverified
  expect_identical(sel$theta, bc$theta1)
  expect_identical(p$candidate[p$winner], "isotropic")
  expect_true(p$converged[p$winner])
  expect_false(any(p$converged[!p$winner]))
  expect_true(sel$fit$ggcv$converged)
  expect_false(sel$budget_verified)
  expect_identical(p$budget_ok[p$winner], FALSE)
  expect_identical(p$score_dev_used[p$winner], Inf)
  expect_true(p$qualifies[p$winner])  # the first choice stands
  expect_identical(sum(sel$search$stage == "budget"), 1L)
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "not converged", fixed = TRUE)
  expect_match(res$warnings, "the GDF change is not finite", fixed = TRUE)
  expect_match(res$warnings, "kept as the only eligible final candidate",
               fixed = TRUE)
  expect_false(grepl("NA%", res$warnings, fixed = TRUE))

  # The limit of three checks: with four final candidates and a GDF 50%
  # higher at every check, each checked candidate loses when rescored and
  # the fourth would be returned unchecked; the first winner is kept
  local_mocked_bindings(.ggcv_eval_batch = .ge_scaled_checks(rep(1.5, 3L), orig))
  res2 <- select(n_final = 4L)
  sel2 <- res2$value
  bc2 <- sel2$budget_checks
  p2 <- sel2$paired
  expect_identical(nrow(p2), 4L)
  expect_identical(nrow(bc2), 3L)
  expect_false(any(bc2$ok))
  expect_identical(sel2$theta, bc2$theta1[1L])
  expect_identical(sum(is.na(p2$budget_ok)), 1L)
  expect_identical(p2$budget_ok[p2$winner], FALSE)
  expect_false(sel2$budget_verified)
  expect_identical(sel2$budget$rel_change, bc2$rel_change[1L])
  expect_length(res2$warnings, 1L)
  expect_match(res2$warnings, "not converged", fixed = TRUE)
  expect_match(res2$warnings, "the limit of 3 checks was reached", fixed = TRUE)
  expect_match(res2$warnings, "the first winner is kept", fixed = TRUE)
})

test_that("budget-check warnings name only eligible runner-ups and never print NA%", {
  skip_if_not_installed("testthat", "3.1.7")  # local_mocked_bindings()
  # Review case (finding 14): a winner kept after a failed check was compared
  # with a "runner-up" that the rule could not choose, e.g. "runner-up
  # smoother at 1461.95" on a 7 x 6 x 6 x 5 Poisson array where the winner
  # was the only converged final candidate. Here the first winner is the only
  # converged candidate and its GDF is 2% higher at 2B: kept, and said to be
  # the only eligible final candidate.
  g <- .ge_grid()
  set.seed(21)
  y <- stats::rpois(nrow(g$X), exp(log(0.8) + sin(2 * pi * g$X[, 1]) +
                                     cos(2 * pi * g$X[, 2])))
  orig <- .ggcv_eval_batch
  local_mocked_bindings(.ggcv_eval_batch = .ge_only_first_converged(1.02, orig))
  res <- .ge_warnings(tt_ggcv(
    y, g$X, rank = 6L, k = 6L, family = stats::poisson(), theta_lower = -3,
    theta_upper = 3, n_grid = 5L, n_refine = 1L, M_search = 3L, M_final = 6L,
    control = tt_control(seed = 1L, compute_edf = FALSE,
                         warn_lambda_boundary = FALSE)
  ))
  sel <- res$value
  expect_identical(nrow(sel$budget_checks), 1L)
  expect_false(sel$budget_checks$ok)
  expect_identical(sel$theta, sel$budget_checks$theta1)
  expect_false(sel$budget_verified)
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "remains the choice", fixed = TRUE)
  expect_match(res$warnings, "as the only eligible final candidate", fixed = TRUE)
  expect_false(grepl("runner-up", res$warnings, fixed = TRUE))
  expect_match(res$warnings, sprintf("GDF changes by %+.1f%%",
                                     100 * sel$budget$rel_change), fixed = TRUE)

  # an unconverged winner whose check at 2B gives no usable GDF: the change
  # reads "not finite" (it printed "NA%")
  set.seed(21)
  y2 <- stats::rpois(nrow(g$X), exp(log(0.2) + sin(2 * pi * g$X[, 1]) +
                                      cos(2 * pi * g$X[, 2])))
  local_mocked_bindings(.ggcv_eval_batch = .ge_scaled_checks(NA_real_, orig))
  res2 <- .ge_warnings(tt_ggcv(
    y2, g$X, rank = 3L, k = 6L, family = stats::poisson(), n_grid = 5L,
    n_refine = 0L, M_search = 2L, M_final = 4L,
    control = tt_control(pirls_maxit = 2L, seed = 1L, compute_edf = FALSE,
                         warn_lambda_boundary = FALSE)
  ))
  expect_true(is.na(res2$value$budget$rel_change))
  expect_length(res2$warnings, 1L)
  expect_match(res2$warnings, "selected lambda is not converged", fixed = TRUE)
  expect_match(res2$warnings, "the GDF change is not finite when the budget doubles",
               fixed = TRUE)
  expect_false(grepl("NA%", res2$warnings, fixed = TRUE))
  # the formatting helpers
  expect_identical(.ggcv_pct_txt(-0.0183), "-1.8%")
  expect_identical(.ggcv_pct_txt(NA_real_), "not finite")
  expect_identical(.ggcv_pct_txt(Inf), "not finite")
  expect_identical(.ggcv_num_txt(Inf, "%.2f"), "not finite")
  expect_identical(.ggcv_num_txt(NA_real_), "not available")
  expect_identical(.ggcv_num_txt(1470.6012, "%.2f"), "1470.60")
})

test_that("budget check: the seed-61 grouped optimum survives a small GDF change", {
  skip_if_not_installed("testthat", "3.1.7")  # local_mocked_bindings()
  # Review case: at fixed budgets where the grouped optimum's GDF moved by
  # -5.6% at twice the budget (worth about 1.4 score units), the fallback
  # dropped it and returned an isotropic candidate it led by about 20 units
  # (the exact criterion agrees, see the parsimony test). The check of the
  # grouped winner is made to fail by -5.6% here.
  dat <- .ge_gauss_xy(61L, .ge_truth_aniso)
  local_mocked_bindings(.ggcv_eval_batch = .ge_scaled_checks(1 / 0.944,
                                                             .ggcv_eval_batch))
  res <- .ge_warnings(tt_ggcv(
    dat$y, dat$X, rank = 6L, k = 6L, groups = 1:2, theta_lower = -3,
    theta_upper = 3, control = tt_control(seed = 1L, compute_edf = FALSE)
  ))
  sel <- res$value
  expect_identical(sel$decision, "grouped")
  expect_gt(sel$theta[2L] - sel$theta[1L], 1)
  bc <- sel$budget_checks
  expect_identical(nrow(bc), 1L)
  expect_false(bc$ok)
  expect_equal(bc$rel_change, -0.056, tolerance = 1e-3)
  expect_false(sel$budget_verified)
  # with its less favourable (2B) score it still leads the isotropic
  # candidate by far more than max(tol, 2 * SE)
  p <- sel$paired
  expect_identical(p$score_dev_used[p$winner], bc$score_dev_2x)
  expect_gt(p$score_dev_used[p$candidate == "isotropic"] - p$score_dev_used[p$winner], 10)
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "-5.6%", fixed = TRUE)
  expect_match(res$warnings, "remains the choice", fixed = TRUE)
})

test_that("adaptive budget: B follows lambda, a binding cap, GDF matches a 4x budget", {
  set.seed(11)
  n <- 200
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  y <- 0.4 + sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + stats::rnorm(n, sd = 0.3)
  select <- function(cap) {
    .ge_quiet_budget(tt_ggcv(
      y, X, rank = 2L, k = 8L, theta_lower = -3, theta_upper = 0, n_grid = 4L,
      n_refine = 1L, M_search = 2L, M_final = 6L,
      control = tt_control(max_sweeps = cap, seed = 1L, compute_edf = FALSE)
    ))
  }
  # B from the pilot: ceiling(1.5 n_iter) + 2 under the cap when converged,
  # the cap otherwise; doubled per probe-validation round
  B_of <- function(at, cap) {
    b0 <- ifelse(at$converged, pmin(cap, ceiling(1.5 * at$n_iter) + 2), cap)
    as.integer(b0 * 2^at$reevals)
  }
  cap <- 300L
  sel <- select(cap)
  at <- sel$search[sel$search$stage != "budget", ]
  # (b) the budget of the estimator, reference and probe fits follows lambda
  expect_gt(length(unique(at$budget)), 1L)
  expect_true(all(at$converged))
  expect_identical(at$budget, B_of(at, cap))
  # the budget check runs at 2B, above the cap if need be
  expect_identical(sel$budget$B_2x, 2L * sel$budget$B)
  expect_true(sel$fit$ggcv$converged)
  # (c) a cap that binds: the slowly mixing evaluations at larger lambda run
  # the whole cap without settling and are unconverged (not eligible); the
  # others settled before it
  cap2 <- 60L
  sel2 <- select(cap2)
  at2 <- sel2$search[sel2$search$stage != "budget", ]
  expect_true(any(!at2$converged))
  expect_true(all(at2$n_iter[!at2$converged] == cap2))
  expect_true(all(at2$n_iter[at2$converged] < cap2))
  expect_identical(at2$budget, B_of(at2, cap2))
  expect_true(all(at2$budget[at2$reevals == 0L] <= cap2))
  p2 <- sel2$paired
  expect_true(any(!p2$converged))
  expect_true(p2$converged[p2$winner])
  expect_false(any(p2$qualifies & !p2$converged))
  # (a) at the selected lambda, the adaptive GDF (budget B) is within 2% of a
  # fixed run at 4B with the same probes; a fixed 20-sweep run is not
  gdf_at <- function(...) {
    tt_gdf(y, X, lambda = sel$lambda, rank = 2L, k = 8L, M = 6L, ...)
  }
  ad <- gdf_at(control = tt_control(max_sweeps = cap, seed = 1L, compute_edf = FALSE))
  expect_equal(ad$gdf, sel$gdf, tolerance = 1e-8)
  expect_identical(ad$budget, sel$budget$B)
  f4 <- gdf_at(budget = "fixed",
               control = tt_control(max_sweeps = 4L * ad$budget, seed = 1L,
                                    compute_edf = FALSE))
  f20 <- gdf_at(budget = "fixed",
                control = tt_control(max_sweeps = 20L, seed = 1L,
                                     compute_edf = FALSE))
  expect_lt(abs(ad$gdf - f4$gdf) / f4$gdf, 0.02)
  expect_gt(abs(f20$gdf - f4$gdf) / f4$gdf, 0.02)
})

test_that("tt_gdf() equals tt_gdf_array() on a grid and the dense EDF at saturating rank", {
  g <- .ge_grid()
  lam <- c(0.5, 2)
  eta <- sin(2 * pi * g$X[, 1]) + cos(2 * pi * g$X[, 2])
  set.seed(23)
  n <- nrow(g$X)
  ctl <- tt_control(max_sweeps = 60L, pirls_maxit = 30L, seed = 1L,
                    compute_edf = FALSE)
  cases <- list(
    list(family = stats::gaussian(), y = eta + stats::rnorm(n, sd = 0.3)),
    list(family = stats::poisson(), y = stats::rpois(n, exp(log(0.5) + eta)))
  )
  bank <- matrix(sample(c(-1, 1), n * 6L, replace = TRUE), n, 6L)
  for (cs in cases) {
    sca <- tt_gdf(cs$y, g$X, lambda = lam, family = cs$family, rank = 3L,
                  k = 6L, M = 6L, control = ctl)
    arr <- tt_gdf_array(array(cs$y, g$ng), lambda = lam, axes = g$axes,
                        family = cs$family, rank = 3L, k = 6L, M = 6L,
                        control = ctl)
    expect_identical(sca$budget, arr$budget)
    expect_equal(sca$contrib, arr$contrib, tolerance = 1e-6)
    expect_equal(sca$gdf, arr$gdf, tolerance = 1e-6)
    expect_equal(sca$deviance, arr$deviance, tolerance = 1e-6)
    pois <- identical(cs$family$family, "poisson")
    expect_identical(sca$n_fits, .ge_adaptive_fits(sca, 6L, b0 = if (pois) 10L else 20L,
                                                   cap = if (pois) 30L else 60L))
    expect_identical(arr$n_fits, sca$n_fits)
    # a user bank equal to the seeded probes reproduces the seeded result
    st <- list(seed = 1L, n_rows = n, support = rep(TRUE, n), n_eff = n)
    seeded <- vapply(1:6, .ggcv_probe, numeric(n), st = st)
    same <- tt_gdf(cs$y, g$X, lambda = lam, family = cs$family, rank = 3L,
                   k = 6L, probes = seeded, control = ctl)
    expect_identical(same$contrib, sca$contrib)
    # another bank: scattered rows and the array still agree
    s2 <- tt_gdf(cs$y, g$X, lambda = lam, family = cs$family, rank = 3L,
                 k = 6L, probes = bank, control = ctl)
    a2 <- tt_gdf_array(array(cs$y, g$ng), lambda = lam, axes = g$axes,
                       family = cs$family, rank = 3L, k = 6L, probes = bank,
                       control = ctl)
    expect_length(s2$contrib, 6L)
    expect_equal(s2$gdf, a2$gdf, tolerance = 1e-6)
    expect_false(isTRUE(all.equal(s2$contrib, sca$contrib)))
  }

  # unit probes at saturating rank: the exact dense EDF
  set.seed(31)
  n <- 80L
  K <- 5L
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  y <- sin(2 * pi * X[, 1]) + cos(2 * pi * X[, 2]) + stats::rnorm(n, sd = 0.3)
  bs <- build_marginal_bases(X, k = K, degree = 3L)
  B <- .ge_row_tensor(bs$basis[[1L]], bs$basis[[2L]])
  A <- crossprod(B) + glam_penalty(c(K, K), lam)
  edf <- sum(diag(solve(A, crossprod(B))))
  u <- tt_gdf(y, X, lambda = lam, rank = K, k = K, probes = "unit")
  expect_true(u$converged)
  expect_equal(u$gdf, edf, tolerance = 1e-6)
  # the pilot, the estimator at B (the base of the differences) unless the
  # last pilot run is it, and one fit per row
  expect_identical(u$n_fits, .ge_adaptive_fits(u, n, b0 = 20L, cap = 400L) - 1L)
  expect_identical(u$gdf_cv, 0)
})

test_that("user probe banks are validated; zero-weight cells are ignored", {
  g <- .ge_grid()
  n <- nrow(g$X)
  set.seed(24)
  y <- sin(2 * pi * g$X[, 1]) + stats::rnorm(n, sd = 0.3)
  ctl <- tt_control(max_sweeps = 20L, seed = 1L, compute_edf = FALSE)
  gdf <- function(probes) {
    tt_gdf(y, g$X, lambda = 1, rank = 2L, k = 5L, probes = probes, control = ctl)
  }
  bank <- matrix(sample(c(-1, 1), n * 3L, replace = TRUE), n, 3L)
  expect_error(gdf(bank[-1L, ]), "one row per observation")
  expect_error(gdf(replace(bank, 5L, 0.5)), "only -1 and \\+1")
  expect_error(gdf(replace(bank, 5L, NA)), "only -1 and \\+1")
  expect_error(gdf(matrix(numeric(0), n, 0L)), "at least one column")
  expect_error(gdf("gaussian"), "should be one of")
  expect_error(tt_gdf(y, g$X, lambda = c(1, 2, 3), rank = 2L, k = 5L,
                      control = ctl), "length d")
  # array with zero-weight cells: their rows of the bank are ignored
  W <- array(1, g$ng)
  W[1:2, 1:2] <- 0
  b0 <- bank
  b0[which(W == 0), ] <- 7
  a1 <- tt_gdf_array(array(y, g$ng), lambda = 1, axes = g$axes, rank = 2L,
                     k = 5L, weights = W, probes = b0, control = ctl)
  a2 <- tt_gdf_array(array(y, g$ng), lambda = 1, axes = g$axes, rank = 2L,
                     k = 5L, weights = W, probes = bank, control = ctl)
  expect_identical(a1$contrib, a2$contrib)
  expect_identical(a1$n_eff, sum(W > 0))
})

# ---------------------------------------------------------------------------
# numerics: the global-mode core solve, probe seeds, UBRE on the AIC scale
# and the Gaussian probe step
# ---------------------------------------------------------------------------

test_that("global-mode core solve: exact when well conditioned, minimum norm when singular", {
  set.seed(3)
  m <- 12L
  solve_core <- function(S, b) {
    update_lambda_fixed(list(S = S, P = diag(m), P0 = matrix(0, m, m), b = b,
                             lambda0 = 0))$g
  }
  A <- matrix(stats::rnorm(40L * m), 40L, m)
  b <- stats::rnorm(m)
  expect_equal(solve_core(crossprod(A), b), as.numeric(solve(crossprod(A), b)),
               tolerance = 1e-10)
  # rank-7 design: S = A'A is singular along 5 gauge directions (zero up to
  # roundoff). Cholesky / LU / ridge used to leave roundoff over a near-zero
  # pivot there, or a ridge bias; the solve is now the minimum-norm one and
  # moves continuously with the data.
  A <- matrix(stats::rnorm(40L * 7L), 40L, 7L) %*% matrix(stats::rnorm(7L * m), 7L, m)
  z <- stats::rnorm(40L)
  sv <- svd(A)
  min_norm <- function(z) {
    as.numeric(sv$v[, 1:7] %*% (crossprod(sv$u[, 1:7], z) / sv$d[1:7]))
  }
  g <- solve_core(crossprod(A), as.numeric(crossprod(A, z)))
  expect_equal(g, min_norm(z), tolerance = 1e-8)
  dz <- 1e-9 * stats::rnorm(40L)
  g2 <- solve_core(crossprod(A), as.numeric(crossprod(A, z + dz)))
  expect_lt(max(abs(g2 - g)), 1e-6)
  expect_error(solve_core(replace(crossprod(A), 1L, NA), b), "non-finite")
})

test_that("GDF at saturating rank under heavy smoothing matches the dense smoother", {
  # The Exp B smoke design (20 x 20 grid, rank = k = 10) with one margin
  # heavily smoothed, and the test-6 point (scattered rows, rank = k = 6).
  # The TT core systems are gauge singular there (rcond ~ 1e-20); the old
  # Cholesky / LU / ridge cascade made the fixed-lambda map jump with
  # roundoff and gave probe contributions off by up to 1e4. Each probe
  # contribution must match r_j' S r_j of the dense tensor-product smoother.
  paired_err <- function(y, X, K, lam, seed, max_sweeps) {
    n <- length(y)
    set.seed(11)
    bank <- matrix(sample(c(-1, 1), n * 8L, replace = TRUE), n, 8L)
    ctl <- tt_control(max_sweeps = max_sweeps, tol = 1e-8, seed = seed,
                      compute_edf = FALSE, warn_lambda_boundary = FALSE)
    g <- tt_gdf(y, X, lambda = lam, rank = K, k = K, probes = bank,
                fit_tol = 1e-8, control = ctl)
    bs <- build_marginal_bases(X, k = K, degree = 3L)
    B <- .ge_row_tensor(bs$basis[[1L]], bs$basis[[2L]])
    S <- B %*% solve(crossprod(B) + glam_penalty(c(K, K), lam), t(B))
    max(abs(g$contrib - colSums(bank * (S %*% bank))))
  }
  ax <- seq(0, 1, length.out = 20L)
  X <- as.matrix(expand.grid(x1 = ax, x2 = ax))
  f <- sin(pi * X[, 1]) * sin(pi * X[, 2]) +
    0.5 * sin(2 * pi * X[, 1]) * cos(pi * X[, 2])
  set.seed(20260818)
  y <- f + stats::rnorm(nrow(X), sd = sqrt(stats::var(f) / 5))
  cells <- list(c(-1, 4), c(2, 4), c(4, 3), c(1, 4), c(3, 4))
  seeds <- c(101L, 102L, 103L, 105L, 108L)
  for (i in seq_along(cells)) {
    expect_lt(paired_err(y, X, 10L, 10^cells[[i]], seeds[i], 25L), 0.5)
  }
  dat <- .ge_gauss_xy(61L, .ge_truth_aniso)
  expect_lt(paired_err(dat$y, dat$X, 6L, 10^c(-2.162954, 0.337046), 1L, 400L), 0.5)
})

test_that("probe seeds above the integer range work; old probes kept; draw errors surface", {
  # seeds >= 21475 used to overflow seed * 100003 + j to NA and stop every
  # route (here a date-like seed, in the probe seed and in control$seed)
  dat <- .ge_scattered()
  ctl <- tt_control(max_sweeps = 30L, seed = 1L, compute_edf = FALSE)
  g <- tt_gdf(dat$y_gauss, dat$X, lambda = 1, rank = 2L, k = 5L, M = 2L,
              seed = 20261002, control = ctl)
  expect_true(is.finite(g$gdf))
  expect_identical(g$M_ok, 2L)
  fit <- .ge_quiet_budget(ttps(
    dat$y_gauss, dat$X, rank = 2L, k = 5L, lambda = "gGCV",
    control = tt_control(seed = 20261002L, compute_edf = FALSE,
                         warn_lambda_boundary = FALSE, ggcv_n_grid = 3L,
                         ggcv_n_refine = 0L, ggcv_M_search = 2L,
                         ggcv_M_final = 2L, ggcv_max_sweeps = 30L)
  ))
  expect_s3_class(fit, "ttpspline")
  expect_identical(fit$lambda_method, "gGCV")
  expect_true(is.finite(fit$ggcv$gdf))
  # below the integer range the seeds (hence the probes) are those of the
  # integer formula; above it the seed wraps into the range
  for (s in c(1L, 7L, 2023L, 21474L, -5L, -21474L)) {
    expect_identical(vapply(1:16, function(j) .ggcv_probe_seed(s, j), integer(1)),
                     s * 100003L + 1:16)
  }
  big <- vapply(1:16, function(j) .ggcv_probe_seed(20261002, j), integer(1))
  expect_false(anyNA(big))
  expect_identical(anyDuplicated(big), 0L)
  # invalid seeds are rejected up front
  for (bad in list(NA, Inf, c(1, 2), "a")) {
    expect_error(tt_gdf(dat$y_gauss, dat$X, lambda = 1, rank = 2L, k = 5L,
                        M = 2L, seed = bad, control = ctl),
                 "`seed` must be one finite number")
  }
  # a failed probe draw stops the batch, also with forked workers (it used
  # to come back as NA contributions from the workers)
  skip_on_os("windows")
  skip_on_cran()
  st <- .ge_state(dat$y_gauss, dat$X, rank = 2L, k = 5L, control = ctl)
  st$seed <- NA_real_
  st$n_cores <- 2L
  expect_error(.ggcv_eval_batch(matrix(0, 1L, 1L), 2L, st, "gdf"), "probe seed")
})

test_that("Gaussian UBRE decisions do not depend on the units of y", {
  # UBRE is compared on the AIC scale, D / scale + 2 GDF, so y and 10 * y with
  # the scale times 100 give the same scores, the same pattern-search moves
  # and the same final choice (same probes). With D + 2 * scale * GDF, tol
  # was in units of y^2: here the anisotropic case came out isotropic at y
  # and grouped at 10 * y.
  set.seed(21)
  n <- 200
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  e <- stats::rnorm(n, sd = 0.2)
  ys <- list(iso = 1 + sin(2 * pi * X[, 1]) + cos(pi * X[, 2]) + e,
             aniso = 1 + sin(2 * pi * X[, 1]) + X[, 2] + e)
  # fixed budget, so that no tolerance stop can flip by roundoff between the
  # two scales; the scores then agree to roundoff (~1e-7 relative here)
  ctl <- tt_control(max_sweeps = 8L, seed = 1L, compute_edf = FALSE,
                    warn_lambda_boundary = FALSE)
  sel <- function(y, units, groups) {
    .ge_quiet_budget(tt_ggcv(units * y, X, rank = 6L, k = 6L, groups = groups,
                             theta_lower = -3, theta_upper = 3,
                             criterion = "ubre", scale = (units * 0.2)^2,
                             budget = "fixed", control = ctl))
  }
  for (case in c("iso", "aniso")) {
    groups <- if (case == "aniso") 1:2 else NULL
    a <- sel(ys[[case]], 1, groups)
    b <- sel(ys[[case]], 10, groups)
    th <- grep("^theta", names(a$search), value = TRUE)
    expect_equal(b$search[, th], a$search[, th], tolerance = 1e-6)
    expect_equal(b$search$score_dev, a$search$score_dev, tolerance = 1e-6)
    expect_equal(b$theta, a$theta, tolerance = 1e-6)
    expect_identical(b$decision, a$decision)
    expect_identical(b$paired$qualifies, a$paired$qualifies)
    expect_identical(b$paired$winner, a$paired$winner)
    # the ubre / score fields keep their definition (units of y^2)
    expect_equal(b$ubre, 100 * a$ubre, tolerance = 1e-6)
    if (case == "aniso") expect_identical(a$decision, "grouped")
  }
})

test_that("Gaussian GDF does not depend on the units of y; the Poisson step is unchanged", {
  # Gaussian step epsilon_rel * RMS(y), without the old floor at 1: the ALS
  # map is scale-equivariant, so y and 1e-3 * y give the same GDF. With the
  # floor, the step was 1.4 x RMS at 1e-3 * y and the secant of the low-rank
  # map overstated the GDF.
  set.seed(5)
  n <- 200
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  y <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.5 * X[, 2] +
    stats::rnorm(n, sd = 0.3)
  ctl <- tt_control(max_sweeps = 30L, seed = 1L, compute_edf = FALSE)
  gdf_of <- function(yy) {
    tt_gdf(yy, X, lambda = 10^-1.5, rank = 2L, k = 8L, M = 8L,
           budget = "fixed", control = ctl)
  }
  a <- gdf_of(y)
  b <- gdf_of(1e-3 * y)
  expect_equal(a$epsilon, 1e-3 * sqrt(mean(y^2)), tolerance = 1e-12)
  expect_equal(b$epsilon, 1e-3 * a$epsilon, tolerance = 1e-12)
  expect_equal(b$gdf, a$gdf, tolerance = 1e-6)
  expect_equal(b$contrib, a$contrib, tolerance = 1e-6)
  # y = 0: the step is epsilon_rel itself
  expect_identical(.ge_state(numeric(n), X, rank = 2L, k = 8L, control = ctl)$eps,
                   1e-3)
  # counts keep epsilon_rel * max(RMS, 1)
  for (mu in c(0.3, 5)) {
    yp <- stats::rpois(n, mu)
    st <- .ge_state(yp, X, family = stats::poisson(), rank = 2L, k = 8L,
                    control = ctl)
    expect_identical(st$eps, 1e-3 * max(sqrt(sum(yp^2) / n), 1))
  }
})

# ---------------------------------------------------------------------------
# the adaptive pilot, probe validation, unstable evaluations, groups, errors
# ---------------------------------------------------------------------------

test_that("adaptive pilot: objective rule, B = ceiling(1.5 n_star) + 2, fixed-B estimator", {
  skip_if_not_installed("testthat", "3.1.7")  # local_mocked_bindings()
  # Review case: the RSS stop of ttps() (tol = 1e-7, floored at 1) ends this
  # low-rank fit after 13 sweeps at a turning point of the RSS, while the
  # penalized objective still changes by 2.5e-5 per sweep (after 3 sweeps
  # when y is in units of 1e-3).
  set.seed(3)
  n <- 400
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  y <- sin(3 * pi * X[, 1]) * cos(3 * pi * X[, 2]) + X[, 1] +
    stats::rnorm(n, sd = 0.1)
  fit_at <- function(yy, sweeps, tol = 0) {
    ttps(yy, X, rank = 3L, k = 10L, lambda = 1e-3, optimizer = "ALS",
         control = tt_control(max_sweeps = sweeps, tol = tol, seed = 1L,
                              compute_edf = FALSE, warn_lambda_boundary = FALSE))
  }
  expect_identical(fit_at(y, 400L, tol = 1e-7)$n_sweeps, 13L)
  # n_star from the objective path of a plain cold run: the first sweep at
  # which the relative change has been <= fit_tol twice in a row
  v <- vapply(fit_at(y, 80L)$history, function(h) h$objective, numeric(1))
  ok <- abs(diff(v)) / abs(v[-length(v)]) <= 1e-7
  n_star <- as.integer(which(ok[-1L] & ok[-length(ok)])[1L] + 2L)
  expect_gt(n_star, 13L)
  B <- as.integer(ceiling(1.5 * n_star) + 2)
  # record every fit of the evaluation: response, budget and tol
  fits <- list()
  map_scattered <- .ggcv_map_scattered
  local_mocked_bindings(.ggcv_map_scattered = function(...) {
    map <- map_scattered(...)
    fit <- map$fit
    map$fit <- function(yvec, lambda, ctrl, init = NULL) {
      fits[[length(fits) + 1L]] <<- list(at_y = identical(yvec, map$y),
                                         budget = as.integer(ctrl$max_sweeps),
                                         tol = ctrl$tol)
      fit(yvec, lambda, ctrl, init)
    }
    map
  })
  g <- tt_gdf(y, X, lambda = 1e-3, rank = 3L, k = 10L, M = 2L)
  expect_true(g$converged)
  expect_identical(g$n_iter, n_star)
  expect_identical(g$budget, B)
  # every fit runs cold with tol = 0: the pilot at 20, 40 and 80 sweeps
  # (doubling until a run reaches n_star), then the estimator, the reference
  # and the probes at B
  at_y <- vapply(fits, function(f) f$at_y, logical(1))
  budgets <- vapply(fits, function(f) f$budget, integer(1))
  expect_true(all(vapply(fits, function(f) f$tol == 0, logical(1))))
  expect_identical(budgets[at_y], c(20L, 40L, 80L, B))
  expect_identical(budgets[!at_y], rep(B, 3L))
  expect_identical(g$n_fits, length(fits))
  # the estimator is that map: the plain cold fit at B
  expect_identical(fitted(g$fit), fitted(fit_at(y, B)))
  # no floor: the same n_star, B and GDF for y in units of 1e-3
  g3 <- tt_gdf(1e-3 * y, X, lambda = 1e-3, rank = 3L, k = 10L, M = 2L)
  expect_identical(g3$n_iter, n_star)
  expect_identical(g3$budget, B)
  expect_equal(g3$gdf, g$gdf, tolerance = 1e-6)

  # Poisson PIRLS with a cap that binds (the paired-SE test's grid): the
  # objective drifts by about 3e-7 per iteration, so the pilot never settles
  # and the evaluation runs the cap unconverged, although the deviance rule
  # of fixed mode (last relative change <= 1e-5) would call it converged
  gr <- .ge_grid()
  set.seed(21)
  yp <- stats::rpois(nrow(gr$X), exp(log(0.8) + sin(2 * pi * gr$X[, 1]) +
                                       cos(2 * pi * gr$X[, 2])))
  g0 <- tt_gdf(yp, gr$X, lambda = 1, family = stats::poisson(), rank = 3L,
               k = 6L, M = 4L,
               control = tt_control(pirls_maxit = 20L, seed = 1L, compute_edf = FALSE))
  expect_identical(g0$n_iter, 20L)
  expect_identical(g0$budget, 20L)
  expect_false(g0$converged)
  expect_lt(g0$last_rel_change, 1e-5)
})

test_that("probe validation re-evaluates an unstable evaluation at 2B, then 4B", {
  skip_if_not_installed("testthat", "3.1.7")  # local_mocked_bindings()
  # Low-rank ALS caught mid-transition: at a fixed 60 sweeps the probe fits
  # at lambda = 1 are not converged and the contributions come out as
  # (-0.6, -18.9, 53.8, 24.4) (per-probe CV 2.1); at 120 sweeps they settle
  # at (6.2, 10.8, 12.1, 6.2), within 0.3 of a 240-sweep run.
  set.seed(11)
  n <- 200
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  y <- 0.4 + sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + stats::rnorm(n, sd = 0.3)
  gdf_at <- function(sweeps) {
    tt_gdf(y, X, lambda = 1, rank = 2L, k = 8L, M = 4L, budget = "fixed",
           control = tt_control(max_sweeps = sweeps, seed = 1L, compute_edf = FALSE))
  }
  a <- gdf_at(60L)
  b <- gdf_at(120L)
  expect_identical(a$reevals, 1L)
  expect_identical(a$budget, 120L)
  expect_true(a$stable)
  expect_identical(b$reevals, 0L)
  # the whole evaluation (estimator, reference, probes) is the one at 2B
  expect_identical(a$contrib, b$contrib)
  expect_identical(a$deviance, b$deviance)
  expect_identical(fitted(a$fit), fitted(b$fit))
  expect_identical(a$n_fits, 2L * (4L + 2L))
  expect_lt(b$gdf_cv, 1)

  # a probe contribution that stays wrong at every budget (probe 2's fits
  # are shifted): two re-evaluations, then flagged unstable and returned
  r2 <- .ge_probes(n, 2L)[, 2L]
  eps <- 1e-3 * sqrt(mean(y^2))
  y2 <- y + eps * (r2 + 1) / 2
  map_scattered <- .ggcv_map_scattered
  local_mocked_bindings(.ggcv_map_scattered = function(...) {
    map <- map_scattered(...)
    fit <- map$fit
    map$fit <- function(yvec, lambda, ctrl, init = NULL) {
      f <- fit(yvec, lambda, ctrl, init)
      if (max(abs(yvec - y2)) < 1e-12) {
        f$fitted.values <- f$fitted.values + 50 * eps * (r2 > 0)
      }
      f
    }
    map
  })
  u <- gdf_at(120L)
  expect_identical(u$reevals, 2L)
  expect_identical(u$budget, 480L)
  expect_false(u$stable)
  expect_gt(abs(u$contrib[2L]), n)
  expect_identical(u$n_fits, 3L * (4L + 2L))
})

test_that("unstable evaluations are never moves, isotropic optima or eligible candidates", {
  skip_if_not_installed("testthat", "3.1.7")  # local_mocked_bindings()
  dat <- .ge_gauss_xy(61L, .ge_truth_aniso)
  select <- function() {
    .ge_warnings(tt_ggcv(
      dat$y, dat$X, rank = 6L, k = 6L, groups = 1:2, theta_lower = -3,
      theta_upper = 3, control = tt_control(seed = 1L, compute_edf = FALSE)
    ))
  }
  orig <- .ggcv_eval_batch
  ref <- select()$value
  expect_identical(ref$decision, "grouped")
  t_iso <- rep(ref$paired$theta1[ref$paired$candidate == "isotropic"], 2L)
  t_grp <- ref$theta
  at <- function(target) function(th) max(abs(th - target)) < 1e-9
  cand <- function(sel) cbind(sel$paired$theta1, sel$paired$theta2)
  # every point with log10(lambda_2) > 1.5 (the grouped optimum among them)
  # flagged unstable: the pattern search evaluates such points, in batches
  # and with doubling steps, but never moves there, so no final candidate
  # and not the winner lies there
  far <- function(th) th[2L] > 1.5
  local_mocked_bindings(.ggcv_eval_batch = .ge_mark_unstable(far, orig))
  sel <- select()$value
  srch <- sel$search[sel$search$stage == "pattern", ]
  expect_true(any(srch$theta2 > 1.5))
  expect_false(any(apply(cand(sel), 1L, far)))
  expect_identical(sel$paired$theta1[sel$paired$candidate == "isotropic"], t_iso[1L])
  expect_false(far(sel$theta))
  expect_true(all(sel$paired$stable))
  # the isotropic optimum flagged unstable: it is not the isotropic optimum
  # of the stage, nor any final candidate
  local_mocked_bindings(.ggcv_eval_batch = .ge_mark_unstable(at(t_iso), orig))
  sel <- select()$value
  expect_false(any(apply(cand(sel), 1L, at(t_iso))))
  expect_true(all(sel$paired$stable))
  # every evaluation unstable: no move is accepted, the choice falls back to
  # the unstable candidates and the winner's instability is reported
  local_mocked_bindings(.ggcv_eval_batch = .ge_mark_unstable(function(th) TRUE, orig))
  res <- select()
  expect_false(res$value$stable)
  expect_identical(res$value$decision, "isotropic")
  expect_false(any(res$value$paired$stable))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "unstable", fixed = TRUE)
})

test_that("probe validation: the CV limit grows at small GDF; outliers stay flagged", {
  # A symmetric smoother with eigenvalues in [0, 1] has a per-probe CV of at
  # most sqrt(2 / GDF); the limit is max(1, 3 * sqrt(2 / max(GDF, 1))).
  expect_equal(.ggcv_cv_limit(4), 3 * sqrt(0.5))
  expect_equal(.ggcv_cv_limit(0.2), 3 * sqrt(2))
  expect_identical(.ggcv_cv_limit(18), 1)
  expect_identical(.ggcv_cv_limit(400), 1)
  # Review case (finding 15): a converged, heavily smoothed evaluation with
  # four probes (contributions of the review's smooth-truth scan, GDF 4.44,
  # CV 1.10) was flagged by the limit 1 and could not be cleared
  cc <- c(1.9432, 11.7565, 1.7971, 2.2748)
  expect_gt(stats::sd(cc) / mean(cc), 1)
  expect_false(.ggcv_unstable(cc, n_eff = 200))
  # the same CV at GDF 20 is flagged: from GDF 18 up the limit is 1
  expect_true(.ggcv_unstable(c(9.5, 9.5, 9.5, 51.5), n_eff = 200))
  # one huge outlier is flagged: by the CV (large GDF, or 16 probes) ...
  expect_true(.ggcv_unstable(c(10, 11, 9, 160), n_eff = 200))
  expect_true(.ggcv_unstable(c(rep(4, 15L), 70), n_eff = 200))
  # ... or by the absolute checks; and the review's mid-transition case
  expect_true(.ggcv_unstable(c(5, 6, 4, 250), n_eff = 200))
  expect_true(.ggcv_unstable(c(5, 6, NA, 4), n_eff = 200))
  expect_true(.ggcv_unstable(c(-9, 1, 2, 1), n_eff = 200))
  expect_true(.ggcv_unstable(c(-0.58, -18.88, 53.84, 24.42), n_eff = 200))

  # End to end, the flagged point of that scan (seed 2, linear truth,
  # rank = k = 6, log10(lambda) = 4): its contributions equal the exact
  # r_j' S r_j of the dense smoother, so the flag was a false positive. It
  # cost two re-evaluations (B up to 44) and kept the point out of the
  # search; now it is stable at its own budget.
  set.seed(2)
  n <- 200
  X <- cbind(x1 = stats::runif(n), x2 = stats::runif(n))
  y <- 1 + X[, 1] - 0.5 * X[, 2] + stats::rnorm(n, sd = 0.3)
  g <- tt_gdf(y, X, lambda = 1e4, rank = 6L, k = 6L, M = 4L)
  expect_true(g$stable)
  expect_identical(g$reevals, 0L)
  expect_gt(g$gdf_cv, 1)
  S <- .ge_dense_smoother(X, 6L)(rep(1e4, 2L))
  R <- .ge_probes(n, 4L)
  expect_equal(g$contrib, colSums(R * (S %*% R)), tolerance = 1e-4)
})

test_that("Poisson smoothing groups resolve alike in every gGCV entry point", {
  g <- .ge_grid()
  set.seed(21)
  y <- stats::rpois(nrow(g$X), exp(log(0.8) + sin(2 * pi * g$X[, 1]) +
                                     cos(2 * pi * g$X[, 2])))
  groups_of <- function(...) {
    ctl <- tt_control(pirls_maxit = 4L, seed = 1L, compute_edf = FALSE,
                      warn_lambda_boundary = FALSE, ggcv_n_grid = 3L,
                      ggcv_n_refine = 0L, ggcv_M_search = 2L, ggcv_M_final = 2L,
                      ggcv_budget_check = FALSE, ggcv_refit = FALSE,
                      ggcv_pirls_maxit = 4L, ...)
    q <- function(expr) paste(.ge_quiet_budget(expr), collapse = "")
    c(q(tt_ggcv_poisson(y, g$X, rank = 2L, k = 5L, control = ctl)$groups),
      q(ttps(y, g$X, family = stats::poisson(), rank = 2L, k = 5L,
             lambda = "gGCV", control = ctl)$ggcv$groups),
      q(ttps(array(y, g$ng), axes = g$axes, array = TRUE,
             family = stats::poisson(), rank = 2L, k = 5L, lambda = "gGCV",
             control = ctl)$ggcv$groups))
  }
  # ggcv_groups first, then ggcv_poisson_anisotropic, on scattered rows, on
  # the array and in tt_ggcv_poisson()
  expect_identical(groups_of(ggcv_poisson_anisotropic = TRUE), rep("12", 3L))
  expect_identical(groups_of(ggcv_groups = c(1L, 1L), ggcv_poisson_anisotropic = TRUE),
                   rep("11", 3L))
  expect_identical(groups_of(ggcv_groups = 1:2, ggcv_poisson_anisotropic = FALSE),
                   rep("12", 3L))
  expect_identical(groups_of(), rep("11", 3L))
})

test_that("fit errors surface with their message", {
  dat <- .ge_scattered()
  ctl <- tt_control(max_sweeps = 10L, pirls_maxit = 5L, seed = 1L, compute_edf = FALSE)
  bad_knots <- list(c(0, 1), c(0, 1))
  # tt_gdf(): an error, not NA
  expect_error(tt_gdf(dat$y_gauss, dat$X, lambda = 1, rank = 2L, k = 5L,
                      knots = bad_knots, control = ctl),
               "'ord' must be positive integer", fixed = TRUE)
  expect_error(tt_gdf(dat$y_gauss, dat$X, lambda = 1, rank = 2L, k = 5L,
                      knots = bad_knots, probes = "unit", control = ctl),
               "'ord' must be positive integer", fixed = TRUE)
  g <- .ge_grid()
  expect_error(tt_gdf_array(array(dat$y_pois[seq_len(72L)], g$ng), lambda = 1,
                            family = stats::poisson(), rank = 2L, k = 2L,
                            control = ctl),
               "k > degree", fixed = TRUE)
  # selections that end with no finite criterion report the first error
  expect_error(tt_ggcv(dat$y_gauss, dat$X, rank = 2L, k = 5L, knots = bad_knots,
                       n_grid = 3L, n_refine = 0L, M_search = 2L, M_final = 2L,
                       control = ctl),
               "first fit error: 'ord' must be positive integer", fixed = TRUE)
  expect_error(ttps(dat$y_pois, dat$X, family = stats::poisson(), rank = 2L,
                    k = 5L, lambda = "gGCV", cyclic = c(TRUE, FALSE),
                    period = c(1, NA), control = ctl),
               "`period` must be a list of length d", fixed = TRUE)
})
