# Unified TT-gGCV engine: one evaluator and one optimizer for every gGCV
# route (array grids, scattered rows, Berman-Turner rows).
#
# Fixed-lambda map. A map (.ggcv_map_array(), .ggcv_map_scattered()) wraps
# the fitting map y |-> mu_hat(y; lambda) of ttps() at a fixed smoothing
# vector: same rows, basis, rank, offset and weights for every response y.
#
# Criterion. With D the deviance at y and GDF = tr(d mu_hat / d y) over the
# n_eff rows carrying positive weight,
#   GCV  = n_eff * D / (n_eff - GDF)^2                   (scale estimated),
#   UBRE = D / n_eff + 2 * scale * GDF / n_eff - scale   (scale known).
# criterion = "auto" takes UBRE for Poisson (scale 1) and GCV for Gaussian,
# as mgcv's GCV.Cp. Every comparison uses the deviance-unit form
#   score_dev = D / scale + 2 * GDF (UBRE, the AIC scale)  or  n_eff * log(GCV),
# whose Monte Carlo noise is g * (GDF noise), g = d score_dev / d GDF; both
# forms are free of the units of y, so tol means the same for any data.
#
# GDF probes. Non-negative Rademacher perturbations against one shared
# reference: r_j ~ Rademacher on the support, probe fit at
# y + eps * (r_j + 1) / 2, reference fit at y + eps / 2, and
#   c_j = (2 / eps) * sum(r_j * (mu_j - mu_ref)),   E[c_j] ~= GDF.
# Responses only move up, so Poisson counts stay valid at y = 0 (no
# clipping), and the reference's own convergence error cancels. The step is
# eps = epsilon_rel * RMS(y) for Gaussian data (the ALS map is
# scale-equivariant) and epsilon_rel * max(RMS(y), 1) for counts. Probe j is
# drawn from seed * 100003 + j (wrapped into the integer range when it
# leaves it): probes are common random numbers across lambda and nested in
# M, so the first M_search probes of a final-stage evaluation are the search
# probes and the cache fits only the new ones. The probes of a batch are
# drawn once, in the parent process, before the fits are dispatched. A user
# probe bank (tt_gdf(), tt_gdf_array(): a matrix of +-1 values) replaces the
# seeded draws column by column.
#
# Iteration budget. All fits of a point differentiate one deterministic
# map: the estimator at y, the reference and the probes run a fixed budget B
# (control$tol = 0, so PIRLS / ALS never stop early) from the same cold
# initialization (control$seed). budget = "adaptive" (default): a pilot at y
# fixes B. It runs cold with tol = 0 for b0, 2 b0, 4 b0, ... iterations
# (b0 = 20 ALS sweeps for Gaussian data, 10 PIRLS iterations for Poisson; the
# last run is the cap, pirls_maxit / max_sweeps) until the relative change
# of the penalized objective, |dObj| / |Obj| from the fit history, has been
# at most fit_tol at two consecutive iterations; n_star is the first
# iteration where that holds. Then converged = (n_star exists) and
# B = min(cap, ceiling(1.5 * n_star) + 2) (B = cap when not converged): the
# GDF settles more slowly than the objective. The rule does not go through
# control$tol: the stop of ttps() tests the RSS (which can stall at a
# turning point while the fit still moves) against max(1, RSS), an absolute
# rule below RSS = 1, so B would depend on the units of y. Cold runs with
# tol = 0 are prefixes of each other, so the doubling only re-runs the same
# map further; the estimator is the cold fit at B (the last pilot run when
# it ran exactly B), so the deviance, the returned fit and the probes share
# one map. A budget that is right at one lambda is wrong at another (slowly
# mixing low-rank ALS needs far more sweeps at small lambda), so B follows
# lambda. budget = "fixed": B is the control budget for every fit, and an
# evaluation counts as converged when the relative change of its estimator's
# last iteration is <= conv_tol.
#
# Probe validation. A pilot that settles at y does not certify the
# neighbourhood that the finite differences probe: probe responses can need
# more iterations from the same cold start. An evaluation whose
# contributions contain a non-finite value or an outlier (|c_j| > n_eff, or
# a mean <= 0, or gdf_cv = sd / mean above max(1, 3 * sqrt(2 / max(GDF, 1))),
# see .ggcv_cv_limit()) is re-evaluated (estimator, reference and all its
# probes, same probes) at 2B and then 4B; if it stays so it is flagged
# unstable (stable = FALSE) and treated like an unconverged one: it is never
# an accepted pattern-search move, the isotropic optimum of a stage (unless
# no stable isotropic point exists) or an eligible final candidate (unless
# nothing else is eligible, with a warning).
#
# Search (.ggcv_search()). (1) Isotropic stage, always: a grid on the shared
# log10(lambda), n_refine halvings of the grid step around the best point
# and one parabolic step. (2) Grouped stage, several groups only: optional
# starts (theta_start, Sobol points), then a coordinate pattern search with
# step 0.5 and then 0.25, extending accepted moves with doubling steps.
# (3) Final stage: the isotropic and grouped optima, the nearest stable
# isotropic evaluation on the smoother side (larger theta) and the next best
# distinct stable points are re-scored with M_final probes.
# Paired rule: b improves on a only if score_dev(a) - score_dev(b) exceeds
# max(tol, 2 * SE), where SE is the standard error of the paired difference
# over the common probes; common random numbers make it far smaller than the
# unpaired SE. Final choice (.ggcv_choose()): among the converged and stable
# candidates, the most regular one that the minimum does not improve on by
# the paired rule wins (isotropic before grouped, larger theta first), as
# in the one-standard-error rule. On Dette d = 8 (Berman-Turner rows) the
# plain argmin sat on the grid boundary with a Monte Carlo GDF of
# 432 +- 334, against 145 +- 10 at the lambda chosen by thinning CV.
# (4) Budget check: the winner's GDF is re-estimated at twice its budget
# with the same probes. When it moves by more than budget_tol, the
# candidate is rescored with the less favourable of its B and 2B scores and
# the choice is redone with that conservative score (budget_fallback): if it
# still wins it is kept (budget_verified = FALSE, with a warning that gives
# the GDF change and the margin), otherwise the new winner gets its own
# check (at most three checks). A GDF change matters only when it could
# reverse the choice: dropping such a candidate outright once returned a
# model about 19 score units worse by the exact criterion, over a GDF change
# worth about 1 unit. The first winner is always checked (the diagnostic of
# the returned fit); no fallback check is spent on an unconverged or
# unstable candidate (the check could not verify its map), so an unconverged
# first winner gets one diagnostic check and no fallback. The returned
# candidate is always a checked one: when the redone choice would return an
# unconverged or unstable candidate (e.g. the failing winner's score at 2B
# is not finite and no other converged, stable candidate is left), or a
# candidate left unchecked because the limit of checks was reached, the
# first winner is kept (budget_verified = FALSE, with a warning).

# ---------------------------------------------------------------------------
# maps
# ---------------------------------------------------------------------------

#' Fixed-lambda map of an array (product-grid) TT fit.
#' @keywords internal
#' @noRd
.ggcv_map_array <- function(Y, axes, family, rank, k, degree, penalty_order,
                            cyclic, period, knots, weights, offset) {
  if (!is.array(Y) || length(dim(Y)) < 2L) {
    stop("`Y` must be a d-way array (d >= 2).", call. = FALSE)
  }
  dims <- dim(Y)
  d <- length(dims)
  n_cells <- prod(dims)
  fam <- normalize_family(family)
  key <- family_key(fam)
  if (!key %in% c("gaussian", "poisson")) {
    stop("tt_ggcv_array supports gaussian() and poisson() only.", call. = FALSE)
  }
  y <- as.numeric(Y)
  if (anyNA(y)) stop("`Y` must not contain NA.", call. = FALSE)
  if (identical(key, "poisson") && any(y < 0)) {
    stop("poisson requires non-negative Y.", call. = FALSE)
  }
  w <- if (is.null(weights)) NULL else as.numeric(weights)
  if (!is.null(w)) {
    if (length(w) != n_cells) stop("`weights` must have one value per cell.", call. = FALSE)
    if (any(!is.finite(w)) || any(w < 0)) {
      stop("`weights` must be finite and non-negative.", call. = FALSE)
    }
    if (all(abs(w - 1) < 1e-12)) w <- NULL
  }
  support <- if (is.null(w)) rep(TRUE, n_cells) else w > 0
  n_eff <- sum(support)
  if (n_eff < 2L) stop("Need at least two cells with positive weight.", call. = FALSE)
  off <- if (is.null(offset)) NULL else as.numeric(offset)
  if (!is.null(off) && length(off) != n_cells && length(off) != 1L) {
    stop("`offset` must be scalar or have one value per cell.", call. = FALSE)
  }
  force(axes); force(rank); force(k); force(degree); force(penalty_order)
  force(cyclic); force(period); force(knots)
  fit <- function(yvec, lambda, ctrl, init = NULL) {
    # ttps() reseeds the global RNG for its core initialization; keep the
    # caller's stream intact.
    .tt_with_preserved_seed(ttps(
      array(yvec, dim = dims), axes = axes, array = TRUE, family = fam,
      rank = rank, k = k, degree = degree, penalty_order = penalty_order,
      lambda = lambda, optimizer = "ALS", init = init, control = ctrl,
      knots = knots, offset = off, weights = w, cyclic = cyclic,
      period = period
    ))
  }
  list(mode = "array", y = y, support = support, n_eff = n_eff, weights = w,
       fam = fam, key = key, d = d, n_rows = n_cells, dims = dims,
       npar = .ggcv_npar(d, rank, k), fit = fit)
}

#' Fixed-lambda map of a scattered-row TT fit (Gaussian or Poisson).
#'
#' `optimizer = "ALS"` gives ALS for Gaussian and PIRLS-ALS for Poisson, the
#' same fitters as the array map. Rows enter with unit weight (scattered gGCV
#' rejects non-uniform weights), so every row is on the probe support.
#' @keywords internal
#' @noRd
.ggcv_map_scattered <- function(y, X, family, rank, k, degree, penalty_order,
                                cyclic, period, knots, offset) {
  fam <- normalize_family(family)
  key <- family_key(fam)
  if (!key %in% c("gaussian", "poisson")) {
    stop("tt_ggcv supports gaussian() and poisson() only.", call. = FALSE)
  }
  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y)
  if (nrow(X) != n) stop("nrow(X) must equal length(y).", call. = FALSE)
  if (ncol(X) < 2L) stop("Need at least d = 2 covariates.", call. = FALSE)
  if (anyNA(y) || anyNA(X)) stop("`y` and `X` must not contain NA.", call. = FALSE)
  if (identical(key, "poisson") && any(y < 0)) {
    stop("poisson requires non-negative y.", call. = FALSE)
  }
  if (n < 2L) stop("Need at least two observations.", call. = FALSE)
  off <- if (is.null(offset)) NULL else as.numeric(offset)
  if (!is.null(off) && length(off) != n && length(off) != 1L) {
    stop("`offset` must be scalar or have one value per row.", call. = FALSE)
  }
  force(rank); force(k); force(degree); force(penalty_order)
  force(cyclic); force(period); force(knots)
  fit <- function(yvec, lambda, ctrl, init = NULL) {
    .tt_with_preserved_seed(ttps(
      y = yvec, X = X, family = fam, rank = rank, k = k, degree = degree,
      penalty_order = penalty_order, lambda = lambda, optimizer = "ALS",
      cyclic = cyclic, period = period, knots = knots, offset = off,
      init = init, control = ctrl
    ))
  }
  list(mode = "scattered", y = y, support = rep(TRUE, n), n_eff = n,
       weights = NULL, fam = fam, key = key, d = ncol(X), n_rows = n,
       npar = .ggcv_npar(ncol(X), rank, k), fit = fit)
}

#' Number of TT parameters (cores plus intercept) of a map's model.
#' @keywords internal
#' @noRd
.ggcv_npar <- function(d, rank, k) {
  r <- rep(as.integer(rank), length.out = max(1L, d - 1L))
  ranks <- c(1L, if (d > 1L) r[seq_len(d - 1L)] else integer(0), 1L)
  K <- rep(as.integer(k), length.out = d)
  sum(ranks[-(d + 1L)] * K * ranks[-1L]) + 1L
}

# ---------------------------------------------------------------------------
# setup, probes, convergence
# ---------------------------------------------------------------------------

#' Selector state: groups, criterion, the fixed-budget control of the
#' estimator, reference and probe fits, the first pilot run of the adaptive
#' budget (`pilot_b0`), probe step and probe bank, worker count, evaluation
#' cache and search log.
#' @keywords internal
#' @noRd
.ggcv_setup <- function(map, groups, criterion, scale, probe_init,
                        probe_budget, epsilon_rel, n_cores, seed, control,
                        conv_tol = 1e-5, budget = "adaptive", fit_tol = 1e-7,
                        probes = NULL, gdf_method = "mc",
                        exact_max_npar = 2000L, exact_max_cells = 5e7) {
  d <- map$d
  key <- map$key
  groups <- if (is.null(groups)) rep(1L, d) else as.integer(groups)
  if (length(groups) != d || anyNA(groups)) {
    stop("`groups` must be an integer vector of length d.", call. = FALSE)
  }
  groups <- match(groups, sort(unique(groups)))
  criterion <- match.arg(criterion, c("auto", "gcv", "ubre"))
  if (identical(criterion, "auto")) {
    criterion <- if (identical(key, "poisson")) "ubre" else "gcv"
  }
  if (!is.null(scale) && length(scale) == 1L && is.na(scale)) scale <- NULL
  if (is.null(scale) && identical(key, "poisson")) scale <- 1
  if (identical(criterion, "ubre") && is.null(scale)) {
    stop("criterion = 'ubre' needs a known `scale` for gaussian().", call. = FALSE)
  }
  if (!is.null(scale)) {
    scale <- as.numeric(scale)
    if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
      stop("`scale` must be one positive number.", call. = FALSE)
    }
  }
  probe_init <- match.arg(probe_init, c("cold", "warm"))
  budget <- match.arg(budget, c("adaptive", "fixed"))
  fit_tol <- as.numeric(fit_tol)
  if (length(fit_tol) != 1L || !is.finite(fit_tol) || fit_tol <= 0) {
    stop("`fit_tol` must be a positive number.", call. = FALSE)
  }
  epsilon_rel <- as.numeric(epsilon_rel)
  if (length(epsilon_rel) != 1L || !is.finite(epsilon_rel) || epsilon_rel <= 0) {
    stop("`epsilon_rel` must be one positive number.", call. = FALSE)
  }
  seed <- suppressWarnings(as.numeric(seed))
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("`seed` must be one finite number.", call. = FALSE)
  }
  ctrl <- control
  if (!inherits(ctrl, "tt_control")) ctrl <- do.call(tt_control, as.list(ctrl))
  ctrl$tol <- 0
  ctrl$compute_edf <- FALSE
  ctrl$trace <- FALSE
  ctrl$monitor <- FALSE
  ctrl$warn_lambda_boundary <- FALSE
  # the iteration count that the budget controls, and its cap
  iter_name <- if (identical(key, "poisson")) "pirls_maxit" else "max_sweeps"
  cap <- as.integer(ctrl[[iter_name]])
  if (length(cap) != 1L || is.na(cap) || cap < 1L) {
    stop(sprintf("control$%s must be a positive integer.", iter_name), call. = FALSE)
  }
  ctrl_probe <- ctrl
  if (identical(probe_init, "warm")) {
    if (identical(key, "poisson")) {
      ctrl_probe$pirls_maxit <- as.integer(probe_budget %||%
                                             max(2L, ceiling(ctrl$pirls_maxit / 3)))
    } else {
      ctrl_probe$max_sweeps <- as.integer(probe_budget %||%
                                            max(2L, ceiling(ctrl$max_sweeps / 3)))
    }
  }
  # finite-difference step: relative to the RMS of y on the support for
  # Gaussian data, with no floor (the ALS map is scale-equivariant, so the
  # GDF does not depend on the units of y; epsilon_rel itself when y = 0);
  # counts keep the floor max(RMS, 1)
  y <- map$y
  rms <- sqrt(sum(y[map$support]^2) / map$n_eff)
  has_rms <- is.finite(rms) && rms > 0
  eps <- epsilon_rel * if (identical(key, "poisson")) {
    if (has_rms) max(rms, 1) else 1
  } else if (has_rms) {
    rms
  } else {
    1
  }
  n_cores <- max(1L, as.integer(n_cores))
  if (n_cores > 1L && (.Platform$OS.type == "windows" ||
                       !requireNamespace("parallel", quietly = TRUE))) {
    n_cores <- 1L
  }
  # GDF method: "exact" = implicit differentiation with the full Hessian
  # (.tt_exact_edf_core(); one fit per lambda polished to a stationary point
  # by Newton, no probes),
  # "mc" = Monte Carlo finite differences of the fixed-budget map. "auto" =
  # exact when the dense Hessian / Jacobian are affordable and no probe bank
  # is given.
  gdf_method <- match.arg(gdf_method, c("auto", "exact", "mc"))
  if (identical(gdf_method, "auto")) {
    gdf_method <- if (is.null(probes) &&
                      map$npar <= as.integer(exact_max_npar) &&
                      as.numeric(map$n_rows) * map$npar <= exact_max_cells) {
      "exact"
    } else {
      "mc"
    }
  }
  ctrl_exact <- control
  if (!inherits(ctrl_exact, "tt_control")) ctrl_exact <- do.call(tt_control, as.list(ctrl_exact))
  ctrl_exact$tol <- fit_tol
  ctrl_exact$compute_edf <- FALSE
  ctrl_exact$trace <- FALSE
  ctrl_exact$monitor <- FALSE
  ctrl_exact$warn_lambda_boundary <- FALSE
  log <- new.env(parent = emptyenv())
  log$rows <- list()
  log$first_error <- NA_character_  # first failed fit of the run
  array_mode <- identical(map$mode, "array")
  list(
    map = map, y = y, support = map$support, n_eff = map$n_eff,
    n_rows = map$n_rows, fam = map$fam, key = key, d = d,
    groups = groups, n_groups = max(groups), criterion = criterion,
    scale = scale, probe_init = probe_init, eps = eps, n_cores = n_cores,
    seed = trunc(seed), control = ctrl, control_probe = ctrl_probe,
    control_user = control, conv_tol = as.numeric(conv_tol), budget = budget,
    gdf_method = gdf_method, control_exact = ctrl_exact,
    fit_tol = fit_tol, iter_name = iter_name, cap = cap,
    # first pilot run of the adaptive budget (doubled up to the cap)
    pilot_b0 = if (identical(key, "poisson")) 10L else 20L,
    probe_budget = if (is.null(probe_budget)) NULL else as.integer(probe_budget),
    probe_bank = if (is.null(probes)) NULL else .ggcv_probe_bank(probes, map),
    label = if (array_mode) "tt_ggcv_array" else "tt_ggcv",
    prefix = if (array_mode) "gGCV-array" else "gGCV",
    verbose = FALSE, cache = new.env(parent = emptyenv()), log = log
  )
}

#' Validate a user probe bank: one row per observation (array: per cell, in
#' `as.numeric(Y)` order), one column per probe, -1 / +1 on the support.
#' Rows off the support are ignored (set to 0).
#' @keywords internal
#' @noRd
.ggcv_probe_bank <- function(probes, map) {
  P <- if (is.data.frame(probes)) as.matrix(probes) else probes
  if (!is.matrix(P) || !is.numeric(P)) {
    stop("`probes` must be \"rademacher\", \"unit\" or a numeric matrix of ",
         "-1 / +1 values.", call. = FALSE)
  }
  what <- if (identical(map$mode, "array")) "cell" else "observation"
  if (nrow(P) != map$n_rows) {
    stop(sprintf("`probes` must have one row per %s (%d), not %d.", what,
                 map$n_rows, nrow(P)), call. = FALSE)
  }
  if (ncol(P) < 1L) stop("`probes` must have at least one column.", call. = FALSE)
  on <- P[map$support, , drop = FALSE]
  if (anyNA(on) || any(on != 1 & on != -1)) {
    stop("`probes` must contain only -1 and +1 on the rows that enter the fit.",
         call. = FALSE)
  }
  bank <- matrix(0, nrow(P), ncol(P))
  bank[map$support, ] <- on
  bank
}

#' Probe j on the support: column j of the user bank, else the seeded
#' Rademacher draw (common random numbers, nested in M).
#' @keywords internal
#' @noRd
.ggcv_probe <- function(j, st) {
  if (!is.null(st$probe_bank)) return(st$probe_bank[, as.integer(j)])
  .tt_with_preserved_seed({
    set.seed(.ggcv_probe_seed(st$seed, j))
    r <- numeric(st$n_rows)
    r[st$support] <- sample(c(-1, 1), st$n_eff, replace = TRUE)
  })
  r
}

#' RNG seed of probe j: seed * 100003 + j in double precision, wrapped into
#' the integer range only when it leaves it, so every seed for which the
#' integer formula fits keeps its probes (seeds above 21474 used to overflow
#' to NA).
#' @keywords internal
#' @noRd
.ggcv_probe_seed <- function(seed, j) {
  s <- as.numeric(seed) * 100003 + as.numeric(j)
  if (length(s) != 1L || !is.finite(s)) {
    stop("Invalid probe seed (the seed must be one finite number).", call. = FALSE)
  }
  if (abs(s) > .Machine$integer.max) s <- s %% .Machine$integer.max
  as.integer(s)
}

#' Message of a failed fit job (the error kept by fit_safe(), or the
#' try-error of a forked worker); `NA` when the job returned a value.
#' @keywords internal
#' @noRd
.ggcv_fit_error <- function(x) {
  if (inherits(x, "ggcv_fit_error")) return(x$message)
  if (inherits(x, "try-error")) {
    cond <- attr(x, "condition")
    return(if (inherits(cond, "condition")) conditionMessage(cond) else as.character(x))
  }
  NA_character_
}

#' " (first fit error: ...)" for error messages when a fit of this run
#' failed, else "".
#' @keywords internal
#' @noRd
.ggcv_error_note <- function(st) {
  msg <- st$log$first_error
  if (is.null(msg) || is.na(msg)) "" else sprintf(" (first fit error: %s)", msg)
}

#' A number for a message: `sprintf(fmt, x)`, "not finite" (Inf, NaN) or
#' "not available" (NA), never "NA" or "Inf".
#' @keywords internal
#' @noRd
.ggcv_num_txt <- function(x, fmt = "%.2g") {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L) return("not available")
  if (is.finite(x)) return(sprintf(fmt, x))
  if (is.nan(x) || is.infinite(x)) "not finite" else "not available"
}

#' Relative GDF change of a budget check for a message: "+1.8%", or "not
#' finite" when the GDF at either budget is not finite (never "NA%").
#' @keywords internal
#' @noRd
.ggcv_pct_txt <- function(rel) {
  rel <- suppressWarnings(as.numeric(rel))
  if (length(rel) == 1L && is.finite(rel)) sprintf("%+.1f%%", 100 * rel) else "not finite"
}

#' "GDF changes by +1.8% when <what> doubles", or "the GDF change is not
#' finite when <what> doubles".
#' @keywords internal
#' @noRd
.ggcv_change_txt <- function(rel, what) {
  pct <- .ggcv_pct_txt(rel)
  if (identical(pct, "not finite")) {
    sprintf("the GDF change is not finite when %s doubles", what)
  } else {
    sprintf("GDF changes by %s when %s doubles", pct, what)
  }
}

#' Relative change of the last iteration of a fit, from its history:
#' PIRLS (data frame: deviance, else objective) or ALS (list of sweeps:
#' objective, else rss). `NA` when no history is recorded.
#' @keywords internal
#' @noRd
.ggcv_last_rel_change <- function(fit) {
  h <- fit$history %||% fit$convergence$history
  v <- NULL
  if (is.data.frame(h)) {
    col <- intersect(c("deviance", "objective"), names(h))
    if (length(col)) v <- as.numeric(h[[col[1L]]])
  } else if (is.list(h) && length(h) && is.list(h[[1L]])) {
    nm <- if (!is.null(h[[1L]]$objective)) "objective" else "rss"
    v <- vapply(h, function(z) as.numeric(z[[nm]] %||% NA_real_)[1L], numeric(1))
  }
  if (length(v) < 2L) return(NA_real_)
  a <- v[length(v) - 1L]
  b <- v[length(v)]
  if (!is.finite(a) || !is.finite(b)) return(Inf)
  abs(b - a) / max(1, abs(a))
}

#' Iterations a fit used: PIRLS iterations (Poisson) or ALS sweeps
#' (Gaussian); `NA` when the fit does not record them.
#' @keywords internal
#' @noRd
.ggcv_n_iter <- function(fit, st) {
  v <- if (identical(st$key, "poisson")) fit$n_pirls else fit$n_sweeps
  v <- suppressWarnings(as.integer(v))
  if (length(v) && is.finite(v[1L])) v[1L] else NA_integer_
}

#' Penalized objective after each iteration of a fit, from its history:
#' PIRLS (data frame, column `objective`) or ALS (list of sweeps, field
#' `objective`); `numeric(0)` when no history is recorded.
#' @keywords internal
#' @noRd
.ggcv_objective_path <- function(fit) {
  h <- fit$history %||% fit$convergence$history
  if (is.data.frame(h)) {
    col <- intersect(c("objective", "deviance"), names(h))
    return(if (length(col)) as.numeric(h[[col[1L]]]) else numeric(0))
  }
  if (is.list(h) && length(h) && is.list(h[[1L]])) {
    return(vapply(h, function(z) as.numeric(z$objective %||% NA_real_)[1L],
                  numeric(1)))
  }
  numeric(0)
}

#' Relative changes |v_i - v_{i-1}| / |v_{i-1}| of an objective path, no
#' floor (0 when both values are 0; Inf when not finite). Element i - 1 is
#' the change at iteration i.
#' @keywords internal
#' @noRd
.ggcv_rel_changes <- function(v) {
  if (length(v) < 2L) return(numeric(0))
  a <- v[-length(v)]
  b <- v[-1L]
  rel <- abs(b - a) / abs(a)
  rel[is.finite(a) & is.finite(b) & a == 0 & b == 0] <- 0
  rel[!is.finite(rel)] <- Inf
  rel
}

#' First iteration n_star at which the relative change of the objective has
#' been at most `tol` at two consecutive iterations (n_star - 1 and n_star);
#' `NA` when that never happens on the path.
#' @keywords internal
#' @noRd
.ggcv_first_settled <- function(v, tol) {
  ok <- .ggcv_rel_changes(v) <= tol
  if (length(ok) < 2L) return(NA_integer_)
  hit <- which(ok[-length(ok)] & ok[-1L])
  if (length(hit)) as.integer(hit[1L] + 2L) else NA_integer_
}

#' Adaptive budget from the pilot's n_star: min(cap, ceiling(1.5 n_star) + 2),
#' the cap when the pilot never settled.
#' @keywords internal
#' @noRd
.ggcv_adaptive_budget <- function(n_star, cap) {
  if (is.na(n_star)) return(as.integer(cap))
  as.integer(min(cap, ceiling(1.5 * n_star) + 2))
}

#' Pilot of the adaptive budget at one lambda: cold fits at y with tol = 0 for
#' b0, 2 b0, 4 b0, ... iterations (the last one the cap) until the objective
#' rule of .ggcv_first_settled() holds. A cold run with tol = 0 repeats the
#' iterations of a shorter one, so each run extends the same path. Returns
#' the last fit (or the failed fit's error object), its budget `b`, `n_star`
#' (`NA` when the rule never held), `n_run` (iterations of the last run) and
#' the number of fits run.
#' @keywords internal
#' @noRd
.ggcv_pilot <- function(st, lambda, fit_safe) {
  b <- min(st$pilot_b0, st$cap)
  n_fits <- 0L
  repeat {
    ctrl <- st$control
    ctrl[[st$iter_name]] <- as.integer(b)
    f <- fit_safe(st$y, lambda, ctrl)
    n_fits <- n_fits + 1L
    if (!inherits(f, "ttpspline")) {
      return(list(fit = f, b = b, n_star = NA_integer_, n_run = NA_integer_,
                  n_fits = n_fits))
    }
    v <- .ggcv_objective_path(f)
    ns <- .ggcv_first_settled(v, st$fit_tol)
    # stop when the rule held, at the cap, or when the fit stopped early for
    # another reason (a longer run would end at the same iteration)
    if (!is.na(ns) || b >= st$cap || length(v) < b) {
      return(list(fit = f, b = b, n_star = ns, n_run = length(v), n_fits = n_fits))
    }
    b <- min(2L * b, st$cap)
  }
}

#' TRUE when Monte Carlo probe contributions do not form a usable GDF: a
#' non-finite value, an outlier |c_j| > n_eff (a smoother's r' J r is at
#' most n_eff), a mean <= 0, or (M >= 2) a per-probe coefficient of
#' variation sd / mean above .ggcv_cv_limit(mean).
#' @keywords internal
#' @noRd
.ggcv_unstable <- function(contrib, n_eff) {
  cc <- as.numeric(contrib)
  if (!length(cc) || any(!is.finite(cc)) || any(abs(cc) > n_eff)) return(TRUE)
  gdf <- mean(cc)
  if (!(gdf > 0)) return(TRUE)
  length(cc) >= 2L && stats::sd(cc) / gdf > .ggcv_cv_limit(gdf)
}

#' Largest per-probe coefficient of variation of a usable Monte Carlo GDF:
#' max(1, 3 * sqrt(2 / max(gdf, 1))).
#'
#' Why it depends on the GDF: for a symmetric smoother S with eigenvalues
#' 0 <= e_i <= 1, a Rademacher probe r gives
#' Var(r' S r) = 2 * sum_{i != j} S_ij^2 <= 2 * sum(e_i^2) <= 2 * sum(e_i)
#' = 2 * GDF, so the per-probe CV is at most sqrt(2 / GDF). A fixed limit of
#' 1 sits above that bound only while GDF > 2, and close to it when the GDF
#' is small: with the M_search = 4 probes of the search, converged and
#' heavily smoothed evaluations (GDF 3.4 to 4.7, bound 0.65 to 0.77) crossed
#' it by chance (11 of 410 search rows on smooth truths at rank = k = 6),
#' and re-evaluation with the same probes cannot clear such a flag. Three
#' times the bound raises the limit only where it exceeds 1 (GDF < 18) and
#' never lowers it below 1. The absolute checks of .ggcv_unstable()
#' (non-finite or |c_j| > n_eff contributions, GDF <= 0) do not depend on
#' it. Trade-off: one outlier among M positive contributions gives a sample
#' CV below sqrt(M), so with M = 4 (CV < 2) a single positive outlier is
#' left to the absolute checks where GDF <= 4.5 (limit >= 2); with
#' M = 16 the CV test still sees it above GDF 1.125. The TT map is not
#' exactly symmetric (ALS / PIRLS-ALS; dmu / dy of a Poisson fit is only
#' similar to a symmetric matrix), one more reason for the factor 3.
#' @keywords internal
#' @noRd
.ggcv_cv_limit <- function(gdf) {
  max(1, 3 * sqrt(2 / max(as.numeric(gdf), 1)))
}

#' Control with the iteration budget multiplied (budget check).
#' @keywords internal
#' @noRd
.ggcv_budget_control <- function(ctrl, mult) {
  mult <- as.integer(mult)
  if (mult == 1L) return(ctrl)
  ctrl$pirls_maxit <- mult * as.integer(ctrl$pirls_maxit)
  ctrl$max_sweeps <- mult * as.integer(ctrl$max_sweeps)
  ctrl
}

#' Cache key of the iteration budget: the control budgets (fixed), or the
#' caps, fit_tol and the multiplier (adaptive).
#' @keywords internal
#' @noRd
.ggcv_budget_key <- function(st, mult) {
  if (identical(st$budget, "adaptive")) {
    return(sprintf("a%d/%d/%g/x%d", as.integer(st$control$pirls_maxit),
                   as.integer(st$control$max_sweeps), st$fit_tol,
                   as.integer(mult)))
  }
  ctrl <- .ggcv_budget_control(st$control, mult)
  sprintf("%d/%d", as.integer(ctrl$pirls_maxit), as.integer(ctrl$max_sweeps))
}

#' Controls of the fixed-budget fits of one point: `est` (estimator at y) and
#' `probe` (reference and probes). `mult` multiplies the control budget
#' (budget check, probe validation). Fixed mode: the control budgets times
#' `mult`. Adaptive: the budgeted iteration count set to `B` (`mult` already
#' included); warm probes get `probe_budget` (times `mult`) or a third of `B`.
#' @keywords internal
#' @noRd
.ggcv_fixed_controls <- function(st, B, mult) {
  if (!identical(st$budget, "adaptive")) {
    return(list(est = .ggcv_budget_control(st$control, mult),
                probe = .ggcv_budget_control(st$control_probe, mult)))
  }
  ctrl <- .ggcv_budget_control(st$control, mult)
  ctrl[[st$iter_name]] <- as.integer(B)
  ctrl_p <- ctrl
  if (identical(st$probe_init, "warm")) {
    ctrl_p[[st$iter_name]] <- if (is.null(st$probe_budget)) {
      max(2L, as.integer(ceiling(B / 3)))
    } else {
      as.integer(mult) * st$probe_budget
    }
  }
  list(est = ctrl, probe = ctrl_p)
}

#' lapply / mclapply over fit jobs.
#' @keywords internal
#' @noRd
.ggcv_apply <- function(jobs, FUN, n_cores) {
  if (n_cores > 1L && length(jobs) > 1L) {
    parallel::mclapply(jobs, FUN, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(jobs, FUN)
  }
}

# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------

#' GCV / UBRE, the deviance-unit score and its slope in GDF.
#' @keywords internal
#' @noRd
.ggcv_criterion <- function(dev, gdf, st) {
  n <- st$n_eff
  valid <- is.finite(gdf) && gdf > 0 && gdf < n && is.finite(dev)
  gcv <- if (valid) n * dev / (n - gdf)^2 else Inf
  ubre <- if (valid && !is.null(st$scale)) {
    dev / n + 2 * st$scale * gdf / n - st$scale
  } else if (valid) NA_real_ else Inf
  if (identical(st$criterion, "ubre")) {
    score <- ubre
    # AIC scale D / scale + 2 GDF (Poisson: D + 2 GDF), so that tol and the
    # paired SE do not depend on the units of y
    score_dev <- if (valid) dev / st$scale + 2 * gdf else Inf
    g <- 2
  } else {
    score <- gcv
    score_dev <- if (valid) n * log(gcv) else Inf
    g <- if (valid) 2 * n / (n - gdf) else NA_real_
  }
  list(valid = valid, score = score, score_dev = score_dev, gcv = gcv,
       ubre = ubre, g = g)
}

#' Monte Carlo GDF from the first M probe contributions, then the criterion.
#' `gdf_cv` = sd(contributions) / GDF, the per-probe coefficient of
#' variation; `stable` = the contributions form a usable GDF (all finite, no
#' |c_j| > n_eff, GDF > 0, `gdf_cv` at most .ggcv_cv_limit(GDF); see
#' .ggcv_unstable()).
#' @keywords internal
#' @noRd
.ggcv_score <- function(dev, contrib, M, st) {
  cc <- contrib[seq_len(min(as.integer(M), length(contrib)))]
  ok <- is.finite(cc)
  gdf <- if (any(ok)) mean(cc[ok]) else NA_real_
  sdc <- if (sum(ok) >= 2L) stats::sd(cc[ok]) else NA_real_
  gdf_se <- if (sum(ok) >= 2L) sdc / sqrt(sum(ok)) else NA_real_
  gdf_cv <- if (is.finite(sdc) && is.finite(gdf) && gdf > 0) sdc / gdf else NA_real_
  c(list(gdf = gdf, gdf_se = gdf_se, gdf_cv = gdf_cv,
         stable = !.ggcv_unstable(cc, st$n_eff), M_ok = sum(ok), contrib = cc),
    .ggcv_criterion(dev, gdf, st))
}

#' Paired comparison of two evaluations on their common probes.
#'
#' `diff = a$score_dev - b$score_dev` (positive: `a` is worse) and the SE of
#' the paired per-probe differences `g_a c_a[j] - g_b c_b[j]`; `improves` is
#' `TRUE` when `b` improves on `a`: `diff > max(tol, 2 * se)`.
#' @keywords internal
#' @noRd
.ggcv_paired <- function(a, b, tol = 1) {
  diff <- a$score_dev - b$score_dev
  M <- min(length(a$contrib), length(b$contrib))
  ca <- a$contrib[seq_len(M)]
  cb <- b$contrib[seq_len(M)]
  ok <- is.finite(ca) & is.finite(cb)
  J <- sum(ok)
  se <- if (J >= 2L && is.finite(a$g) && is.finite(b$g)) {
    stats::sd(a$g * ca[ok] - b$g * cb[ok]) / sqrt(J)
  } else {
    NA_real_
  }
  thr <- if (is.finite(se)) max(tol, 2 * se) else tol
  list(diff = diff, se = se, J = J, improves = isTRUE(diff > thr))
}

#' Final choice among candidates scored on common probes (pure function).
#'
#' Each candidate carries `theta`, `iso` (logical), `score_dev`, `g`,
#' `contrib`, `converged` and `stable` (missing = `TRUE`). Eligible: finite
#' `score_dev`, not in `exclude`, converged and stable; when none is, the
#' stable ones; when none of those either, all of them. `best` is the
#' eligible argmin of `score_dev`; an eligible candidate qualifies when
#' `best` does not improve on it under the paired rule,
#' `diff <= max(tol, 2 * SE)` (`best` qualifies trivially). Regularity:
#' isotropic before grouped, among isotropic candidates larger theta first.
#' The winner is the qualifying isotropic candidate with the largest theta,
#' else `best`. Returns `winner`, `best` (indices, `NA` when nothing is
#' eligible), `eligible`, `qualifies`, the paired `diff` / `se` against
#' `best`, `none_converged` (no converged candidate available) and
#' `none_stable` (no stable one).
#' @keywords internal
#' @noRd
.ggcv_choose <- function(cands, tol = 1, exclude = integer(0)) {
  n <- length(cands)
  sdf <- vapply(cands, function(r) as.numeric(r$score_dev)[1L], numeric(1))
  conv <- vapply(cands, function(r) isTRUE(r$converged), logical(1))
  stab <- vapply(cands, function(r) !isFALSE(r$stable), logical(1))
  iso <- vapply(cands, function(r) isTRUE(r$iso), logical(1))
  t1 <- vapply(cands, function(r) as.numeric(r$theta)[1L], numeric(1))
  avail <- is.finite(sdf) & !(seq_len(n) %in% exclude)
  none_conv <- !any(avail & conv)
  none_stab <- !any(avail & stab)
  elig <- avail & conv & stab
  if (!any(elig)) elig <- avail & stab
  if (!any(elig)) elig <- avail
  out <- list(winner = NA_integer_, best = NA_integer_, eligible = elig,
              qualifies = rep(FALSE, n), diff = rep(NA_real_, n),
              se = rep(NA_real_, n), none_converged = none_conv,
              none_stable = none_stab)
  if (!any(elig)) return(out)
  ie <- which(elig)
  b <- ie[which.min(sdf[ie])]
  pr <- lapply(cands, function(r) .ggcv_paired(r, cands[[b]], tol))
  out$diff <- vapply(pr, function(p) p$diff, numeric(1))
  out$se <- vapply(pr, function(p) p$se, numeric(1))
  out$qualifies <- elig & !vapply(pr, function(p) p$improves, logical(1))
  out$qualifies[b] <- TRUE
  iq <- which(out$qualifies & iso)
  out$best <- b
  out$winner <- if (length(iq)) iq[which.max(t1[iq])] else b
  out
}

# ---------------------------------------------------------------------------
# evaluation
# ---------------------------------------------------------------------------

#' Score several points (rows of `thetas`, group scale) in one batch.
#'
#' Fits only what the cache lacks at fidelity `M`: per point the estimator
#' (y), the reference (y + eps / 2) and the missing probes, all at the
#' point's fixed budget B (tol = 0, cold unless `probe_init = "warm"`).
#' Adaptive budget: the pilots of the new points run first, one parallel
#' batch (each pilot doubles its own budget, see .ggcv_pilot()), and fix B;
#' the estimators (unless a pilot run is it), references and probes of all
#' points then form a second batch. Fixed budget: all jobs go to one
#' `lapply` / `mclapply` call. Warm probes: estimators first. Probe
#' validation: a point whose first M contributions do not form a usable GDF
#' (.ggcv_unstable()) is re-evaluated whole (estimator, reference, probes
#' 1..M, same probes) at 2B and then 4B, one batch per round. The budget
#' check (`budget_mult = 2`) runs every fit at twice the point's budget and
#' is not validated. The cache is keyed by log10(lambda) and the budget and
#' keeps deviance, per-probe contributions, B, the convergence flag, the
#' number of re-evaluations and the first fit error, never fits.
#' `keep_fit = TRUE` returns the estimator fits (final stage). Each computed
#' point appends one row to the search log.
#' @keywords internal
#' @noRd
.ggcv_eval_batch <- function(thetas, M, st, stage, budget_mult = 1L,
                             keep_fit = FALSE) {
  if (identical(st$gdf_method, "exact")) {
    return(.ggcv_eval_batch_exact(thetas, M, st, stage, budget_mult, keep_fit))
  }
  t0 <- proc.time()[["elapsed"]]
  G <- st$n_groups
  thetas <- matrix(as.numeric(thetas), ncol = G)
  n_pts <- nrow(thetas)
  stage <- rep(as.character(stage), length.out = n_pts)
  M <- as.integer(M)
  mult <- as.integer(budget_mult)
  adaptive <- identical(st$budget, "adaptive")
  if (!is.null(st$probe_bank) && M > ncol(st$probe_bank)) {
    stop(sprintf("The probe bank has %d columns; %d probes requested.",
                 ncol(st$probe_bank), M), call. = FALSE)
  }
  key_of <- function(th, m) {
    paste(c(sprintf("%.6f", th[st$groups]), .ggcv_budget_key(st, m)),
          collapse = "|")
  }
  keys <- vapply(seq_len(n_pts), function(i) key_of(thetas[i, ], mult),
                 character(1))

  # One work item per distinct key that misses probes at this M (or, with
  # keep_fit, needs its estimator). B = fixed budget of its fits: from the
  # cache entry, twice the evaluated point's B (budget check), the control
  # budget (fixed), or this batch's pilot (NA). `mt` multiplies the control
  # budget in fixed mode and the warm probe budget: budget check times
  # probe validation.
  work <- list()
  for (i in which(!duplicated(keys))) {
    ent <- st$cache[[keys[i]]]
    J <- if (is.null(ent)) 0L else length(ent$contrib)
    if (J >= M && !isTRUE(keep_fit)) next
    base <- if (is.null(ent) && mult > 1L) st$cache[[key_of(thetas[i, ], 1L)]]
    B <- if (!is.null(ent)) {
      ent$budget
    } else if (!is.null(base)) {
      mult * base$budget
    } else if (!adaptive) {
      mult * st$cap
    } else {
      NA_integer_
    }
    th_d <- if (is.null(ent)) thetas[i, st$groups] else ent$theta_d
    new <- if (J < M) seq.int(J + 1L, M) else integer(0)
    work[[length(work) + 1L]] <- list(
      row = i, key = keys[i], theta_d = th_d, lambda = 10^th_d, new = new,
      old = if (J > 0L) ent$contrib else numeric(0),
      cnew = rep(NA_real_, length(new)), B = as.integer(B),
      mt = as.integer(if (is.null(ent)) mult * (base$mt %||% 1L) else ent$mt %||% 1L),
      vmult = as.integer(if (!is.null(ent)) ent$vmult %||% 1L else 1L),
      cached = !is.null(ent), conv_pilot = NA, n_iter = ent$n_iter %||% NA_integer_,
      est = NULL, failed = FALSE, n_fits = 0L, reevals = 0L
    )
  }

  fits <- list()
  if (length(work)) {
    # a failed fit returns its error message (kept for the caller's error
    # messages) instead of stopping the batch
    fit_safe <- function(yy, lambda, ctrl, init = NULL) {
      tryCatch(st$map$fit(yy, lambda, ctrl, init = init), error = function(e) {
        structure(list(message = conditionMessage(e)), class = "ggcv_fit_error")
      })
    }
    errs <- rep(NA_character_, length(work))
    note_error <- function(w, x) {
      msg <- if (is.null(x)) "a worker returned no result" else .ggcv_fit_error(x)
      if (is.na(msg)) return(invisible(NULL))
      if (is.na(errs[w])) errs[w] <<- msg
      if (is.na(st$log$first_error)) st$log$first_error <- msg
    }
    # probes of this batch, drawn here in the parent before any fit is
    # dispatched (probe j is the same at every lambda): a failed draw stops
    # the batch with its error instead of turning into NA contributions
    # inside forked workers
    probes <- new.env(parent = emptyenv())
    draw_probes <- function(js) {
      for (j in setdiff(js, as.integer(ls(probes)))) {
        assign(as.character(j), .ggcv_probe(j, st), envir = probes)
      }
    }
    probe_of <- function(j) get(as.character(j), envir = probes)

    # 1. pilots of the new points (adaptive budget): they fix B, and the last
    # pilot run is the estimator when it ran exactly B
    pw <- which(vapply(work, function(wk) is.na(wk$B), logical(1)))
    if (length(pw)) {
      res_p <- .ggcv_apply(pw, function(w) {
        .ggcv_pilot(st, work[[w]]$lambda, fit_safe)
      }, st$n_cores)
      for (q in seq_along(pw)) {
        w <- pw[q]
        p <- res_p[[q]]
        if (!is.list(p) || inherits(p, "try-error") || is.null(p$n_fits)) {
          p <- list(fit = p, b = NA_integer_, n_star = NA_integer_,
                    n_run = NA_integer_, n_fits = 1L)
        }
        B0 <- .ggcv_adaptive_budget(p$n_star, st$cap)
        work[[w]]$B <- as.integer(mult * B0)
        work[[w]]$mt <- mult
        work[[w]]$n_fits <- p$n_fits
        work[[w]]$conv_pilot <- !is.na(p$n_star)
        work[[w]]$n_iter <- if (!is.na(p$n_star)) p$n_star else p$n_run
        if (!inherits(p$fit, "ttpspline")) {
          note_error(w, p$fit)
          work[[w]]$failed <- TRUE  # the fixed-budget fits would fail alike
        } else if (mult == 1L && identical(as.integer(p$b), B0)) {
          work[[w]]$est <- p$fit
        }
      }
    }

    # 2. fixed-budget fits of the work items `ws`: the estimator (unless it
    # is known), the reference and the probes new_of(w), at the item's B
    run_fixed <- function(ws, new_of) {
      draw_probes(sort(unique(unlist(lapply(ws, new_of)))))
      ctl <- vector("list", length(work))
      for (w in ws) ctl[[w]] <- .ggcv_fixed_controls(st, work[[w]]$B, work[[w]]$mt)
      # job = c(w, j): j = -1 estimator, 0 reference, >= 1 probe j
      est_jobs <- list()
      probe_jobs <- list()
      for (w in ws) {
        if (is.null(work[[w]]$est)) est_jobs[[length(est_jobs) + 1L]] <- c(w = w, j = -1L)
        jn <- new_of(w)
        if (length(jn)) {
          probe_jobs <- c(probe_jobs, list(c(w = w, j = 0L)),
                          lapply(jn, function(j) c(w = w, j = j)))
        }
      }
      inits <- vector("list", length(work))
      run <- function(job) {
        w <- job[["w"]]
        j <- job[["j"]]
        wk <- work[[w]]
        if (j < 0L) return(fit_safe(st$y, wk$lambda, ctl[[w]]$est))
        yy <- if (j == 0L) {
          st$y + 0.5 * st$eps * as.numeric(st$support)
        } else {
          st$y + st$eps * (probe_of(j) + 1) / 2 * as.numeric(st$support)
        }
        f <- fit_safe(yy, wk$lambda, ctl[[w]]$probe, init = inits[[w]])
        if (inherits(f, "ttpspline")) as.numeric(f$fitted.values) else f
      }
      set_est <- function(jobs, res) {
        for (q in seq_along(jobs)) {
          w <- jobs[[q]][["w"]]
          work[[w]]$n_fits <<- work[[w]]$n_fits + 1L
          if (inherits(res[[q]], "ttpspline")) {
            work[[w]]$est <<- res[[q]]
          } else {
            note_error(w, res[[q]])
          }
        }
      }
      if (identical(st$probe_init, "warm")) {
        # probes restart from the estimator's cores: estimators first
        if (length(est_jobs)) set_est(est_jobs, .ggcv_apply(est_jobs, run, st$n_cores))
        for (w in ws) {
          if (inherits(work[[w]]$est, "ttpspline")) {
            inits[[w]] <- .tt_clone_cores(work[[w]]$est$cores)
          }
        }
        probe_jobs <- Filter(function(jb) !is.null(inits[[jb[["w"]]]]), probe_jobs)
        res <- if (length(probe_jobs)) {
          .ggcv_apply(probe_jobs, run, st$n_cores)
        } else {
          list()
        }
      } else {
        all_jobs <- c(est_jobs, probe_jobs)
        res_all <- if (length(all_jobs)) {
          .ggcv_apply(all_jobs, run, st$n_cores)
        } else {
          list()
        }
        ne <- length(est_jobs)
        if (ne) set_est(est_jobs, res_all[seq_len(ne)])
        res <- res_all[seq_along(probe_jobs) + ne]
      }
      # collect: references and probe fits -> contributions
      out <- lapply(seq_along(work), function(w) list(ref = NULL, mus = list()))
      for (q in seq_along(probe_jobs)) {
        w <- probe_jobs[[q]][["w"]]
        j <- probe_jobs[[q]][["j"]]
        val <- res[[q]]
        work[[w]]$n_fits <<- work[[w]]$n_fits + 1L
        if (is.numeric(val) && length(val) == st$n_rows) {
          if (j == 0L) out[[w]]$ref <- val else out[[w]]$mus[[as.character(j)]] <- val
        } else {
          note_error(w, val)
        }
      }
      for (w in ws) {
        jn <- new_of(w)
        ref <- out[[w]]$ref
        work[[w]]$cnew <<- vapply(jn, function(j) {
          m <- out[[w]]$mus[[as.character(j)]]
          if (is.null(ref) || is.null(m)) return(NA_real_)
          sum(probe_of(j) * (m - ref)) * 2 / st$eps
        }, numeric(1))
      }
      invisible(NULL)
    }
    contrib_of <- function(w) c(work[[w]]$old, work[[w]]$cnew)

    ws0 <- which(!vapply(work, function(wk) wk$failed, logical(1)))
    if (length(ws0)) run_fixed(ws0, function(w) work[[w]]$new)

    # 3. probe validation (not for the budget check): a point whose first M
    # contributions do not form a usable GDF is re-evaluated at 2B, then 4B
    # (points whose estimator failed have no score to validate)
    if (mult == 1L) {
      repeat {
        redo <- Filter(function(w) {
          work[[w]]$vmult < 4L && inherits(work[[w]]$est, "ttpspline") &&
            .ggcv_unstable(contrib_of(w)[seq_len(M)], st$n_eff)
        }, ws0)
        if (!length(redo)) break
        for (w in redo) {
          work[[w]]$B <- 2L * work[[w]]$B
          work[[w]]$mt <- 2L * work[[w]]$mt
          work[[w]]$vmult <- 2L * work[[w]]$vmult
          work[[w]]$reevals <- work[[w]]$reevals + 1L
          work[[w]]["est"] <- list(NULL)
          work[[w]]$old <- numeric(0)
          work[[w]]$cnew <- rep(NA_real_, M)
        }
        run_fixed(redo, function(w) seq_len(M))
      }
    }

    # 4. cache entries: estimator -> deviance and convergence
    for (w in seq_along(work)) {
      wk <- work[[w]]
      f0 <- wk$est
      ent <- st$cache[[wk$key]]
      ni <- wk$n_iter
      if (inherits(f0, "ttpspline")) {
        mu0 <- as.numeric(f0$fitted.values)
        dev <- glm_deviance(st$fam, st$y, mu0, weights = st$map$weights)
        rel <- .ggcv_last_rel_change(f0)
        if (adaptive && mult == 1L) {
          # adaptive: converged iff the pilot settled (kept from the cache
          # when the point was evaluated before)
          conv <- if (wk$cached) isTRUE(ent$converged) else {
            isTRUE(wk$conv_pilot) && !isFALSE(f0$converged)
          }
        } else {
          ni <- .ggcv_n_iter(f0, st)
          conv <- is.na(rel) || rel <= st$conv_tol
        }
      } else {
        dev <- NA_real_
        rel <- NA_real_
        conv <- FALSE
      }
      st$cache[[wk$key]] <- list(
        theta_d = wk$theta_d, deviance = dev, contrib = contrib_of(w),
        converged = conv, last_rel_change = rel, budget = wk$B, mt = wk$mt,
        vmult = wk$vmult, n_iter = ni,
        reevals = (if (wk$cached) ent$reevals %||% 0L else 0L) + wk$reevals,
        # first failed fit of this point (estimator first), NA if none
        error = if (!is.na(ent$error %||% NA_character_)) ent$error else errs[w]
      )
      if (isTRUE(keep_fit)) fits[[wk$key]] <- f0
    }
  }

  recs <- lapply(seq_len(n_pts), function(i) {
    ent <- st$cache[[keys[i]]]
    sc <- .ggcv_score(ent$deviance, ent$contrib, M, st)
    rec <- c(list(theta = thetas[i, ], theta_d = ent$theta_d,
                  lambda = 10^ent$theta_d, key = keys[i], M = M,
                  deviance = ent$deviance, converged = ent$converged,
                  last_rel_change = ent$last_rel_change, budget = ent$budget,
                  n_iter = ent$n_iter, reevals = ent$reevals %||% 0L,
                  error = ent$error), sc)
    rec$n_fits <- 0L
    if (isTRUE(keep_fit)) rec$fit <- fits[[keys[i]]]
    rec
  })

  # search log: one row per computed point, wall time shared equally
  elapsed <- proc.time()[["elapsed"]] - t0
  for (wk in work) {
    i <- wk$row
    recs[[i]]$n_fits <- wk$n_fits
    rec <- recs[[i]]
    st$log$rows[[length(st$log$rows) + 1L]] <- data.frame(
      stage = stage[i], t(stats::setNames(rec$theta, paste0("theta", seq_len(G)))),
      score = rec$score, score_dev = rec$score_dev, gcv = rec$gcv,
      ubre = rec$ubre, gdf = rec$gdf, gdf_se = rec$gdf_se,
      gdf_cv = rec$gdf_cv, stable = rec$stable, deviance = rec$deviance,
      M = M, M_ok = rec$M_ok, converged = rec$converged,
      last_rel_change = rec$last_rel_change, budget = rec$budget,
      n_iter = rec$n_iter, reevals = wk$reevals, n_fits = wk$n_fits,
      time_s = elapsed / length(work), stringsAsFactors = FALSE
    )
    if (isTRUE(st$verbose)) {
      message(sprintf(
        "%s %-12s theta=(%s) %s=%.6g [%.2f] GDF=%.2f (se %.2f) dev=%.6g M=%d B=%d fits=%d%s%s [%.1fs]",
        st$prefix, stage[i], paste(sprintf("%.3f", rec$theta), collapse = ","),
        st$criterion, rec$score, rec$score_dev, rec$gdf, rec$gdf_se,
        rec$deviance, M, rec$budget, wk$n_fits,
        if (isTRUE(rec$converged)) "" else " unconverged",
        if (isTRUE(rec$stable)) "" else " unstable", elapsed / length(work)
      ))
    }
  }
  recs
}

#' Score several points with the exact GDF (`gdf_method = "exact"`).
#'
#' One job per uncached point, all points in one parallel batch: the
#' estimator stopped at `fit_tol` (the control budget is the cap), then
#' .tt_exact_edf_core(), which polishes the cores to a stationary point by
#' Newton and takes the trace there. Deviance and fit are the polished ones.
#' The point counts as converged when Newton reached a stationary point with
#' no negative curvature (a local minimum). There are no probes: `contrib` is
#' empty, `gdf_se` 0, so paired comparisons reduce to the tolerance. The
#' budget check (`budget_mult = 2`) refits at `fit_tol / 100` with twice the
#' cap before polishing; a different GDF there means another local minimum. Cache keys ignore `M`, so the final stage reuses the search fits'
#' scores (fits are kept only when `keep_fit = TRUE`, which refits the
#' estimator when the point is cached).
#' @keywords internal
#' @noRd
.ggcv_eval_batch_exact <- function(thetas, M, st, stage, budget_mult = 1L,
                                   keep_fit = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  G <- st$n_groups
  thetas <- matrix(as.numeric(thetas), ncol = G)
  n_pts <- nrow(thetas)
  stage <- rep(as.character(stage), length.out = n_pts)
  mult <- as.integer(budget_mult)
  keys <- vapply(seq_len(n_pts), function(i) {
    paste(c(sprintf("%.6f", thetas[i, st$groups]), "exact", mult), collapse = "|")
  }, character(1))
  ctrl <- st$control_exact
  if (mult > 1L) {
    ctrl$tol <- ctrl$tol / 100
    ctrl$pirls_maxit <- mult * as.integer(ctrl$pirls_maxit)
    ctrl$max_sweeps <- mult * as.integer(ctrl$max_sweeps)
  }
  todo <- which(!duplicated(keys) &
                  vapply(keys, function(k) is.null(st$cache[[k]]) || isTRUE(keep_fit),
                         logical(1)))
  # cached points that only need their fit (keep_fit) skip the Hessian
  need_edf <- vapply(todo, function(i) is.null(st$cache[[keys[i]]]), logical(1))
  one <- function(q) {
    i <- todo[q]
    th_d <- thetas[i, st$groups]
    f <- tryCatch(st$map$fit(st$y, 10^th_d, ctrl),
                  error = function(e) structure(list(message = conditionMessage(e)),
                                                class = "ggcv_fit_error"))
    if (!inherits(f, "ttpspline")) return(list(fit = f, error = .ggcv_fit_error(f)))
    yy <- as.numeric(f$y %||% st$y)
    e <- tryCatch({
      basis <- eval_marginal_bases(as.matrix(f$X), f$knots, f$degree, cyclic = f$cyclic)
      .tt_exact_edf_core(f$cores, f$intercept, basis, yy, f$family,
                         offset = f$offset, weights = f$weights, lambda = f$lambda,
                         penalty_order = f$penalty_order %||% 2L)
    }, error = function(err) list(edf = NA_real_, error = conditionMessage(err)))
    if (is.finite(e$edf %||% NA_real_)) f <- .tt_fit_polished(f, e, yy)
    if (!need_edf[q]) return(list(fit = f, edf = NULL, error = NA_character_))
    # the trace is the derivative of the polished (stationary) map; its
    # deviance and convergence are the ones scored
    list(fit = f, edf = e, error = e$error %||% NA_character_)
  }
  res <- .ggcv_apply(seq_along(todo), one, st$n_cores)
  fits <- list()
  for (q in seq_along(todo)) {
    i <- todo[q]
    r <- res[[q]]
    f <- if (is.list(r)) r$fit else NULL
    ent_old <- st$cache[[keys[i]]]
    if (!need_edf[q]) {
      if (isTRUE(keep_fit)) fits[[keys[i]]] <- f
      next
    }
    if (inherits(f, "ttpspline")) {
      dev <- glm_deviance(st$fam, st$y, as.numeric(f$fitted.values),
                          weights = st$map$weights)
      ni <- .ggcv_n_iter(f, st)
      # converged = Newton reached a stationary point that is a local minimum
      # (no negative curvature); the ALS cap no longer decides it
      conv <- isTRUE(r$edf$converged) && identical(as.integer(r$edf$n_negative), 0L)
      ent <- list(theta_d = thetas[i, st$groups], deviance = dev,
                  gdf = as.numeric(r$edf$edf), contrib = numeric(0),
                  converged = conv, last_rel_change = .ggcv_last_rel_change(f),
                  budget = ni, n_iter = ni, reevals = 0L,
                  grad_norm = r$edf$grad_norm %||% NA_real_,
                  n_negative = r$edf$n_negative %||% NA_integer_,
                  edf_time = r$edf$time_s %||% NA_real_,
                  error = if (is.na(r$error)) NA_character_ else r$error)
      if (!is.na(r$error) && is.na(st$log$first_error)) st$log$first_error <- r$error
    } else {
      msg <- if (is.list(r)) r$error else "a worker returned no result"
      if (is.na(st$log$first_error)) st$log$first_error <- msg
      ent <- list(theta_d = thetas[i, st$groups], deviance = NA_real_, gdf = NA_real_,
                  contrib = numeric(0), converged = FALSE, last_rel_change = NA_real_,
                  budget = NA_integer_, n_iter = NA_integer_, reevals = 0L,
                  grad_norm = NA_real_, n_negative = NA_integer_, edf_time = NA_real_,
                  error = msg)
    }
    if (is.null(ent_old) || !isTRUE(keep_fit)) st$cache[[keys[i]]] <- ent
    if (isTRUE(keep_fit)) fits[[keys[i]]] <- f
  }
  recs <- lapply(seq_len(n_pts), function(i) {
    ent <- st$cache[[keys[i]]]
    cr <- .ggcv_criterion(ent$deviance, ent$gdf, st)
    rec <- c(list(theta = thetas[i, ], theta_d = ent$theta_d,
                  lambda = 10^ent$theta_d, key = keys[i], M = 0L,
                  deviance = ent$deviance, converged = ent$converged,
                  last_rel_change = ent$last_rel_change, budget = ent$budget,
                  n_iter = ent$n_iter, reevals = 0L, error = ent$error,
                  gdf = ent$gdf, gdf_se = 0, gdf_cv = 0,
                  stable = is.finite(ent$gdf), M_ok = 0L, contrib = numeric(0),
                  grad_norm = ent$grad_norm), cr)
    rec$n_fits <- 0L
    if (isTRUE(keep_fit)) rec$fit <- fits[[keys[i]]]
    rec
  })
  elapsed <- proc.time()[["elapsed"]] - t0
  for (q in seq_along(todo)) {
    i <- todo[q]
    recs[[i]]$n_fits <- 1L
    rec <- recs[[i]]
    st$log$rows[[length(st$log$rows) + 1L]] <- data.frame(
      stage = stage[i], t(stats::setNames(rec$theta, paste0("theta", seq_len(G)))),
      score = rec$score, score_dev = rec$score_dev, gcv = rec$gcv,
      ubre = rec$ubre, gdf = rec$gdf, gdf_se = 0, gdf_cv = 0,
      stable = rec$stable, deviance = rec$deviance, M = 0L, M_ok = 0L,
      converged = rec$converged, last_rel_change = rec$last_rel_change,
      budget = rec$budget, n_iter = rec$n_iter, reevals = 0L, n_fits = 1L,
      time_s = elapsed / length(todo), stringsAsFactors = FALSE
    )
    if (isTRUE(st$verbose)) {
      message(sprintf(
        "%s %-12s theta=(%s) %s=%.6g [%.2f] GDF=%.2f (exact) dev=%.6g iter=%s%s [%.1fs]",
        st$prefix, stage[i], paste(sprintf("%.3f", rec$theta), collapse = ","),
        st$criterion, rec$score, rec$score_dev, rec$gdf, rec$deviance,
        rec$n_iter, if (isTRUE(rec$converged)) "" else " unconverged",
        elapsed / length(todo)
      ))
    }
  }
  recs
}

#' Exact finite-difference trace: one unit perturbation per support row
#' (small problems only; used by `probes = "unit"` in tt_gdf_array() and
#' tt_gdf()). Adaptive budget: the pilot fixes B and the estimator (the
#' cold fit at B, the base of the differences) and the unit fits run at B.
#' @keywords internal
#' @noRd
.ggcv_eval_unit <- function(theta, st) {
  t0 <- proc.time()[["elapsed"]]
  lambda <- 10^as.numeric(theta)[st$groups]
  eps <- st$eps
  cells <- which(st$support)
  adaptive <- identical(st$budget, "adaptive")
  # failed fits return their error message (see .ggcv_eval_batch())
  fit_safe <- function(yy, ctrl, init = NULL) {
    tryCatch(st$map$fit(yy, lambda, ctrl, init = init), error = function(e) {
      structure(list(message = conditionMessage(e)), class = "ggcv_fit_error")
    })
  }
  fit_cell <- function(i, ctrl, init) {
    yy <- st$y
    yy[i] <- yy[i] + eps
    f <- fit_safe(yy, ctrl, init = init)
    if (inherits(f, "ttpspline")) as.numeric(f$fitted.values)[i] else f
  }
  fit_base <- function(ctrl, init) {
    f <- fit_safe(st$y, ctrl, init = init)
    if (inherits(f, "ttpspline")) as.numeric(f$fitted.values) else f
  }
  res <- list()
  fit0 <- NULL
  # first failed fit, the estimator first (fit0 and res are set below)
  first_error <- function() {
    msg <- .ggcv_fit_error(fit0)
    for (x in res) {
      if (!is.na(msg)) break
      msg <- .ggcv_fit_error(x)
    }
    if (!is.na(msg) && is.na(st$log$first_error)) st$log$first_error <- msg
    msg
  }
  as_num <- function(v) if (is.numeric(v) && length(v) == 1L) v else NA_real_
  n_iter <- NA_integer_
  n_fits <- 0L
  conv_pilot <- FALSE
  B <- st$cap
  if (adaptive) {
    p <- .ggcv_pilot(st, lambda, function(yy, lam, ctrl, init = NULL) {
      fit_safe(yy, ctrl, init = init)
    })
    n_fits <- p$n_fits
    B <- .ggcv_adaptive_budget(p$n_star, st$cap)
    n_iter <- if (!is.na(p$n_star)) p$n_star else p$n_run
    conv_pilot <- !is.na(p$n_star)
    # the last pilot run is the estimator when it ran exactly B; a failed
    # pilot is kept as the estimator's error
    if (!inherits(p$fit, "ttpspline") || identical(as.integer(p$b), B)) fit0 <- p$fit
  }
  ctl <- .ggcv_fixed_controls(st, B, 1L)
  if (adaptive && is.null(fit0)) {
    fit0 <- fit_safe(st$y, ctl$est)
    n_fits <- n_fits + 1L
  }
  base <- NULL
  if (identical(st$probe_init, "warm")) {
    if (!adaptive) {
      fit0 <- fit_safe(st$y, st$control)
      n_fits <- n_fits + 1L
    }
    if (inherits(fit0, "ttpspline")) {
      init <- .tt_clone_cores(fit0$cores)
      res <- .ggcv_apply(c(0L, cells), function(i) {
        if (i > 0L) fit_cell(i, ctl$probe, init) else fit_base(ctl$probe, init)
      }, st$n_cores)
      base <- res[[1L]]
      if (!is.numeric(base) || length(base) != st$n_rows) base <- NULL
      n_fits <- n_fits + 1L + length(cells)
    }
  } else if (adaptive) {
    # cold: the estimator (cold fit at B) is the base of the unit fits at B
    if (inherits(fit0, "ttpspline")) {
      res <- c(list(NULL), .ggcv_apply(cells, function(i) {
        fit_cell(i, ctl$probe, NULL)
      }, st$n_cores))
      base <- as.numeric(fit0$fitted.values)
      n_fits <- n_fits + length(cells)
    }
  } else {
    # cold, fixed: the unit fits do not depend on the estimator -> one batch
    res <- .ggcv_apply(c(0L, cells), function(i) {
      if (i == 0L) fit_safe(st$y, st$control) else fit_cell(i, st$control_probe, NULL)
    }, st$n_cores)
    fit0 <- res[[1L]]
    base <- if (inherits(fit0, "ttpspline")) as.numeric(fit0$fitted.values) else NULL
    n_fits <- 1L + length(cells)
  }
  err <- first_error()
  if (!inherits(fit0, "ttpspline")) {
    return(c(list(theta = theta, deviance = NA_real_, contrib = numeric(0),
                  gdf = NA_real_, gdf_se = NA_real_, gdf_cv = NA_real_,
                  stable = FALSE, M_ok = 0L, fit = NULL, converged = FALSE,
                  last_rel_change = NA_real_, budget = B, n_iter = n_iter,
                  reevals = 0L, n_fits = n_fits, error = err,
                  time_s = proc.time()[["elapsed"]] - t0),
             .ggcv_criterion(NA_real_, NA_real_, st)))
  }
  if (is.null(base)) base <- rep(NA_real_, st$n_rows)
  mi <- if (length(res) > 1L) {
    vapply(res[-1L], as_num, numeric(1))
  } else {
    rep(NA_real_, length(cells))
  }
  contrib <- (mi - base[cells]) / eps
  ok <- is.finite(contrib)
  gdf <- if (all(ok)) sum(contrib) else NA_real_
  dev <- glm_deviance(st$fam, st$y, as.numeric(fit0$fitted.values),
                      weights = st$map$weights)
  rel <- .ggcv_last_rel_change(fit0)
  conv <- if (adaptive) {
    conv_pilot && !isFALSE(fit0$converged)
  } else {
    is.na(rel) || rel <= st$conv_tol
  }
  if (!adaptive) n_iter <- .ggcv_n_iter(fit0, st)
  c(list(theta = theta, deviance = dev, contrib = contrib, gdf = gdf,
         gdf_se = 0, gdf_cv = 0, stable = all(ok), M_ok = sum(ok), fit = fit0,
         converged = conv, last_rel_change = rel, budget = B, n_iter = n_iter,
         reevals = 0L, n_fits = n_fits, error = err,
         time_s = proc.time()[["elapsed"]] - t0),
    .ggcv_criterion(dev, gdf, st))
}

#' Search log as a data frame.
#' @keywords internal
#' @noRd
.ggcv_log_table <- function(st) {
  rows <- st$log$rows
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' GDF and scores of one map at one smoothing vector (tt_gdf(),
#' tt_gdf_array()): seeded Rademacher probes, unit probes (exact trace) or
#' a user probe bank (matrix; M = its number of columns).
#' @keywords internal
#' @noRd
.ggcv_gdf_point <- function(map, lambda, M, probes, scale, probe_init,
                            probe_budget, epsilon_rel, n_cores, seed, budget,
                            fit_tol, control, gdf_method = "mc") {
  t0 <- proc.time()[["elapsed"]]
  bank <- NULL
  if (is.character(probes)) {
    probes <- match.arg(probes, c("rademacher", "unit"))
  } else {
    bank <- probes
    probes <- "matrix"
  }
  d <- map$d
  lam <- as.numeric(lambda)
  if (!length(lam) || anyNA(lam) || any(!is.finite(lam)) || any(lam < 0) ||
      !length(lam) %in% c(1L, d)) {
    stop(sprintf("`lambda` must be non-negative, scalar or of length d = %d.", d),
         call. = FALSE)
  }
  st <- .ggcv_setup(
    map, groups = seq_len(d), criterion = "gcv", scale = scale,
    probe_init = probe_init, probe_budget = probe_budget,
    epsilon_rel = epsilon_rel, n_cores = n_cores, seed = seed,
    control = control, budget = budget, fit_tol = fit_tol, probes = bank,
    gdf_method = if (identical(gdf_method, "exact")) "exact" else "mc"
  )
  theta <- log10(rep(lam, length.out = d))
  ev <- if (identical(st$gdf_method, "exact")) {
    .ggcv_eval_batch(matrix(theta, nrow = 1L), 0L, st, stage = "gdf",
                     keep_fit = TRUE)[[1L]]
  } else if (identical(probes, "unit")) {
    .ggcv_eval_unit(theta, st)
  } else {
    MM <- if (is.null(st$probe_bank)) as.integer(M) else ncol(st$probe_bank)
    if (length(MM) != 1L || is.na(MM) || MM < 1L) {
      stop("`M` must be a positive integer.", call. = FALSE)
    }
    .ggcv_eval_batch(matrix(theta, nrow = 1L), MM, st, stage = "gdf",
                     keep_fit = TRUE)[[1L]]
  }
  # a failed estimator (e.g. invalid model arguments) is an error, not NA;
  # failed perturbation fits leave NA contributions and a warning
  lab <- if (identical(map$mode, "array")) "tt_gdf_array" else "tt_gdf"
  msg <- ev$error %||% NA_character_
  if (!inherits(ev$fit, "ttpspline")) {
    stop(sprintf("%s: the fit at lambda = (%s) failed: %s", lab,
                 paste(signif(rep(lam, length.out = d), 4), collapse = ", "),
                 if (is.na(msg)) "no fit returned" else msg), call. = FALSE)
  }
  if (!is.na(msg)) {
    warning(sprintf("%s: %d of %d perturbation fits failed (first error: %s).",
                    lab, sum(!is.finite(ev$contrib)), length(ev$contrib), msg),
            call. = FALSE)
  }
  list(gdf = ev$gdf, gdf_se = ev$gdf_se, deviance = ev$deviance, gcv = ev$gcv,
       ubre = ev$ubre, n_eff = st$n_eff, contrib = ev$contrib, M_ok = ev$M_ok,
       converged = ev$converged, last_rel_change = ev$last_rel_change,
       budget = ev$budget, n_iter = ev$n_iter, gdf_cv = ev$gdf_cv,
       stable = isTRUE(ev$stable), reevals = ev$reevals %||% 0L, fit = ev$fit,
       epsilon = st$eps, n_fits = ev$n_fits,
       time_s = proc.time()[["elapsed"]] - t0)
}

# ---------------------------------------------------------------------------
# Sobol starts (optional n_global starts of the grouped stage)
# ---------------------------------------------------------------------------

#' Joe-Kuo direction-number tables for Sobol dims 2..8 (dim 1 = van der Corput).
#' @keywords internal
#' @noRd
.ggcv_sobol_params <- function() {
  # Each entry: s = degree, a = polynomial coefficients bitmask, m = initial m_i
  list(
    list(s = 1L, a = 0L, m = 1L),
    list(s = 2L, a = 1L, m = c(1L, 3L)),
    list(s = 3L, a = 1L, m = c(1L, 3L, 1L)),
    list(s = 3L, a = 2L, m = c(1L, 1L, 1L)),
    list(s = 4L, a = 1L, m = c(1L, 1L, 3L, 3L)),
    list(s = 4L, a = 4L, m = c(1L, 3L, 5L, 13L)),
    list(s = 5L, a = 2L, m = c(1L, 1L, 5L, 5L, 17L)),
    list(s = 5L, a = 4L, m = c(1L, 1L, 5L, 5L, 5L))
  )
}

#' Sobol points in (0,1)^d (skip origin). Builtin; no Suggests dependency.
#' @keywords internal
#' @noRd
.ggcv_sobol_unit <- function(n, d, skip = 1L) {
  n <- as.integer(n)
  d <- as.integer(d)
  skip <- as.integer(skip)
  stopifnot(n >= 1L, d >= 1L, d <= 8L, skip >= 0L)
  maxbit <- 30L
  params <- .ggcv_sobol_params()
  V <- matrix(0, nrow = d, ncol = maxbit)

  # Dimension 1: van der Corput base-2 direction numbers
  for (j in seq_len(maxbit)) {
    V[1L, j] <- bitwShiftL(1L, maxbit - j)
  }

  if (d >= 2L) {
    for (dim in 2:d) {
      p <- params[[dim - 1L]]
      s <- p$s
      a <- p$a
      m <- as.integer(p$m)
      stopifnot(length(m) >= s)
      mm <- integer(maxbit)
      mm[seq_len(s)] <- m[seq_len(s)]
      for (k in (s + 1L):maxbit) {
        mm[k] <- bitwXor(mm[k - s], bitwShiftL(mm[k - s], s))
        for (i in seq_len(s - 1L)) {
          if (bitwAnd(bitwShiftR(a, s - 1L - i), 1L) == 1L) {
            mm[k] <- bitwXor(mm[k], bitwShiftL(mm[k - i], i))
          }
        }
      }
      for (j in seq_len(maxbit)) {
        V[dim, j] <- bitwShiftL(mm[j], maxbit - j)
      }
    }
  }

  out <- matrix(NA_real_, nrow = n, ncol = d)
  # Gray-code Sobol
  X <- integer(d)
  # advance skip points
  for (i in seq_len(skip)) {
    c <- 1L
    value <- i
    while (bitwAnd(value, 1L) == 0L) {
      value <- bitwShiftR(value, 1L)
      c <- c + 1L
    }
    for (k in seq_len(d)) {
      X[k] <- bitwXor(X[k], V[k, c])
    }
  }
  for (i in seq_len(n)) {
    idx <- skip + i
    c <- 1L
    value <- idx
    while (bitwAnd(value, 1L) == 0L) {
      value <- bitwShiftR(value, 1L)
      c <- c + 1L
    }
    for (k in seq_len(d)) {
      X[k] <- bitwXor(X[k], V[k, c])
      out[i, k] <- X[k] / (2^maxbit)
    }
  }
  # Keep off exact 0/1 for log-lambda safety if ever mapped without padding
  out <- pmin(pmax(out, .Machine$double.eps), 1 - .Machine$double.eps)
  out
}

#' Map unit cube Sobol to a hyper-rectangle.
#' @keywords internal
#' @noRd
.ggcv_sobol_box <- function(n, lower, upper, skip = 1L) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  d <- length(lower)
  stopifnot(length(upper) == d, all(upper > lower))
  U <- .ggcv_sobol_unit(n, d, skip = skip)
  sweep(U, 2L, upper - lower, `*`) +
    matrix(lower, nrow = n, ncol = d, byrow = TRUE)
}

# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------

#' Isotropic-first gGCV search with paired-SE moves, the noise-aware final
#' choice and a decision-aware budget check.
#'
#' See the header of this file. `lo`, `hi`: box in log10(lambda), scalar or
#' one value per group. Returns the selector value of [tt_ggcv_array()].
#' @keywords internal
#' @noRd
.ggcv_search <- function(st, lo, hi, n_grid = 9L, n_refine = 2L,
                         n_global = 0L, n_final = 3L, M_search = 4L,
                         M_final = 16L, tol = 1, theta_start = NULL,
                         budget_check = TRUE, budget_tol = 0.01,
                         budget_fallback = TRUE, verbose = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  st$verbose <- isTRUE(verbose)
  G <- st$n_groups
  lab <- st$label
  lo <- rep(as.numeric(lo), length.out = G)
  hi <- rep(as.numeric(hi), length.out = G)
  if (anyNA(lo) || anyNA(hi) || any(hi <= lo)) {
    stop("theta_upper must exceed theta_lower.", call. = FALSE)
  }
  M_search <- as.integer(M_search)
  if (length(M_search) != 1L || is.na(M_search) || M_search < 1L) {
    stop("`M_search` must be a positive integer.", call. = FALSE)
  }
  M_final <- max(as.integer(M_final), M_search)
  n_refine <- max(0L, as.integer(n_refine))
  n_final <- max(1L, as.integer(n_final))
  n_global <- max(0L, as.integer(n_global %||% 0L))
  tol <- as.numeric(tol)
  starts <- if (is.null(theta_start)) NULL else {
    m <- as.matrix(theta_start)
    if (ncol(m) != G && nrow(m) == G) m <- t(m)
    if (ncol(m) != G) stop("theta_start must have one column per group.", call. = FALSE)
    m
  }
  clamp <- function(th) pmin(pmax(as.numeric(th), lo), hi)
  sdev <- function(recs) vapply(recs, function(r) r$score_dev, numeric(1))
  stable_of <- function(recs) vapply(recs, function(r) isTRUE(r$stable), logical(1))
  # best record by score_dev among the stable ones (with `fallback`, among
  # all of them when none is stable); NA when there is none
  best_of <- function(recs, fallback = FALSE) {
    s <- sdev(recs)
    ok <- stable_of(recs) & is.finite(s)
    if (!any(ok) && fallback) ok <- is.finite(s)
    if (!any(ok)) return(NA_integer_)
    which(ok)[which.min(s[ok])]
  }
  # b improves on a: b has a stable GDF and beats a under the paired rule
  improves <- function(a, b) isTRUE(b$stable) && .ggcv_paired(a, b, tol)$improves
  pool <- list()  # search-fidelity evaluations (final-stage candidates)

  # ---- 1. isotropic stage ------------------------------------------------
  ilo <- max(lo)
  ihi <- min(hi)
  if (!(ihi > ilo)) {
    stop(sprintf(paste0("%s: the isotropic stage needs a common range, ",
                        "max(theta_lower) < min(theta_upper)."), lab), call. = FALSE)
  }
  ng <- max(3L, as.integer(n_grid))
  step <- (ihi - ilo) / (ng - 1L)
  iso_t <- numeric(0)
  iso_r <- list()
  add_iso <- function(t, stage) {
    t <- pmin(pmax(as.numeric(t), ilo), ihi)
    recs <- .ggcv_eval_batch(matrix(rep(t, G), ncol = G), M_search, st, stage)
    iso_t <<- c(iso_t, t)
    iso_r <<- c(iso_r, recs)
    pool <<- c(pool, recs)
  }
  # the isotropic optimum: unstable evaluations only when nothing else exists
  best_iso <- function() best_of(iso_r, fallback = TRUE)
  t_grid <- seq(ilo, ihi, length.out = ng)
  stg <- rep("iso_grid", ng)
  if (G == 1L && !is.null(starts)) {
    t_grid <- c(t_grid, starts[, 1L])
    stg <- c(stg, rep("start", nrow(starts)))
  }
  add_iso(t_grid, stg)
  if (!any(is.finite(sdev(iso_r)))) {
    stop(sprintf("%s: no finite criterion value in the isotropic stage%s.", lab,
                 .ggcv_error_note(st)), call. = FALSE)
  }
  for (level in seq_len(n_refine)) {
    h <- step / 2^level
    tc <- iso_t[best_iso()] + c(-h, h)
    tc <- tc[tc >= ilo - 1e-9 & tc <= ihi + 1e-9]
    tc <- tc[vapply(tc, function(t) min(abs(t - iso_t)) > 1e-6, logical(1))]
    if (length(tc)) add_iso(tc, "iso_refine")
  }
  if (n_refine > 0L) {
    # one parabolic step through the best point and its nearest stable
    # neighbours
    h_last <- step / 2^n_refine
    b <- best_iso()
    tb <- iso_t[b]
    usable <- stable_of(iso_r) & is.finite(sdev(iso_r))
    left <- which(iso_t < tb - 1e-9 & usable)
    right <- which(iso_t > tb + 1e-9 & usable)
    if (isTRUE(usable[b]) && length(left) && length(right)) {
      l <- left[which.max(iso_t[left])]
      r <- right[which.min(iso_t[right])]
      tl <- iso_t[l]
      tr <- iso_t[r]
      s3 <- sdev(iso_r[c(l, b, r)])
      sl <- s3[1L]
      sb <- s3[2L]
      sr <- s3[3L]
      curv <- ((sr - sb) / (tr - tb) - (sb - sl) / (tb - tl)) / (tr - tl)
      if (is.finite(curv) && curv > 0) {
        v <- tb - 0.5 * ((tb - tl)^2 * (sb - sr) - (tb - tr)^2 * (sb - sl)) /
          ((tb - tl) * (sb - sr) - (tb - tr) * (sb - sl))
        v <- min(max(v, tb - h_last, ilo), tb + h_last, ihi)
        if (is.finite(v) && min(abs(v - iso_t)) > 1e-3) add_iso(v, "iso_parabola")
      }
    }
  }
  rec_iso <- iso_r[[best_iso()]]
  if (isTRUE(st$verbose)) {
    message(sprintf("%s isotropic optimum: theta=%.3f %s=%.6g [%.2f]", st$prefix,
                    rec_iso$theta[1L], st$criterion, rec_iso$score, rec_iso$score_dev))
  }

  # ---- 2. grouped stage --------------------------------------------------
  # Moves go only to evaluations with a stable GDF.
  cur <- rec_iso  # theta_cur = rep(theta_iso, G)
  if (G > 1L) {
    sp <- if (is.null(starts)) NULL else t(apply(starts, 1L, clamp))
    if (n_global > 0L) sp <- rbind(sp, .ggcv_sobol_box(n_global, lo, hi, skip = 1L))
    if (!is.null(sp) && nrow(sp)) {
      recs <- .ggcv_eval_batch(sp, M_search, st, "start")
      pool <- c(pool, recs)
      j <- best_of(recs)
      if (!is.na(j) && improves(cur, recs[[j]])) cur <- recs[[j]]
    }
    s <- 0.5
    passes <- 0L
    while (passes < 4L) {
      passes <- passes + 1L
      accepted <- FALSE
      for (g in seq_len(G)) {
        e <- replace(numeric(G), g, 1)
        pts <- rbind(clamp(cur$theta - s * e), clamp(cur$theta + s * e))
        dirs <- c(-1, 1)
        keep <- which(apply(pts, 1L, function(p) max(abs(p - cur$theta)) > 1e-9))
        if (!length(keep)) next
        recs <- .ggcv_eval_batch(pts[keep, , drop = FALSE], M_search, st, "pattern")
        pool <- c(pool, recs)
        j <- best_of(recs)
        if (is.na(j) || !improves(cur, recs[[j]])) next
        cur <- recs[[j]]
        accepted <- TRUE
        # keep going in the accepted direction with doubling steps
        dir <- dirs[keep[j]]
        len <- s
        repeat {
          len <- 2 * len
          nxt <- clamp(cur$theta + dir * len * e)
          if (max(abs(nxt - cur$theta)) <= 1e-9) break
          rec <- .ggcv_eval_batch(matrix(nxt, nrow = 1L), M_search, st, "pattern")[[1L]]
          pool <- c(pool, list(rec))
          if (!improves(cur, rec)) break
          cur <- rec
        }
      }
      if (!accepted) {
        if (s > 0.25) s <- 0.25 else break
      }
    }
    if (isTRUE(st$verbose)) {
      message(sprintf("%s grouped optimum: theta=(%s) %s=%.6g [%.2f]", st$prefix,
                      paste(sprintf("%.3f", cur$theta), collapse = ","),
                      st$criterion, cur$score, cur$score_dev))
    }
  }

  # ---- 3. final stage ----------------------------------------------------
  # Candidates: the isotropic optimum, the grouped optimum, the nearest
  # stable isotropic evaluation on the smoother side (larger theta), so that
  # the choice can move to a more regular point that the minimum does not
  # beat, and the next best stable search points more than 0.25 apart
  # (n_final in all; one group: the smoother neighbour comes on top of
  # n_final).
  is_iso <- function(th) diff(range(th)) < 1e-9
  in_list <- function(th) {
    any(vapply(cands, function(cj) max(abs(cj - th)) < 1e-9, logical(1)))
  }
  cands <- list(rec_iso$theta)
  labels <- "isotropic"
  if (G > 1L && !in_list(cur$theta)) {
    cands <- c(cands, list(cur$theta))
    labels <- c(labels, "grouped")
  }
  t_iso <- rec_iso$theta[1L]
  smoother <- NULL
  for (r in pool) {
    if (!is.finite(r$score_dev) || !isTRUE(r$stable) || !is_iso(r$theta) ||
        r$theta[1L] <= t_iso + 1e-9) next
    if (is.null(smoother) || r$theta[1L] < smoother[1L]) smoother <- r$theta
  }
  add_smoother <- function() {
    if (!is.null(smoother) && !in_list(smoother)) {
      cands <<- c(cands, list(smoother))
      labels <<- c(labels, "smoother")
    }
  }
  if (G > 1L) add_smoother()
  n_alt <- 0L
  for (i in order(sdev(pool))) {
    if (length(cands) >= n_final || !is.finite(pool[[i]]$score_dev)) break
    if (!isTRUE(pool[[i]]$stable)) next
    ti <- pool[[i]]$theta
    far <- all(vapply(cands, function(cj) sqrt(sum((cj - ti)^2)) > 0.25, logical(1)))
    if (far) {
      n_alt <- n_alt + 1L
      cands <- c(cands, list(ti))
      labels <- c(labels, paste0("alt", n_alt))
    }
  }
  if (G == 1L) add_smoother()
  fin <- .ggcv_eval_batch(do.call(rbind, cands), M_final, st, "final", keep_fit = TRUE)
  sdf <- sdev(fin)
  if (!any(is.finite(sdf))) {
    stop(sprintf("%s: no finite criterion value at final fidelity%s.", lab,
                 .ggcv_error_note(st)), call. = FALSE)
  }
  for (q in seq_along(fin)) fin[[q]]$iso <- is_iso(fin[[q]]$theta)
  conv <- vapply(fin, function(r) isTRUE(r$converged), logical(1))
  stab <- stable_of(fin)
  ch1 <- .ggcv_choose(fin, tol)
  none_conv <- isTRUE(ch1$none_converged)
  rule <- sprintf(paste0(
    "most regular final candidate that the minimum does not beat by more ",
    "than max(tol = %g, 2 * paired SE) in score_dev units (%s; isotropic ",
    "before grouped, larger theta first; converged candidates with a stable ",
    "GDF%s)"), tol,
    if (identical(st$criterion, "ubre")) "deviance / scale + 2 GDF" else "n log GCV",
    if (isTRUE(budget_check) && isTRUE(budget_fallback)) {
      paste0("; a candidate whose GDF fails the budget check is rescored with ",
             "the less favourable of its scores at B and 2B")
    } else "")

  # ---- 4. budget check: same probes at twice the winner's budget ---------
  # A GDF that moves by more than budget_tol is not the GDF of the converged
  # map. The first winner is always checked (the diagnostic of the returned
  # fit). With budget_fallback, a converged and stable winner that fails is
  # rescored with the less favourable of its B and 2B scores and the choice
  # is redone with that conservative score: a candidate that still wins is
  # kept (its GDF change cannot reverse the choice); otherwise the new winner
  # is checked (at most max_checks checks). No fallback check is spent on an
  # unconverged or unstable candidate: the check could not verify its map
  # (and an unconverged first winner means no converged candidate is left).
  # The returned candidate is always a checked one, converged and stable
  # whenever such a candidate exists: when the redone choice returns no
  # converged, stable candidate (.ggcv_choose() falls back to the others
  # once every converged, stable one has a non-finite rescored score), or
  # an unchecked one at the limit of checks, the first winner is kept
  # (`kept_first`, budget_verified = FALSE) and the first choice stands.
  max_checks <- 3L
  fin_c <- fin  # the candidates as scored in the choice
  ch <- ch1
  budget_ok <- rep(NA, length(fin))
  checks <- list()
  limit_hit <- FALSE
  kept_first <- NA_character_  # "ineligible" or "limit": first winner kept
  if (isTRUE(budget_check)) {
    repeat {
      w <- ch$winner
      if (length(checks)) {
        # the choice redone after a failed check of a converged, stable
        # winner: move only to a converged, stable candidate that is checked
        # or can still be checked
        if (is.na(w) || !conv[w] || !stab[w]) {
          kept_first <- "ineligible"
        } else if (is.na(budget_ok[w]) && length(checks) >= max_checks) {
          limit_hit <- TRUE
          kept_first <- "limit"
        }
        if (!is.na(kept_first)) {
          ch <- ch1
          break
        }
      }
      if (is.na(w) || !is.na(budget_ok[w])) break
      tb <- proc.time()[["elapsed"]]
      ev2 <- .ggcv_eval_batch(matrix(fin[[w]]$theta, nrow = 1L), M_final, st,
                              "budget", budget_mult = 2L)[[1L]]
      rel <- (fin[[w]]$gdf - ev2$gdf) / ev2$gdf
      ok <- is.finite(rel) && abs(rel) <= budget_tol
      budget_ok[w] <- ok
      checks[[length(checks) + 1L]] <- list(w = as.integer(w), ev2 = ev2, rel = rel,
                                            ok = ok,
                                            time_s = proc.time()[["elapsed"]] - tb)
      if (isTRUE(st$verbose)) {
        message(sprintf("%s budget check: %s theta=(%s) GDF %.3f at B=%d, %.3f at B=%d (%s)%s",
                        st$prefix, labels[w],
                        paste(sprintf("%.3f", fin[[w]]$theta), collapse = ","),
                        fin[[w]]$gdf, fin[[w]]$budget, ev2$gdf, ev2$budget,
                        .ggcv_pct_txt(rel), if (ok) "" else " failed"))
      }
      if (ok || !isTRUE(budget_fallback) || !conv[w] || !stab[w]) break
      # the less favourable of the B and 2B evaluations (same probes)
      if (!is.finite(ev2$score_dev) || ev2$score_dev > fin[[w]]$score_dev) {
        flds <- c("score_dev", "g", "contrib", "gdf", "gdf_se")
        fin_c[[w]][flds] <- ev2[flds]
      }
      ch <- .ggcv_choose(fin_c, tol)
    }
  }
  w <- ch$winner
  ok_checks <- vapply(checks, function(z) z$ok, logical(1))
  verified <- if (length(checks)) isTRUE(budget_ok[w]) else NA
  ev <- fin[[w]]
  theta <- ev$theta
  decision <- if (isTRUE(ev$iso)) "isotropic" else "grouped"
  chk <- Filter(function(z) z$w == w, checks)
  budget <- list(checked = FALSE, gdf_2x = NA_real_, rel_change = NA_real_,
                 score_2x = NA_real_, time_s = 0, mode = st$budget,
                 B = ev$budget, B_2x = NA_integer_, n_iter = ev$n_iter)
  if (length(chk)) {
    z <- chk[[1L]]
    budget[c("checked", "gdf_2x", "rel_change", "score_2x", "time_s", "B_2x")] <-
      list(TRUE, z$ev2$gdf, z$rel, z$ev2$score, z$time_s, z$ev2$budget)
  }
  th_of <- function(idx) {
    m <- do.call(rbind, lapply(idx, function(q) fin[[q]]$theta))
    colnames(m) <- paste0("theta", seq_len(G))
    m
  }
  budget_checks <- if (length(checks)) {
    ic <- vapply(checks, function(z) z$w, integer(1))
    data.frame(
      candidate = labels[ic], th_of(ic),
      gdf = vapply(ic, function(q) fin[[q]]$gdf, numeric(1)),
      gdf_2x = vapply(checks, function(z) z$ev2$gdf, numeric(1)),
      rel_change = vapply(checks, function(z) z$rel, numeric(1)),
      score_dev = vapply(ic, function(q) fin[[q]]$score_dev, numeric(1)),
      score_dev_2x = vapply(checks, function(z) z$ev2$score_dev, numeric(1)),
      ok = ok_checks, stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
  pr <- lapply(fin, function(r) .ggcv_paired(r, ev, tol))
  idx <- seq_along(fin)
  paired <- data.frame(
    candidate = labels, th_of(idx), score_dev = sdf,
    score_dev_used = sdev(fin_c),
    gdf = vapply(fin, function(r) r$gdf, numeric(1)),
    gdf_se = vapply(fin, function(r) r$gdf_se, numeric(1)),
    gdf_cv = vapply(fin, function(r) r$gdf_cv, numeric(1)),
    stable = stab, converged = conv, diff_vs_best = ch$diff,
    se_vs_best = ch$se, qualifies = ch$qualifies, best = idx == ch$best,
    winner = idx == w, budget_ok = budget_ok,
    diff_vs_winner = vapply(pr, function(p) p$diff, numeric(1)),
    se_vs_winner = vapply(pr, function(p) p$se, numeric(1)),
    stringsAsFactors = FALSE
  )
  if (isTRUE(st$verbose)) {
    message(sprintf("%s winner: %s theta=(%s) %s=%.6g [%.2f] decision=%s (minimum: %s)",
                    st$prefix, labels[w], paste(sprintf("%.3f", theta), collapse = ","),
                    st$criterion, ev$score, ev$score_dev, decision, labels[ch$best]))
  }

  # ---- 5. boundary, convergence, stability and budget warnings -----------
  boundary <- any(abs(theta - lo) < 0.05 | abs(theta - hi) < 0.05)
  if (isTRUE(boundary) && isTRUE(st$control_user$warn_lambda_boundary %||% TRUE)) {
    warning(lab, ": selected log10(lambda) is at the search boundary ",
            "(", paste(sprintf("%.2f", theta), collapse = ","), ").", call. = FALSE)
  }
  # one warning for one cause: an unconverged winner, an unstable winner, a
  # budget-sensitive winner that was kept, and a fall back past candidates
  # that failed the budget check. A GDF change that is not finite (no usable
  # GDF at 2B) reads "not finite", never "NA%".
  adaptive <- identical(st$budget, "adaptive")
  check_txt <- function(zs) {
    paste(vapply(zs, function(z) {
      sprintf("%s at theta=(%s)", .ggcv_pct_txt(z$rel),
              paste(sprintf("%.2f", fin[[z$w]]$theta), collapse = ","))
    }, character(1)), collapse = ", ")
  }
  failed <- checks[!ok_checks]
  n_fail <- length(failed)
  others <- Filter(function(z) z$w != w, failed)
  own_txt <- if (isTRUE(budget$checked)) {
    paste0("; ", .ggcv_change_txt(budget$rel_change, "the budget"))
  } else ""
  others_txt <- if (length(others)) {
    sprintf("; %d other final candidate(s) failed the budget check (GDF %s)",
            length(others), check_txt(others))
  } else ""
  advice <- if (adaptive) {
    paste0("Lower fit_tol, or raise control$pirls_maxit (Poisson) / ",
           "control$max_sweeps (Gaussian) if the cap was reached")
  } else {
    "Increase control$pirls_maxit (Poisson) or control$max_sweeps (Gaussian)"
  }
  unstable <- !isTRUE(ev$stable)
  cv_txt <- if (is.finite(ev$gdf_cv)) sprintf("%.2f", ev$gdf_cv) else "undefined"
  if (!isTRUE(ev$converged)) {
    reason <- if (adaptive) {
      rc <- .ggcv_rel_changes(.ggcv_objective_path(ev$fit))
      sprintf(paste0("the penalized objective did not settle to fit_tol = %.2g ",
                     "within the cap of %d iterations; relative change %s at ",
                     "the last one"), st$fit_tol, st$cap,
              .ggcv_num_txt(if (length(rc)) rc[length(rc)] else NA_real_))
    } else {
      sprintf("relative change %s at the last iteration, tolerance %.2g",
              .ggcv_num_txt(ev$last_rel_change), st$conv_tol)
    }
    warning(sprintf(paste0(
      "%s: the fit at the selected lambda is not converged at the iteration ",
      "budget (%s%s%s%s%s). ",
      "Increase control$pirls_maxit (Poisson) or control$max_sweeps (Gaussian)."),
      lab, reason,
      if (none_conv) "; no final candidate converged" else "",
      own_txt, others_txt,
      if (unstable) {
        sprintf(paste0("; its Monte Carlo GDF is unstable (per-probe CV %s) and ",
                       "no stable final candidate is eligible"), cv_txt)
      } else ""
    ), call. = FALSE)
  } else if (unstable) {
    warning(sprintf(paste0(
      "%s: Monte Carlo GDF is unstable at the selected lambda (per-probe CV %s, ",
      "or an outlying or non-finite probe contribution) and no stable final ",
      "candidate is eligible%s%s; increase M_final or the iteration budget."), lab,
      cv_txt, own_txt, others_txt), call. = FALSE)
  } else if (isFALSE(budget_ok[w])) {
    # checked, failed and kept: rescored with its less favourable score it
    # still wins, the fallback is off, or no converged, stable candidate
    # could replace it (kept_first: the first winner is kept)
    remains <- paste0("; rescored with the less favourable of its scores at B ",
                      "and 2B it remains the choice")
    margin <- if (identical(kept_first, "ineligible")) {
      sprintf(paste0("; rescored with the less favourable of its scores at B ",
                     "and 2B (score_dev %s) it would leave no converged final ",
                     "candidate with a stable GDF, and it is kept as the only ",
                     "eligible final candidate"),
              .ggcv_num_txt(fin_c[[w]]$score_dev, "%.2f"))
    } else if (identical(kept_first, "limit")) {
      sprintf(paste0("; the limit of %d checks was reached with every checked ",
                     "candidate failing, and the first winner is kept"),
              max_checks)
    } else if (isTRUE(budget_fallback)) {
      # the runner-up is the rule's winner without w, named only when it is
      # eligible (converged, with a stable GDF); .ggcv_choose() would fall
      # back to an unconverged or unstable candidate
      r <- .ggcv_choose(fin_c, tol, exclude = w)$winner
      if (is.na(r) || !conv[r] || !stab[r]) {
        sprintf(paste0("%s (score_dev %s) as the only eligible final candidate ",
                       "(converged, with a stable GDF)"), remains,
                .ggcv_num_txt(fin_c[[w]]$score_dev, "%.2f"))
      } else {
        sprintf("%s: score_dev %s, runner-up %s at %s (paired SE %s)", remains,
                .ggcv_num_txt(fin_c[[w]]$score_dev, "%.2f"), labels[r],
                .ggcv_num_txt(fin_c[[r]]$score_dev, "%.2f"),
                .ggcv_num_txt(.ggcv_paired(fin_c[[r]], fin_c[[w]], tol)$se, "%.2f"))
      }
    } else ""
    warning(sprintf(paste0(
      "%s: %s at the selected lambda; the fixed-lambda map is not converged ",
      "there%s%s, so budget_verified = FALSE. %s."), lab,
      .ggcv_change_txt(budget$rel_change, "the iteration budget"), margin,
      others_txt, advice), call. = FALSE)
  } else if (n_fail > 0L && isTRUE(verified)) {
    warning(sprintf(paste0(
      "%s: the fixed-lambda map is not converged at %d final candidate(s) (GDF %s ",
      "when the iteration budget doubles); rescored with the less favourable of ",
      "their scores at B and 2B they lose, and the selection moved to the next ",
      "candidate by the same rule, whose GDF is budget-verified."), lab, n_fail,
      check_txt(failed)), call. = FALSE)
  } else if (n_fail > 0L) {
    # not reached by the current loop (the returned candidate is always a
    # checked one); kept so that an unchecked winner is never silent
    warning(sprintf(paste0(
      "%s: the fixed-lambda map is not converged at %d final candidate(s) (GDF ",
      "%s when the iteration budget doubles)%s; the returned candidate's GDF ",
      "was not checked, so budget_verified = FALSE. %s."), lab, n_fail,
      check_txt(failed),
      if (limit_hit) {
        sprintf(" and the limit of %d checks was reached", max_checks)
      } else "",
      advice), call. = FALSE)
  }

  # ---- value ---------------------------------------------------------------
  search <- .ggcv_log_table(st)
  n_fits <- sum(search$n_fits)
  n_evals <- sum(search$stage != "budget")
  best_theta <- fin[[ch$best]]$theta
  fit <- ev$fit
  fit$lambda_method <- "gGCV"
  fit$ggcv <- list(
    mode = st$map$mode, criterion = st$criterion, score = ev$score,
    score_dev = ev$score_dev, gcv = ev$gcv, ubre = ev$ubre, gdf = ev$gdf,
    gdf_se = ev$gdf_se, gdf_mc_se = ev$gdf_se, gdf_cv = ev$gdf_cv,
    stable = ev$stable, deviance = ev$deviance, n_eff = st$n_eff,
    theta = theta, groups = st$groups, decision = decision, rule = rule,
    best_theta = best_theta, boundary = boundary, M = M_final,
    probe_init = st$probe_init, epsilon = st$eps, budget = budget,
    budget_checks = budget_checks, budget_verified = verified,
    converged = ev$converged, last_rel_change = ev$last_rel_change,
    paired = paired, search = search
  )
  list(
    lambda = ev$lambda, theta = theta, score = ev$score,
    score_dev = ev$score_dev, gcv = ev$gcv, ubre = ev$ubre, gdf = ev$gdf,
    gdf_se = ev$gdf_se, gdf_cv = ev$gdf_cv, stable = ev$stable,
    deviance = ev$deviance, n_eff = st$n_eff, criterion = st$criterion,
    groups = st$groups, decision = decision, rule = rule,
    best_theta = best_theta, boundary = boundary, budget = budget,
    budget_checks = budget_checks, budget_verified = verified,
    paired = paired, fit = fit, search = search, n_evals = n_evals,
    n_fits = n_fits, method = "gGCV", elapsed = proc.time()[["elapsed"]] - t0
  )
}
