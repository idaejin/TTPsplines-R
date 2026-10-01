# Global TT-gGCV for array (product-grid) data, ttps(array = TRUE).
# Fitting regimes (naming per 2026-08-19 DECISION): K-ALS for unweighted
# Gaussian (Kronecker Gram); weighted GLAM row-tensor PIRLS for Poisson or
# observation weights. The selector treats both as one fixed-lambda map.
#
# Criterion at a fixed smoothing vector lambda:
#   gGCV(lambda) = n * D(y, mu_hat) / (n - GDF)^2,
#   GDF = tr(d mu_hat / d y) of the fixed-lambda fitting map
#         y |-> ttps(Y, array = TRUE, lambda = lambda),
# with n = number of cells carrying positive weight. Optional UBRE form
# D / n + 2 * scale * GDF / n - scale (Poisson: scale = 1 known).
# criterion = "auto" (default) follows mgcv's GCV.Cp rule: UBRE when the
# scale is known (Poisson), GCV when it is estimated (Gaussian). GCV on a
# Poisson deviance estimates the scale from D / n, which goes wrong when the
# response is under- or near-degenerate (STATS19 occupied-only pilot).
#
# GDF is estimated by Hutchinson / Monte Carlo finite differences with
# NON-NEGATIVE perturbations, so Poisson responses stay valid on zero-heavy
# counts (no clipping at 0):
#   r ~ Rademacher on the weighted cells, z = (r + 1) / 2 in {0, 1},
#   J r = 2 J (z - 1/2) ~= (2 / eps) * (mu(y + eps z) - mu(y + eps / 2)),
#   E[r' J r] = tr J.
# One shared reference fit at y + eps / 2 replaces the base fit in the
# difference, so each probe costs one fit and the base fit's own
# convergence error cancels.
#
# Determinism: every fit (estimator, reference, probes) runs the SAME fixed
# iteration budget (control$tol is forced to 0, so neither PIRLS nor ALS stops
# early) from the SAME initialization ("cold", default). The GDF is then the
# derivative of one smooth map, which is the regime in which Hutchinson GDF
# tracks the joint EDF (manuscript appendix on cGCV versus gGCV).
# probe_init = "warm" restarts probes from the estimator's cores with a
# shorter budget; cheaper, but it under-responds when ALS mixes slowly.

#' Global TT-gGCV smoothing selection for array data (TT-GLAM, `array = TRUE`)
#'
#' Selects \eqn{\lambda} for a product-grid TT P-spline fitted with
#' [ttps()]`(array = TRUE)` by minimising the global GCV score
#' \eqn{n D(\lambda) / (n - \mathrm{GDF}(\lambda))^2}, where GDF is the
#' Monte Carlo trace of the derivative of the fixed-\eqn{\lambda} array fit:
#' K-ALS for unweighted Gaussian data, weighted GLAM row-tensor PIRLS for
#' Poisson responses or observation weights. Offset arrays are supported.
#'
#' @section Perturbations:
#' Probes are non-negative (\eqn{y + \epsilon (r+1)/2} against a shared
#' reference \eqn{y + \epsilon/2}), so Poisson counts are never clipped at
#' zero and each probe costs one refit. Probes are drawn on cells with
#' positive weight only and are identical across \eqn{\lambda} (common random
#' numbers); the first `M_search` probes of the final bank equal the search
#' bank.
#'
#' @section Iteration budget:
#' `control$tol` is set to `0`, so every fit runs exactly
#' `control$max_sweeps` ALS sweeps (Gaussian) or `control$pirls_maxit` PIRLS
#' iterations (Poisson). Choose a budget at which the fixed-\eqn{\lambda} fit
#' has converged; the returned `fit$ggcv$last_rel_change` reports the relative
#' deviance change of the last iteration when available.
#'
#' @param Y d-way response array (d >= 2).
#' @param axes List of grid coordinates, one per margin (default: unit grids).
#' @param family `gaussian()` or `poisson()`.
#' @param rank,k,degree,penalty_order,cyclic,period,knots Passed to [ttps()].
#' @param weights Optional non-negative weight array (same dim as `Y`).
#' @param offset Optional offset array (same dim as `Y`, linear-predictor scale).
#' @param groups Integer vector of length d mapping margins to smoothing
#'   groups that share one \eqn{\lambda}. `NULL` = isotropic (one group);
#'   `seq_len(d)` = fully anisotropic.
#' @param theta_lower,theta_upper Search box in \eqn{\log_{10}\lambda}
#'   (scalar or one value per group; default from `control$lambda_bounds`).
#' @param theta_start Optional matrix (or vector) of extra starting points in
#'   group space (\eqn{\log_{10}\lambda}), evaluated in the global stage.
#' @param n_grid Coarse grid size for the isotropic (one-group) search.
#' @param n_global Sobol points for the multi-group search
#'   (default `8 * n_groups`).
#' @param n_refine Brent refinement (one group) or coordinate-descent passes
#'   (several groups). `0` disables refinement.
#' @param n_final Number of distinct best candidates re-scored with `M_final`
#'   probes; the winner is chosen at final fidelity only.
#' @param M_search,M_final Rademacher probes in the search / final stages.
#' @param criterion `"auto"` (default): UBRE with known scale for Poisson
#'   and GCV for Gaussian, as `mgcv`'s `"GCV.Cp"`; or force `"gcv"` /
#'   `"ubre"`. For overdispersed counts use `"gcv"` or pass an estimated
#'   `scale` with `"ubre"`.
#' @param scale Scale for UBRE. Poisson default `1`; required for Gaussian
#'   UBRE.
#' @param probe_init `"cold"` (default) or `"warm"`.
#' @param probe_budget Iteration budget for warm probes (PIRLS iterations for
#'   Poisson, ALS sweeps for Gaussian); ignored when `probe_init = "cold"`.
#' @param epsilon_rel Finite-difference step relative to the RMS of `Y`
#'   (floored at `epsilon_rel` itself).
#' @param n_cores Parallel workers over probes (`parallel::mclapply`; forked
#'   processes, so `1` on Windows).
#' @param seed Probe seed.
#' @param control [tt_control()] for every fixed-\eqn{\lambda} fit.
#' @param budget_check If `TRUE` (default), re-estimate the winner's GDF with
#'   twice the iteration budget and the same probes. A relative change above
#'   `budget_tol` means the fixed-\eqn{\lambda} map is not converged at the
#'   current budget (seen for low-rank Poisson fits at 10 PIRLS iterations);
#'   a warning suggests raising `control$pirls_maxit` / `control$max_sweeps`.
#'   The doubling change is a lower bound on the distance to the converged
#'   map when ALS converges slowly (2 to 4 PIRLS iterations moved GDF by 1.5%
#'   while the converged value was 10% away), hence the tight default.
#' @param budget_tol Relative GDF change that triggers the warning (`0.01`).
#' @param verbose Print one line per evaluation.
#' @return A list with `lambda` (length d), `theta` (group scale), `score`,
#'   `gcv`, `ubre`, `gdf`, `gdf_se`, `deviance`, `n_eff`, `criterion`,
#'   `boundary`, `budget` (GDF at twice the budget and its relative change),
#'   `fit` (the [ttps()] fit at the winner, with `$ggcv`), `search`
#'   (all evaluations), `n_fits` and `elapsed`.
#' @seealso [tt_gdf_array()] to score a given \eqn{\lambda}.
#' @export
tt_ggcv_array <- function(Y,
                          axes = NULL,
                          family = stats::gaussian(),
                          rank = 3L,
                          k = 8L,
                          degree = 3L,
                          penalty_order = 2L,
                          cyclic = NULL,
                          period = NULL,
                          knots = NULL,
                          weights = NULL,
                          offset = NULL,
                          groups = NULL,
                          theta_lower = NULL,
                          theta_upper = NULL,
                          theta_start = NULL,
                          n_grid = 9L,
                          n_global = NULL,
                          n_refine = 2L,
                          n_final = 3L,
                          M_search = 4L,
                          M_final = 16L,
                          criterion = c("auto", "gcv", "ubre"),
                          scale = NULL,
                          probe_init = c("cold", "warm"),
                          probe_budget = NULL,
                          epsilon_rel = 1e-3,
                          n_cores = 1L,
                          seed = 1L,
                          control = tt_control(max_sweeps = 10L,
                                               pirls_maxit = 20L,
                                               compute_edf = FALSE,
                                               seed = 1L),
                          budget_check = TRUE,
                          budget_tol = 0.01,
                          verbose = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  criterion <- match.arg(criterion)
  probe_init <- match.arg(probe_init)
  st <- .ggcv_arr_setup(
    Y = Y, axes = axes, family = family, rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, cyclic = cyclic, period = period,
    knots = knots, weights = weights, offset = offset, groups = groups,
    criterion = criterion, scale = scale, probe_init = probe_init,
    probe_budget = probe_budget, epsilon_rel = epsilon_rel,
    n_cores = n_cores, seed = seed, control = control
  )
  criterion <- st$criterion
  G <- st$n_groups
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  lo <- rep(as.numeric(theta_lower %||% log10(bounds[1L])), length.out = G)
  hi <- rep(as.numeric(theta_upper %||% log10(bounds[2L])), length.out = G)
  if (any(hi <= lo)) stop("theta_upper must exceed theta_lower.", call. = FALSE)
  M_search <- as.integer(M_search)
  M_final <- max(as.integer(M_final), M_search)

  cache <- new.env(parent = emptyenv())
  evals <- list()
  evaluate <- function(theta, M, stage) {
    theta <- pmin(pmax(as.numeric(theta), lo), hi)
    key <- paste(c(sprintf("%.6f", theta), M), collapse = "|")
    if (!is.null(cache[[key]])) return(cache[[key]])
    ev <- .ggcv_arr_eval(theta, M = M, st = st)
    row <- data.frame(
      stage = stage, t(setNames(theta, paste0("theta", seq_len(G)))),
      score = ev$score, gcv = ev$gcv, ubre = ev$ubre, gdf = ev$gdf,
      gdf_se = ev$gdf_se, deviance = ev$deviance, M = M, M_ok = ev$M_ok,
      time_s = ev$time_s, stringsAsFactors = FALSE
    )
    evals[[length(evals) + 1L]] <<- row
    if (isTRUE(verbose)) {
      message(sprintf(
        "gGCV-array %-7s theta=(%s) %s=%.6g GDF=%.2f (se %.2f) dev=%.6g [%.1fs]",
        stage, paste(sprintf("%.2f", theta), collapse = ","), criterion,
        ev$score, ev$gdf, ev$gdf_se, ev$deviance, ev$time_s
      ))
    }
    cache[[key]] <- ev
    ev
  }

  # ---- global stage --------------------------------------------------------
  starts <- if (is.null(theta_start)) NULL else {
    m <- as.matrix(theta_start)
    if (ncol(m) != G && nrow(m) == G) m <- t(m)
    if (ncol(m) != G) stop("theta_start must have one column per group.", call. = FALSE)
    m
  }
  if (G == 1L) {
    pts <- matrix(seq(lo, hi, length.out = max(3L, as.integer(n_grid))), ncol = 1L)
  } else {
    ng <- as.integer(n_global %||% (8L * G))
    pts <- rbind((lo + hi) / 2, .tt_lab_sobol_box(ng, lo, hi, skip = 1L))
  }
  pts <- rbind(starts, pts)
  for (i in seq_len(nrow(pts))) evaluate(pts[i, ], M_search, "global")

  best_theta <- function() {
    tab <- do.call(rbind, evals)
    tab <- tab[tab$M == M_search & is.finite(tab$score), , drop = FALSE]
    if (!nrow(tab)) return(NULL)
    as.numeric(tab[which.min(tab$score), paste0("theta", seq_len(G))])
  }
  th <- best_theta()
  if (is.null(th)) {
    stop("tt_ggcv_array: no finite criterion value in the global stage.",
         call. = FALSE)
  }

  # ---- refinement ----------------------------------------------------------
  if (as.integer(n_refine) > 0L) {
    if (G == 1L) {
      step <- if (nrow(pts) > 1L) diff(range(pts[, 1L])) / (nrow(pts) - 1L) else 1
      iv <- c(max(lo, th - step), min(hi, th + step))
      tryCatch(stats::optimize(function(t) {
        s <- evaluate(t, M_search, "refine")$score
        if (is.finite(s)) s else .Machine$double.xmax
      }, interval = iv, tol = 0.02), error = function(e) NULL)
    } else {
      for (pass in seq_len(as.integer(n_refine))) {
        for (g in seq_len(G)) {
          iv <- c(max(lo[g], th[g] - 1.5), min(hi[g], th[g] + 1.5))
          f1 <- function(t) {
            tt <- th; tt[g] <- t
            s <- evaluate(tt, M_search, "refine")$score
            if (is.finite(s)) s else .Machine$double.xmax
          }
          tryCatch(stats::optimize(f1, interval = iv, tol = 0.05),
                   error = function(e) NULL)
          th <- best_theta()
        }
      }
    }
    th <- best_theta()
  }

  # ---- final stage: re-score distinct best candidates with M_final ---------
  tab <- do.call(rbind, evals)
  tab <- tab[tab$M == M_search & is.finite(tab$score), , drop = FALSE]
  tab <- tab[order(tab$score), , drop = FALSE]
  cand <- list()
  for (i in seq_len(nrow(tab))) {
    ti <- as.numeric(tab[i, paste0("theta", seq_len(G))])
    far <- all(vapply(cand, function(cj) sqrt(sum((cj - ti)^2)) > 0.25, logical(1)))
    if (!length(cand) || far) cand[[length(cand) + 1L]] <- ti
    if (length(cand) >= as.integer(n_final)) break
  }
  fin <- lapply(cand, function(ti) list(theta = ti, ev = evaluate(ti, M_final, "final")))
  sc <- vapply(fin, function(z) z$ev$score, numeric(1))
  if (!any(is.finite(sc))) {
    stop("tt_ggcv_array: no finite criterion value at final fidelity.", call. = FALSE)
  }
  win <- fin[[which.min(ifelse(is.finite(sc), sc, Inf))]]
  theta <- win$theta
  ev <- win$ev
  lambda <- 10^theta[st$groups]
  boundary <- any(abs(theta - lo) < 0.05 | abs(theta - hi) < 0.05)
  if (isTRUE(boundary) && isTRUE(control$warn_lambda_boundary %||% TRUE)) {
    warning("tt_ggcv_array: selected log10(lambda) is at the search boundary ",
            "(", paste(sprintf("%.2f", theta), collapse = ","), ").", call. = FALSE)
  }
  budget <- list(checked = FALSE, gdf_2x = NA_real_, rel_change = NA_real_)
  if (isTRUE(budget_check)) {
    st2 <- st
    for (nm in c("control", "control_probe")) {
      st2[[nm]]$pirls_maxit <- 2L * as.integer(st[[nm]]$pirls_maxit)
      st2[[nm]]$max_sweeps <- 2L * as.integer(st[[nm]]$max_sweeps)
    }
    ev2 <- .ggcv_arr_eval(theta, M = M_final, st = st2)  # same probes (CRN)
    rel <- (ev$gdf - ev2$gdf) / ev2$gdf
    budget <- list(checked = TRUE, gdf_2x = ev2$gdf, rel_change = rel,
                   score_2x = ev2$score, time_s = ev2$time_s)
    if (isTRUE(verbose)) {
      message(sprintf("gGCV-array budget check: GDF %.3f at budget, %.3f at 2x (%+.1f%%)",
                      ev$gdf, ev2$gdf, 100 * rel))
    }
    if (is.finite(rel) && abs(rel) > budget_tol) {
      warning(sprintf(paste0(
        "tt_ggcv_array: GDF changes by %+.1f%% when the iteration budget doubles; ",
        "the fixed-lambda fit is not converged. Increase control$pirls_maxit ",
        "(Poisson) or control$max_sweeps (Gaussian)."), 100 * rel), call. = FALSE)
    }
  }
  search <- do.call(rbind, evals)
  n_fits <- sum(search$M + 2L)  # estimator + reference + M probes per evaluation
  if (isTRUE(budget$checked)) n_fits <- n_fits + M_final + 2L
  fit <- ev$fit
  fit$lambda_method <- "gGCV"
  fit$ggcv <- list(
    mode = "array", criterion = criterion, score = ev$score, gcv = ev$gcv,
    ubre = ev$ubre, gdf = ev$gdf, gdf_se = ev$gdf_se, gdf_mc_se = ev$gdf_se,
    deviance = ev$deviance, n_eff = st$n_eff, theta = theta,
    groups = st$groups, boundary = boundary, M = M_final,
    probe_init = probe_init, epsilon = st$eps, budget = budget,
    last_rel_change = .ggcv_arr_last_rel_change(fit), search = search
  )
  list(
    lambda = lambda, theta = theta, score = ev$score, gcv = ev$gcv,
    ubre = ev$ubre, gdf = ev$gdf, gdf_se = ev$gdf_se, deviance = ev$deviance,
    n_eff = st$n_eff, criterion = criterion, groups = st$groups,
    boundary = boundary, budget = budget, fit = fit, search = search,
    n_fits = n_fits,
    method = "gGCV", elapsed = proc.time()[["elapsed"]] - t0
  )
}

#' Global GDF and gGCV score of an array TT fit at a fixed lambda
#'
#' Scores one smoothing vector with the same estimator as [tt_ggcv_array()]:
#' the fixed-\eqn{\lambda} array fit (K-ALS or weighted PIRLS), its deviance,
#' and a Monte Carlo
#' GDF from non-negative Rademacher perturbations.
#'
#' @inheritParams tt_ggcv_array
#' @param lambda Smoothing vector (scalar or length d).
#' @param M Number of probes.
#' @param probes `"rademacher"` (default) or `"unit"`. `"unit"` perturbs one
#'   weighted cell at a time and returns the exact finite-difference trace
#'   (one fit per cell; small arrays only).
#' @return A list with `gdf`, `gdf_se`, `deviance`, `gcv`, `ubre`, `n_eff`,
#'   `contrib`, `fit`, `epsilon`.
#' @export
tt_gdf_array <- function(Y,
                         lambda,
                         axes = NULL,
                         family = stats::gaussian(),
                         rank = 3L,
                         k = 8L,
                         degree = 3L,
                         penalty_order = 2L,
                         cyclic = NULL,
                         period = NULL,
                         knots = NULL,
                         weights = NULL,
                         offset = NULL,
                         M = 16L,
                         probes = c("rademacher", "unit"),
                         scale = NULL,
                         probe_init = c("cold", "warm"),
                         probe_budget = NULL,
                         epsilon_rel = 1e-3,
                         n_cores = 1L,
                         seed = 1L,
                         control = tt_control(max_sweeps = 10L,
                                              pirls_maxit = 20L,
                                              compute_edf = FALSE,
                                              seed = 1L)) {
  probes <- match.arg(probes)
  probe_init <- match.arg(probe_init)
  d <- length(dim(Y))
  st <- .ggcv_arr_setup(
    Y = Y, axes = axes, family = family, rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, cyclic = cyclic, period = period,
    knots = knots, weights = weights, offset = offset, groups = seq_len(d),
    criterion = "gcv", scale = scale, probe_init = probe_init,
    probe_budget = probe_budget, epsilon_rel = epsilon_rel,
    n_cores = n_cores, seed = seed, control = control
  )
  theta <- log10(rep(as.numeric(lambda), length.out = d))
  ev <- .ggcv_arr_eval(theta, M = as.integer(M), st = st, unit = identical(probes, "unit"))
  list(gdf = ev$gdf, gdf_se = ev$gdf_se, deviance = ev$deviance, gcv = ev$gcv,
       ubre = ev$ubre, n_eff = st$n_eff, contrib = ev$contrib, M_ok = ev$M_ok,
       fit = ev$fit, epsilon = st$eps, time_s = ev$time_s)
}

# ---------------------------------------------------------------------------
# internals
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.ggcv_arr_setup <- function(Y, axes, family, rank, k, degree, penalty_order,
                            cyclic, period, knots, weights, offset, groups,
                            criterion, scale, probe_init, probe_budget,
                            epsilon_rel, n_cores, seed, control) {
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
  groups <- if (is.null(groups)) rep(1L, d) else as.integer(groups)
  if (length(groups) != d || anyNA(groups)) {
    stop("`groups` must be an integer vector of length d.", call. = FALSE)
  }
  groups <- match(groups, sort(unique(groups)))
  if (identical(criterion, "auto")) {
    criterion <- if (identical(key, "poisson")) "ubre" else "gcv"
  }
  if (is.null(scale) && identical(key, "poisson")) scale <- 1
  if (identical(criterion, "ubre") && is.null(scale)) {
    stop("criterion = 'ubre' needs a known `scale` for gaussian().", call. = FALSE)
  }
  ctrl <- control
  if (!inherits(ctrl, "tt_control")) ctrl <- do.call(tt_control, as.list(ctrl))
  ctrl$tol <- 0
  ctrl$compute_edf <- FALSE
  ctrl$trace <- FALSE
  ctrl$monitor <- FALSE
  ctrl$warn_lambda_boundary <- FALSE
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
  rms <- sqrt(sum(y[support]^2) / n_eff)
  eps <- as.numeric(epsilon_rel) * if (is.finite(rms) && rms > 0) max(rms, 1) else 1
  n_cores <- max(1L, as.integer(n_cores))
  if (n_cores > 1L && (.Platform$OS.type == "windows" ||
                       !requireNamespace("parallel", quietly = TRUE))) {
    n_cores <- 1L
  }
  list(
    dims = dims, d = d, n_cells = n_cells, y = y, fam = fam, key = key,
    axes = axes, rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, cyclic = cyclic, period = period,
    knots = knots, weights = w, offset = off, support = support,
    n_eff = n_eff, groups = groups, n_groups = max(groups),
    criterion = criterion, scale = scale, probe_init = probe_init,
    eps = eps, n_cores = n_cores, seed = as.integer(seed),
    control = ctrl, control_probe = ctrl_probe
  )
}

#' One fixed-lambda array fit.
#' @keywords internal
#' @noRd
.ggcv_arr_fit <- function(yvec, lambda, st, init = NULL, ctrl = st$control) {
  # ttps() reseeds the global RNG for its core initialization; keep the
  # caller's stream intact.
  .tt_with_preserved_seed(ttps(
    array(yvec, dim = st$dims), axes = st$axes, array = TRUE, family = st$fam,
    rank = st$rank, k = st$k, degree = st$degree,
    penalty_order = st$penalty_order, lambda = lambda, optimizer = "ALS",
    init = init, control = ctrl, knots = st$knots, offset = st$offset,
    weights = st$weights, cyclic = st$cyclic, period = st$period
  ))
}

#' Rademacher probe j on the weighted cells (common random numbers).
#' @keywords internal
#' @noRd
.ggcv_arr_probe <- function(j, st) {
  .tt_with_preserved_seed({
    set.seed(st$seed * 100003L + as.integer(j))
    r <- numeric(st$n_cells)
    r[st$support] <- sample(c(-1, 1), st$n_eff, replace = TRUE)
  })
  r
}

#' Score one theta (group scale): estimator fit, deviance, Monte Carlo GDF.
#' @keywords internal
#' @noRd
.ggcv_arr_eval <- function(theta, M, st, unit = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  lambda <- 10^as.numeric(theta)[st$groups]
  fit0 <- tryCatch(.ggcv_arr_fit(st$y, lambda, st), error = function(e) e)
  if (inherits(fit0, "error")) {
    return(list(score = Inf, gcv = Inf, ubre = Inf, gdf = NA_real_,
                gdf_se = NA_real_, deviance = NA_real_, M_ok = 0L,
                contrib = numeric(0), fit = NULL, time_s = 0,
                reason = conditionMessage(fit0)))
  }
  mu0 <- as.numeric(fit0$fitted.values)
  dev <- glm_deviance(st$fam, st$y, mu0, weights = st$weights)
  eps <- st$eps
  init <- if (identical(st$probe_init, "warm")) .tt_clone_cores(fit0$cores) else NULL
  ctrl_p <- st$control_probe
  refit_mu <- function(yy) {
    f <- tryCatch(.ggcv_arr_fit(yy, lambda, st, init = init, ctrl = ctrl_p),
                  error = function(e) NULL)
    if (is.null(f)) return(NULL)
    as.numeric(f$fitted.values)
  }
  apply_jobs <- function(jobs, FUN) {
    if (st$n_cores > 1L && length(jobs) > 1L) {
      parallel::mclapply(jobs, FUN, mc.cores = st$n_cores, mc.preschedule = FALSE)
    } else {
      lapply(jobs, FUN)
    }
  }
  if (isTRUE(unit)) {
    # exact finite-difference trace: one non-negative unit perturbation per cell
    base <- if (is.null(init)) mu0 else refit_mu(st$y)
    cells <- which(st$support)
    contrib <- unlist(apply_jobs(cells, function(i) {
      yy <- st$y
      yy[i] <- yy[i] + eps
      m <- refit_mu(yy)
      if (is.null(m)) NA_real_ else (m[i] - base[i]) / eps
    }))
    ok <- is.finite(contrib)
    gdf <- if (all(ok)) sum(contrib) else NA_real_
    gdf_se <- 0
  } else {
    jobs <- c(0L, seq_len(M))
    mus <- apply_jobs(jobs, function(j) {
      if (j == 0L) {
        yy <- st$y + 0.5 * eps * as.numeric(st$support)
      } else {
        r <- .ggcv_arr_probe(j, st)
        yy <- st$y + eps * (r + 1) / 2 * as.numeric(st$support)
      }
      refit_mu(yy)
    })
    mu_ref <- mus[[1L]]
    contrib <- rep(NA_real_, M)
    if (!is.null(mu_ref)) {
      for (j in seq_len(M)) {
        m <- mus[[j + 1L]]
        if (is.null(m)) next
        r <- .ggcv_arr_probe(j, st)
        contrib[j] <- sum(r * (m - mu_ref)) * 2 / eps
      }
    }
    ok <- is.finite(contrib)
    gdf <- if (any(ok)) mean(contrib[ok]) else NA_real_
    gdf_se <- if (sum(ok) >= 2L) stats::sd(contrib[ok]) / sqrt(sum(ok)) else NA_real_
  }
  n <- st$n_eff
  valid <- is.finite(gdf) && gdf > 0 && gdf < n && is.finite(dev)
  gcv <- if (valid) n * dev / (n - gdf)^2 else Inf
  ubre <- if (valid && !is.null(st$scale)) {
    dev / n + 2 * st$scale * gdf / n - st$scale
  } else if (valid) NA_real_ else Inf
  score <- if (identical(st$criterion, "ubre")) ubre else gcv
  list(score = score, gcv = gcv, ubre = ubre, gdf = gdf, gdf_se = gdf_se,
       deviance = dev, M_ok = sum(ok), contrib = contrib, fit = fit0,
       time_s = proc.time()[["elapsed"]] - t0)
}

#' Relative deviance change at the last iteration of a fit (if recorded).
#' @keywords internal
#' @noRd
.ggcv_arr_last_rel_change <- function(fit) {
  h <- fit$history %||% fit$convergence$history %||% NULL
  if (!is.data.frame(h) || nrow(h) < 2L) return(NA_real_)
  col <- intersect(c("deviance", "rss", "objective"), names(h))
  if (!length(col)) return(NA_real_)
  v <- as.numeric(h[[col[1L]]])
  v <- v[is.finite(v)]
  if (length(v) < 2L) return(NA_real_)
  abs(v[length(v)] - v[length(v) - 1L]) / max(1, abs(v[length(v) - 1L]))
}

#' Resolve `lambda = "gGCV"` inside [ttps()] when `array = TRUE`.
#'
#' Search knobs are read from `control` under array-specific names (the
#' scattered `ggcv_M_*` defaults of [tt_control()] are far too expensive on a
#' full grid): `ggcv_array_groups`, `ggcv_array_n_global`,
#' `ggcv_array_n_refine`, `ggcv_array_M_search`, `ggcv_array_M_final`,
#' `ggcv_array_criterion`, `ggcv_array_probe_init`, `ggcv_array_n_cores`,
#' `ggcv_array_budget_check`.
#' Set them on a [tt_control()] object with `$<-`.
#' @keywords internal
#' @noRd
.ttps_dispatch_ggcv_array <- function(Y, axes, family, rank, k, degree,
                                      penalty_order, cyclic, period, knots,
                                      weights, offset, control, cl = NULL) {
  w <- if (is.null(weights) || all(abs(as.numeric(weights) - 1) < 1e-12)) {
    NULL
  } else {
    weights
  }
  off <- if (is.null(offset) || all(as.numeric(offset) == 0)) NULL else offset
  opt <- tt_ggcv_array(
    Y = Y, axes = axes, family = family, rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, cyclic = cyclic, period = period,
    knots = knots, weights = w, offset = off,
    groups = control$ggcv_array_groups %||% NULL,
    n_global = control$ggcv_array_n_global %||% NULL,
    n_refine = as.integer(control$ggcv_array_n_refine %||% 2L),
    M_search = as.integer(control$ggcv_array_M_search %||% 4L),
    M_final = as.integer(control$ggcv_array_M_final %||% 16L),
    criterion = control$ggcv_array_criterion %||% "auto",
    probe_init = control$ggcv_array_probe_init %||% "cold",
    n_cores = as.integer(control$ggcv_array_n_cores %||% 1L),
    seed = as.integer(control$seed %||% 1L),
    control = control,
    budget_check = !isFALSE(control$ggcv_array_budget_check),
    verbose = isTRUE(control$trace)
  )
  fit <- opt$fit
  if (!is.null(cl)) fit$call <- cl
  fit
}
