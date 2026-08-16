# Experimental adaptive joint λ optimizer (v1).
#
# NOT exported. Does not replace operational cGCV (KEEP cGCV).
# Builds on v0: adaptive box expansion, near-optimal regions, λ
# identifiability labels, batched Sobol, MC-tie rule.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Near-optimal region from evaluated θ points (1+α relative to Q_min).
#' @keywords internal
#' @noRd
.tt_lab_near_optimal_region <- function(theta, q, alpha = 0.01) {
  theta <- as.matrix(theta)
  q <- as.numeric(q)
  ok <- is.finite(q) & apply(theta, 1L, function(r) all(is.finite(r)))
  if (!any(ok)) {
    return(list(
      mask = ok, q_min = NA_real_, alpha = alpha,
      intervals = NULL, n_in = 0L, theta_in = NULL, q_in = NULL
    ))
  }
  q_min <- min(q[ok])
  mask <- ok & (q <= (1 + alpha) * q_min)
  th_in <- theta[mask, , drop = FALSE]
  d <- ncol(theta)
  intervals <- matrix(NA_real_, nrow = d, ncol = 2L,
                      dimnames = list(paste0("theta", seq_len(d)),
                                      c("lo", "hi")))
  if (nrow(th_in)) {
    for (j in seq_len(d)) {
      intervals[j, ] <- range(th_in[, j], na.rm = TRUE)
    }
  }
  list(
    mask = mask,
    q_min = q_min,
    alpha = alpha,
    intervals = intervals,
    n_in = sum(mask),
    theta_in = th_in,
    q_in = q[mask]
  )
}

#' Classify per-margin identifiability from profile + region interval.
#' @keywords internal
#' @noRd
.tt_lab_classify_lambda_id <- function(intervals, lower, upper,
                                      profile_list,
                                      boundary_tol = 0.5,
                                      narrow_width = 1.0,
                                      flat_rel = 0.005) {
  d <- length(lower)
  labels <- character(d)
  details <- vector("list", d)
  for (j in seq_len(d)) {
    I <- intervals[j, ]
    width <- I[2] - I[1]
    touches_lo <- is.finite(I[1]) && (I[1] - lower[j] < boundary_tol)
    touches_hi <- is.finite(I[2]) && (upper[j] - I[2] < boundary_tol)
    prof <- profile_list[[j]]
    flat <- FALSE
    flat_lo_tail <- FALSE
    rises_both <- FALSE
    edge_still_improving <- FALSE
    if (!is.null(prof) && nrow(prof) >= 3L && any(is.finite(prof$q))) {
      qv <- prof$q
      tv <- prof$theta_j
      q_rng <- diff(range(qv, na.rm = TRUE))
      q_med <- stats::median(qv, na.rm = TRUE)
      flat <- is.finite(q_rng) && q_rng <= flat_rel * max(abs(q_med), 1e-8)
      # Flatness on the near-optimal interval / low-λ tail (not the whole box)
      in_I <- is.finite(qv) & is.finite(I[1]) & is.finite(I[2]) &
        tv >= I[1] - 1e-9 & tv <= I[2] + 1e-9
      if (sum(in_I) >= 2L) {
        q_I <- qv[in_I]
        flat_I <- diff(range(q_I, na.rm = TRUE)) <=
          flat_rel * max(abs(stats::median(q_I, na.rm = TRUE)), 1e-8)
      } else {
        flat_I <- FALSE
      }
      # Low-λ plateau: among θ_j ≤ max(I_hi, -4), Q barely changes
      hi_cut <- if (is.finite(I[2])) max(I[2], -4) else -4
      lo_tail <- is.finite(qv) & tv <= hi_cut + 1e-9
      if (sum(lo_tail) >= 2L) {
        q_lo <- qv[lo_tail]
        flat_lo_tail <- diff(range(q_lo, na.rm = TRUE)) <=
          flat_rel * max(abs(stats::median(q_lo, na.rm = TRUE)), 1e-8)
      }
      flat <- flat || flat_I || flat_lo_tail
      # rise both sides of interval midpoint
      mid <- mean(I)
      left <- qv[tv < mid - 1e-9]
      right <- qv[tv > mid + 1e-9]
      center <- qv[tv >= mid - 0.25 & tv <= mid + 0.25]
      c0 <- if (length(center) && any(is.finite(center))) {
        min(center, na.rm = TRUE)
      } else {
        min(qv, na.rm = TRUE)
      }
      rises_both <- length(left) > 0L && length(right) > 0L &&
        mean(left, na.rm = TRUE) > c0 * (1 + flat_rel) &&
        mean(right, na.rm = TRUE) > c0 * (1 + flat_rel)
      # Still improving at the lower face? (true unresolved boundary)
      i_lo <- which.min(tv)
      i_next <- order(tv)[min(2L, length(tv))]
      if (length(i_lo) && is.finite(qv[i_lo]) && is.finite(qv[i_next])) {
        edge_still_improving <- qv[i_lo] < qv[i_next] * (1 - flat_rel)
      }
    }
    plateau_low <- flat_lo_tail || (flat && is.finite(I[2]) && I[2] < -3.5)
    if (plateau_low && touches_lo) {
      lab <- "effectively_unpenalized"
    } else if (plateau_low) {
      lab <- "weakly_identified"
    } else if ((touches_lo || touches_hi) && width > narrow_width &&
               edge_still_improving) {
      lab <- "boundary_unresolved"
    } else if (flat || width > narrow_width) {
      lab <- "weakly_identified"
    } else if (rises_both && width <= narrow_width && !touches_lo && !touches_hi) {
      lab <- "well_identified"
    } else if (width <= narrow_width && !touches_lo && !touches_hi) {
      lab <- "well_identified"
    } else if (touches_lo && width > narrow_width) {
      lab <- "weakly_identified"
    } else {
      lab <- "weakly_identified"
    }
    labels[j] <- lab
    details[[j]] <- list(
      width = width, touches_lo = touches_lo, touches_hi = touches_hi,
      flat = flat, flat_lo_tail = flat_lo_tail, rises_both = rises_both,
      edge_still_improving = edge_still_improving
    )
  }
  names(labels) <- paste0("lambda", seq_len(d))
  list(labels = labels, details = details)
}

#' 1-D profile of Q along margin j, other θ fixed at theta0.
#' @keywords internal
#' @noRd
.tt_lab_profile_margin <- function(evaluator, theta0, j, lower, upper,
                                   M, n_grid = 17L) {
  d <- length(theta0)
  j <- as.integer(j)
  grid_j <- seq(lower[j], upper[j], length.out = as.integer(n_grid))
  rows <- vector("list", length(grid_j))
  for (i in seq_along(grid_j)) {
    th <- theta0
    th[j] <- grid_j[i]
    ev <- evaluator$eval(th, M = M, mc_bank_id = "profile", return_fit = FALSE)
    rows[[i]] <- data.frame(
      margin = j, theta_j = grid_j[i], q = ev$global_gcv,
      gdf = ev$gdf, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' MC-tie: candidates indistinguishable if |Qa-Qb| < 2 sqrt(SEa^2+SEb^2).
#' Falls back to relative 1% if SE missing.
#' @keywords internal
#' @noRd
.tt_lab_gcv_tied <- function(qa, qb, sea = NA_real_, seb = NA_real_,
                            alpha = 0.01) {
  if (!is.finite(qa) || !is.finite(qb)) return(FALSE)
  if (is.finite(sea) && is.finite(seb) && (sea + seb) > 0) {
    return(abs(qa - qb) < 2 * sqrt(sea^2 + seb^2))
  }
  qmin <- min(qa, qb)
  abs(qa - qb) <= alpha * max(abs(qmin), 1e-12)
}

#' Edge probes: is Q at the box face still competitive with the best?
#'
#' Holding θ_{-j} at `theta_best`, evaluate Q at θ_j = L_j and U_j.
#' A face is "active" if its Q is ≤ (1+α) Q_best (flat / truncated valley).
#' @keywords internal
#' @noRd
.tt_lab_edge_probe <- function(evaluator, theta_best, lower, upper,
                               M = 12L, alpha = 0.01) {
  d <- length(theta_best)
  q_best <- evaluator$eval(theta_best, M = M, return_fit = FALSE)$global_gcv
  force_lo <- rep(FALSE, d)
  force_hi <- rep(FALSE, d)
  q_lo <- rep(NA_real_, d)
  q_hi <- rep(NA_real_, d)
  thr <- if (is.finite(q_best)) (1 + alpha) * q_best else Inf
  for (j in seq_len(d)) {
    th_lo <- theta_best
    th_lo[j] <- lower[j]
    th_hi <- theta_best
    th_hi[j] <- upper[j]
    ev_lo <- evaluator$eval(th_lo, M = M, return_fit = FALSE)
    ev_hi <- evaluator$eval(th_hi, M = M, return_fit = FALSE)
    q_lo[j] <- ev_lo$global_gcv
    q_hi[j] <- ev_hi$global_gcv
    force_lo[j] <- is.finite(q_lo[j]) && q_lo[j] <= thr
    force_hi[j] <- is.finite(q_hi[j]) && q_hi[j] <= thr
  }
  list(q_best = q_best, q_lo = q_lo, q_hi = q_hi,
       force_lo = force_lo, force_hi = force_hi)
}

#' Expand box dims that touch the boundary (directed).
#'
#' Expands margin j when:
#' * the **best** θ_j is within `soft_tol` of an edge, or
#' * the 1% region touches an edge and the best is within `2*soft_tol`, or
#' * an edge probe marks that face as still competitive (`force_lo`/`force_hi`).
#' @keywords internal
#' @noRd
.tt_lab_expand_box <- function(lower, upper, theta_best, region_intervals,
                               boundary_tol = 0.5,
                               soft_tol = 1.0,
                               expand_step = 3,
                               hard_lo = -12, hard_hi = 8,
                               force_lo = NULL,
                               force_hi = NULL) {
  d <- length(lower)
  lower2 <- lower
  upper2 <- upper
  expanded <- integer(0)
  soft_tol <- max(as.numeric(soft_tol), as.numeric(boundary_tol))
  if (is.null(force_lo)) force_lo <- rep(FALSE, d)
  if (is.null(force_hi)) force_hi <- rep(FALSE, d)
  for (j in seq_len(d)) {
    thb <- theta_best[j]
    I <- if (!is.null(region_intervals)) region_intervals[j, ] else c(thb, thb)
    region_lo <- is.finite(I[1]) && (I[1] - lower[j] <= boundary_tol)
    region_hi <- is.finite(I[2]) && (upper[j] - I[2] <= boundary_tol)
    hit_lo <- isTRUE(force_lo[j]) ||
      (thb - lower[j] <= boundary_tol) ||
      (region_lo && (thb - lower[j] <= soft_tol))
    hit_hi <- isTRUE(force_hi[j]) ||
      (upper[j] - thb <= boundary_tol) ||
      (region_hi && (upper[j] - thb <= soft_tol))
    if (hit_lo) {
      lower2[j] <- max(hard_lo, lower[j] - expand_step)
      if (lower2[j] < lower[j] - 1e-12) expanded <- c(expanded, j)
    }
    if (hit_hi) {
      upper2[j] <- min(hard_hi, upper[j] + expand_step)
      if (upper2[j] > upper[j] + 1e-12) expanded <- c(expanded, j)
    }
  }
  list(lower = lower2, upper = upper2,
       expanded_dimensions = unique(as.integer(expanded)))
}
#' Sobol points restricted to the newly added slab(s) after a box expand.
#' @keywords internal
#' @noRd
.tt_lab_sobol_new_slab <- function(n, lower_old, upper_old, lower_new, upper_new,
                                   skip = 1L) {
  d <- length(lower_new)
  # Sample in new box; keep points that fall outside the old box
  raw <- .tt_lab_sobol_box(max(n * 4L, n + 10L), lower_new, upper_new, skip = skip)
  keep <- apply(raw, 1L, function(th) {
    any(th < lower_old - 1e-12 | th > upper_old + 1e-12)
  })
  pts <- raw[keep, , drop = FALSE]
  if (!nrow(pts)) {
    # fallback: fill expanded margins with a line of points
    pts <- matrix(NA_real_, nrow = 0L, ncol = d)
  }
  if (nrow(pts) > n) pts <- pts[seq_len(n), , drop = FALSE]
  if (nrow(pts) < n) {
    # pad with Sobol in new box
    more <- .tt_lab_sobol_box(n - nrow(pts), lower_new, upper_new,
                             skip = skip + 1000L)
    pts <- rbind(pts, more)
  }
  pts[seq_len(min(n, nrow(pts))), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# v1 optimizer
# ---------------------------------------------------------------------------

#' Adaptive experimental joint λ optimizer (v1).
#'
#' Extends [tt_global_lambda_optimize()] with batched Sobol, directed box
#' expansion, near-optimal regions, and λ identifiability labels.
#' Not exported. Productive path remains `lambda = "cGCV"`.
#'
#' @inheritParams tt_global_lambda_optimize
#' @param adaptive_box Expand only margins that hit the boundary.
#' @param max_expansions Maximum directed expansions.
#' @param expand_boundary_tol δ for edge detection (default 0.5).
#' @param near_optimal_tol α for R_α (default 0.01).
#' @param sobol_batches Integer vector of incremental Sobol sizes.
#' @param profile_lambda If TRUE, build 1-D profiles for ID labels.
#' @param classify_identifiability If TRUE, label each λ.
#' @param expand_step How much to widen a bound when expanding.
#' @param stop_on_stable_region If TRUE, stop Sobol batches early when the
#'   near-optimal region intervals stabilize (default TRUE). Prefer FALSE
#'   for higher-d validation gates.
#' @param min_batches_before_stop Minimum batches before early stop can fire.
#' @return Diagnostic list including region and identifiability.
#' @keywords internal
tt_global_lambda_optimize_v1 <- function(y = NULL,
                                         X = NULL,
                                         formula = NULL,
                                         data = NULL,
                                         rank,
                                         family = stats::gaussian(),
                                         theta_lower = -5,
                                         theta_upper = 5,
                                         n_refine = 5L,
                                         n_diverse = 10L,
                                         M_search = 15L,
                                         M_final = 40L,
                                         core_starts_final = 3L,
                                         min_dist = 0.75,
                                         boundary_tol = 0.05,
                                         seed = 1L,
                                         k = 8L,
                                         degree = 3L,
                                         penalty_order = 2L,
                                         control = tt_control(max_sweeps = 40L,
                                                              tol = 1e-8,
                                                              compute_edf = FALSE,
                                                              seed = 1L),
                                         include_cgcv_anchor = TRUE,
                                         extra_theta = NULL,
                                         epsilon_rel = 1e-3,
                                         scheme = c("forward", "central"),
                                         adaptive_box = TRUE,
                                         max_expansions = 2L,
                                         expand_boundary_tol = 0.5,
                                         near_optimal_tol = 0.01,
                                         sobol_batches = c(64L, 64L, 128L),
                                         profile_lambda = TRUE,
                                         classify_identifiability = TRUE,
                                         expand_step = 3,
                                         stop_on_stable_region = TRUE,
                                         min_batches_before_stop = 2L,
                                         verbose = FALSE) {
  t_wall0 <- proc.time()[["elapsed"]]
  scheme <- match.arg(scheme)
  fam_key <- family_key(normalize_family(family))
  if (!identical(fam_key, "gaussian")) {
    stop("tt_global_lambda_optimize_v1: Gaussian only.", call. = FALSE)
  }

  if (is.null(y) || is.null(X)) {
    if (is.null(formula) || is.null(data)) {
      stop("Supply y+X or formula+data.", call. = FALSE)
    }
    mf <- stats::model.frame(formula, data = data)
    y <- stats::model.response(mf)
    X <- stats::model.matrix(formula, data = data)
    if ("(Intercept)" %in% colnames(X)) {
      X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
    }
  }
  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y)
  d <- ncol(X)
  stopifnot(nrow(X) == n)

  lower <- rep(as.numeric(theta_lower), length.out = d)
  upper <- rep(as.numeric(theta_upper), length.out = d)
  sobol_batches <- as.integer(sobol_batches)
  max_expansions <- as.integer(max_expansions)
  M_search <- as.integer(M_search)
  M_final <- as.integer(M_final)
  seed <- as.integer(seed)
  k <- as.integer(k)

  ctrl <- .tt_lab_refit_control(control)
  ctrl$seed <- seed
  common_init <- .tt_with_preserved_seed({
    tt_initialize(d = d, rank = rank, k = k, seed = seed)
  })
  M_bank <- max(M_search, M_final)
  probes_bank <- .tt_lab_rademacher_probes(n, M_bank, probe_seed = seed + 11L)
  probes_search <- probes_bank[, seq_len(M_search), drop = FALSE]
  probes_final <- probes_bank[, seq_len(M_final), drop = FALSE]
  probes_alt <- .tt_lab_rademacher_probes(n, M_final, probe_seed = seed + 911L)

  # Anchors once (cGCV + extras)
  cgcv <- list(ok = FALSE)
  if (isTRUE(include_cgcv_anchor)) {
    if (verbose) message("[v1] cGCV anchor...")
    cgcv <- .tt_lab_cgcv_anchor(
      y, X, rank, common_init, ctrl, k, degree, penalty_order
    )
  }
  anchors <- .tt_lab_theta_anchors(
    d, lower, upper,
    theta_cgcv = if (isTRUE(cgcv$ok)) cgcv$theta else NULL
  )
  an_names <- rownames(anchors)
  an_names[an_names == "cgcv"] <- "cgcv_anchor"
  an_names[an_names == "centre"] <- "center"
  rownames(anchors) <- an_names
  if (!is.null(extra_theta)) {
    et <- as.matrix(extra_theta)
    if (ncol(et) != d && nrow(et) == d) et <- t(et)
    rownames(et) <- paste0("prior_solution_", seq_len(nrow(et)))
    anchors <- rbind(anchors, et)
  }

  search_history <- list()
  all_theta <- matrix(numeric(0), nrow = 0L, ncol = d)
  all_q <- numeric(0)
  all_source <- character(0)
  all_gdf_se <- numeric(0)
  expanded_dims_all <- integer(0)
  sobol_skip <- 1L
  region_stable_streak <- 0L
  prev_intervals <- NULL

  eval_points <- function(theta_mat, source_lab, lower_use, upper_use) {
    evaluator <- .tt_lab_make_theta_evaluator(
      y = y, X = X, rank = rank, common_init = common_init,
      probes = probes_search, control = ctrl, k = k, degree = degree,
      penalty_order = penalty_order, epsilon_rel = epsilon_rel, scheme = scheme,
      lower = lower_use, upper = upper_use
    )
    n_pt <- nrow(theta_mat)
    th_out <- matrix(NA_real_, n_pt, d)
    q_out <- rep(Inf, n_pt)
    se_out <- rep(NA_real_, n_pt)
    src_out <- character(n_pt)
    for (i in seq_len(n_pt)) {
      th <- as.numeric(theta_mat[i, ])
      # project into current box for evaluation
      thc <- pmin(pmax(th, lower_use), upper_use)
      ev <- evaluator$eval(thc, M = M_search, mc_bank_id = "bank0",
                           return_fit = FALSE)
      th_out[i, ] <- thc
      q_out[i] <- ev$global_gcv
      se_out[i] <- ev$gdf_mc_se %||% NA_real_
      src_out[i] <- source_lab[i]
    }
    list(theta = th_out, q = q_out, gdf_se = se_out, source = src_out,
         cost = evaluator$cost)
  }

  # ---- Batched Sobol exploration ------------------------------------------
  batch_id <- 0L
  for (nb in sobol_batches) {
    batch_id <- batch_id + 1L
    if (verbose) {
      message(sprintf("[v1] Sobol batch %d: +%d (box lo=%s hi=%s)",
                      batch_id, nb,
                      paste(sprintf("%.1f", lower), collapse = ","),
                      paste(sprintf("%.1f", upper), collapse = ",")))
    }
    sob <- .tt_lab_sobol_box(nb, lower, upper, skip = sobol_skip)
    sobol_skip <- sobol_skip + nb
    src <- rep(paste0("sobol_batch", batch_id), nrow(sob))
    # include anchors only on first batch
    if (batch_id == 1L && nrow(anchors)) {
      sob <- rbind(sob, anchors)
      src <- c(src, rownames(anchors))
    }
    evb <- eval_points(sob, src, lower, upper)
    all_theta <- rbind(all_theta, evb$theta)
    all_q <- c(all_q, evb$q)
    all_source <- c(all_source, evb$source)
    all_gdf_se <- c(all_gdf_se, evb$gdf_se)

    reg <- .tt_lab_near_optimal_region(all_theta, all_q, near_optimal_tol)
    search_history[[batch_id]] <- list(
      type = "sobol_batch", n = nb, q_min = reg$q_min,
      n_in_region = reg$n_in, intervals = reg$intervals,
      lower = lower, upper = upper
    )
    # Stability of region intervals
    if (!is.null(prev_intervals) && !is.null(reg$intervals)) {
      delta_I <- max(abs(reg$intervals - prev_intervals), na.rm = TRUE)
      if (is.finite(delta_I) && delta_I < 0.25) {
        region_stable_streak <- region_stable_streak + 1L
      } else {
        region_stable_streak <- 0L
      }
    }
    prev_intervals <- reg$intervals
    min_b <- max(2L, as.integer(min_batches_before_stop))
    if (isTRUE(stop_on_stable_region) &&
        batch_id >= min_b && region_stable_streak >= 1L &&
        is.finite(reg$q_min)) {
      if (verbose) message("[v1] near-optimal region stable — stop batches")
      break
    }
  }

  # ---- Refine elites ------------------------------------------------------
  refine_once <- function(theta_mat, q_vec, source_vec, lower_use, upper_use) {
    elite <- .tt_lab_diverse_top(theta_mat, q_vec, n_keep = n_diverse,
                                min_dist = min_dist)
    if (!length(elite)) elite <- which.min(q_vec)
    starts <- elite[seq_len(min(as.integer(n_refine), length(elite)))]
    evaluator <- .tt_lab_make_theta_evaluator(
      y = y, X = X, rank = rank, common_init = common_init,
      probes = probes_search, control = ctrl, k = k, degree = degree,
      penalty_order = penalty_order, epsilon_rel = epsilon_rel, scheme = scheme,
      lower = lower_use, upper = upper_use
    )
    out_th <- matrix(NA_real_, length(starts), d)
    out_q <- rep(Inf, length(starts))
    out_src <- character(length(starts))
    for (ii in seq_along(starts)) {
      i0 <- starts[ii]
      th0 <- as.numeric(theta_mat[i0, ])
      obj <- function(th) {
        # Bypass cache so nlminb finite differences are not collapsed
        ev <- evaluator$eval(th, M = M_search, return_fit = FALSE,
                             use_cache = FALSE)
        if (!is.finite(ev$global_gcv)) 1e30 else ev$global_gcv
      }
      opt <- tryCatch(
        stats::nlminb(th0, obj, lower = lower_use, upper = upper_use,
                      control = list(iter.max = 60L)),
        error = function(e) list(par = th0)
      )
      th_r <- as.numeric(opt$par)
      evr <- evaluator$eval(th_r, M = M_search, return_fit = FALSE)
      out_th[ii, ] <- th_r
      out_q[ii] <- evr$global_gcv
      src0 <- if (length(source_vec) >= i0) source_vec[i0] else "unknown"
      out_src[ii] <- paste0("refined_from_", src0)
    }
    list(theta = out_th, q = out_q, source = out_src)
  }

  if (verbose) message("[v1] refining elites...")
  ref <- refine_once(all_theta, all_q, all_source, lower, upper)
  all_theta <- rbind(all_theta, ref$theta)
  all_q <- c(all_q, ref$q)
  all_source <- c(all_source, ref$source)
  all_gdf_se <- c(all_gdf_se, rep(NA_real_, length(ref$q)))

  # ---- Directed expansions ------------------------------------------------
  # Edge probes detect truncated flat valleys (Sobol best may sit interior
  # while Q at the face remains within the 1% band — e.g. strong_aniso θ2).
  n_exp <- 0L
  while (isTRUE(adaptive_box) && n_exp < max_expansions) {
    i_best <- which.min(all_q)
    th_best <- as.numeric(all_theta[i_best, ])
    reg <- .tt_lab_near_optimal_region(all_theta, all_q, near_optimal_tol)
    evaluator_edge <- .tt_lab_make_theta_evaluator(
      y = y, X = X, rank = rank, common_init = common_init,
      probes = probes_search, control = ctrl, k = k, degree = degree,
      penalty_order = penalty_order, epsilon_rel = epsilon_rel, scheme = scheme,
      lower = lower, upper = upper
    )
    edge <- .tt_lab_edge_probe(
      evaluator_edge, th_best, lower, upper,
      M = M_search, alpha = near_optimal_tol
    )
    if (verbose) {
      message(sprintf(
        "[v1] edge probe Q_best=%.6g  force_lo=%s  force_hi=%s",
        edge$q_best,
        paste(as.integer(edge$force_lo), collapse = ","),
        paste(as.integer(edge$force_hi), collapse = ",")
      ))
    }
    exp <- .tt_lab_expand_box(
      lower, upper, th_best, reg$intervals,
      boundary_tol = expand_boundary_tol,
      soft_tol = max(expand_boundary_tol, 1.0),
      expand_step = expand_step,
      force_lo = edge$force_lo,
      force_hi = edge$force_hi
    )
    if (!length(exp$expanded_dimensions)) break
    n_exp <- n_exp + 1L
    if (verbose) {
      message(sprintf("[v1] expansion %d on margins %s → lo=%s",
                      n_exp,
                      paste(exp$expanded_dimensions, collapse = ","),
                      paste(sprintf("%.1f", exp$lower), collapse = ",")))
    }
    lower_old <- lower
    upper_old <- upper
    lower <- exp$lower
    upper <- exp$upper
    expanded_dims_all <- unique(c(expanded_dims_all, exp$expanded_dimensions))

    # New Sobol only in added slab
    n_new <- max(64L, as.integer(sum(sobol_batches) / 2L))
    slab <- .tt_lab_sobol_new_slab(
      n_new, lower_old, upper_old, lower, upper, skip = sobol_skip
    )
    sobol_skip <- sobol_skip + n_new
    src <- rep(paste0("sobol_expand", n_exp), nrow(slab))
    evb <- eval_points(slab, src, lower, upper)
    all_theta <- rbind(all_theta, evb$theta)
    all_q <- c(all_q, evb$q)
    all_source <- c(all_source, evb$source)
    all_gdf_se <- c(all_gdf_se, evb$gdf_se)

    ref2 <- refine_once(all_theta, all_q, all_source, lower, upper)
    all_theta <- rbind(all_theta, ref2$theta)
    all_q <- c(all_q, ref2$q)
    all_source <- c(all_source, ref2$source)
    all_gdf_se <- c(all_gdf_se, rep(NA_real_, length(ref2$q)))

    search_history[[length(search_history) + 1L]] <- list(
      type = "expansion", expanded = exp$expanded_dimensions,
      lower = lower, upper = upper
    )
  }

  # ---- Near-optimal region + profiles -------------------------------------
  reg <- .tt_lab_near_optimal_region(all_theta, all_q, near_optimal_tol)
  i_best <- which.min(all_q)
  theta_best_search <- as.numeric(all_theta[i_best, ])
  winner_source_search <- all_source[i_best]

  profiles <- vector("list", d)
  if (isTRUE(profile_lambda)) {
    if (verbose) message("[v1] profiling margins...")
    evaluator_p <- .tt_lab_make_theta_evaluator(
      y = y, X = X, rank = rank, common_init = common_init,
      probes = probes_search, control = ctrl, k = k, degree = degree,
      penalty_order = penalty_order, epsilon_rel = epsilon_rel, scheme = scheme,
      lower = lower, upper = upper
    )
    for (j in seq_len(d)) {
      profiles[[j]] <- .tt_lab_profile_margin(
        evaluator_p, theta_best_search, j, lower, upper, M_search, n_grid = 15L
      )
      # fold profile points into region pool
      all_theta <- rbind(all_theta, {
        m <- matrix(rep(theta_best_search, nrow(profiles[[j]])),
                    nrow = nrow(profiles[[j]]), ncol = d, byrow = TRUE)
        m[, j] <- profiles[[j]]$theta_j
        m
      })
      all_q <- c(all_q, profiles[[j]]$q)
      all_source <- c(all_source, rep(paste0("profile_m", j), nrow(profiles[[j]])))
      all_gdf_se <- c(all_gdf_se, rep(NA_real_, nrow(profiles[[j]])))
    }
    reg <- .tt_lab_near_optimal_region(all_theta, all_q, near_optimal_tol)
  }

  id <- NULL
  if (isTRUE(classify_identifiability)) {
    id <- .tt_lab_classify_lambda_id(
      reg$intervals, lower, upper, profiles,
      boundary_tol = expand_boundary_tol
    )
  }

  # ---- Final re-eval with MC-tie rule -------------------------------------
  # Candidates: best search, diverse region points, cGCV
  cand_th <- list(theta_best_search)
  cand_src <- winner_source_search
  if (!is.null(reg$theta_in) && nrow(reg$theta_in)) {
    # up to 5 diverse points from region
    idx_r <- .tt_lab_diverse_top(reg$theta_in, reg$q_in, n_keep = 5L,
                                min_dist = min_dist)
    for (ii in idx_r) {
      cand_th[[length(cand_th) + 1L]] <- as.numeric(reg$theta_in[ii, ])
      cand_src <- c(cand_src, "region_point")
    }
  }
  if (isTRUE(cgcv$ok)) {
    cand_th[[length(cand_th) + 1L]] <- as.numeric(cgcv$theta)
    cand_src <- c(cand_src, "cgcv_anchor")
  }
  key <- vapply(cand_th, function(th) paste(sprintf("%.4f", th), collapse = ","),
                character(1))
  keep <- !duplicated(key)
  cand_th <- cand_th[keep]
  cand_src <- cand_src[keep]

  if (verbose) {
    message(sprintf("[v1] final re-eval %d candidates...", length(cand_th)))
  }
  ctrl_strict <- ctrl
  ctrl_strict$tol <- min(as.numeric(ctrl$tol %||% 1e-8), 1e-10)
  ctrl_strict$max_sweeps <- max(as.integer(ctrl$max_sweeps %||% 40L), 60L)

  final_rows <- vector("list", length(cand_th))
  final_fits <- vector("list", length(cand_th))
  for (j in seq_along(cand_th)) {
    re1 <- .tt_lab_reeval_multistart(
      theta = cand_th[[j]], y = y, X = X, rank = rank,
      probes = probes_final, M = M_final,
      n_starts = core_starts_final, seed = seed,
      control = ctrl_strict, k = k, degree = degree,
      penalty_order = penalty_order, epsilon_rel = epsilon_rel,
      scheme = scheme, lower = lower, upper = upper
    )
    re2 <- .tt_lab_reeval_multistart(
      theta = cand_th[[j]], y = y, X = X, rank = rank,
      probes = probes_alt, M = M_final,
      n_starts = max(1L, min(3L, core_starts_final)), seed = seed + 17L,
      control = ctrl_strict, k = k, degree = degree,
      penalty_order = penalty_order, epsilon_rel = epsilon_rel,
      scheme = scheme, lower = lower, upper = upper
    )
    final_fits[[j]] <- re1$fit
    # Rough SE proxy from GDF MC se via delta method not available; store gdf se
    final_rows[[j]] <- data.frame(
      cand = j, source = cand_src[j],
      matrix(cand_th[[j]], nrow = 1L,
             dimnames = list(NULL, paste0("theta", seq_len(d)))),
      gcv_bank0 = re1$global_gcv,
      gcv_bank1 = re2$global_gcv,
      gdf_se = re1$gdf_mc_se,
      basin_index = re1$basin_index,
      ok = re1$ok,
      stringsAsFactors = FALSE
    )
  }
  final_tab <- do.call(rbind, final_rows)
  # Tie-aware selection: all within tie of the best stay in region
  j0 <- which.min(final_tab$gcv_bank0)
  q0 <- final_tab$gcv_bank0[j0]
  se0 <- final_tab$gdf_se[j0]
  tied <- vapply(seq_len(nrow(final_tab)), function(j) {
    .tt_lab_gcv_tied(q0, final_tab$gcv_bank0[j], se0, final_tab$gdf_se[j],
                     near_optimal_tol)
  }, logical(1))
  outer_search_stable <- isTRUE(
    identical(which.min(final_tab$gcv_bank0), which.min(final_tab$gcv_bank1))
  ) || all(tied[c(which.min(final_tab$gcv_bank0), which.min(final_tab$gcv_bank1))])

  # Inner basin: majority of multistart basins agree on best cand (proxy)
  inner_basin_stable <- isTRUE(final_tab$ok[j0])

  theta_best <- as.numeric(final_tab[j0, paste0("theta", seq_len(d))])
  gcv_best <- as.numeric(final_tab$gcv_bank0[j0])
  best_fit <- final_fits[[j0]]
  winner_source <- as.character(final_tab$source[j0])

  near_lo <- abs(theta_best - lower) <= boundary_tol
  near_hi <- abs(theta_best - upper) <= boundary_tol
  # Also flag expand-tol boundary for reporting
  near_expand <- abs(theta_best - lower) <= expand_boundary_tol |
    abs(theta_best - upper) <= expand_boundary_tol

  elapsed <- proc.time()[["elapsed"]] - t_wall0
  list(
    theta = theta_best,
    lambda = 10^theta_best,
    gcv = gcv_best,
    fit = best_fit,
    near_optimal_region = list(
      alpha = near_optimal_tol,
      q_min = reg$q_min,
      n_points = reg$n_in,
      intervals = reg$intervals,
      theta = reg$theta_in,
      q = reg$q_in,
      tied_final_candidates = which(tied)
    ),
    profile_intervals = reg$intervals,
    profiles = profiles,
    lambda_identifiability = if (is.null(id)) NULL else id$labels,
    lambda_identifiability_details = if (is.null(id)) NULL else id$details,
    expanded_dimensions = as.integer(expanded_dims_all),
    search_history = search_history,
    outer_search_stable = outer_search_stable,
    inner_basin_stable = inner_basin_stable,
    winner_source = winner_source,
    boundary = any(near_lo | near_hi),
    boundary_expand_tol = any(near_expand),
    box = list(lower = lower, upper = upper),
    final_candidates = final_tab,
    explore = data.frame(
      source = all_source,
      matrix(all_theta, ncol = d,
             dimnames = list(NULL, paste0("theta", seq_len(d)))),
      global_gcv = all_q,
      stringsAsFactors = FALSE
    ),
    cgcv = list(ok = isTRUE(cgcv$ok), theta = cgcv$theta, lambda = cgcv$lambda),
    diagnostics = list(
      version = "v1",
      adaptive_box = isTRUE(adaptive_box),
      n_expansions = n_exp,
      sobol_batches_used = batch_id,
      region_stable_streak = region_stable_streak,
      seed = seed,
      M_search = M_search,
      M_final = M_final,
      note = "Experimental adaptive global-λ; not a replacement for cGCV."
    ),
    elapsed = elapsed
  )
}
