# Experimental joint global λ optimizer (box-Sobol + local refine + SAA TT-gGCV).
#
# NOT exported. NOT a product replacement for cGCV.
# Productive path remains lambda = "cGCV" (CLOSED — KEEP cGCV).
# "Global" here = extensive search inside a compact θ-box + multistart
# refinement + rigorous re-evaluation — not mathematical certification.

# ---------------------------------------------------------------------------
# Low-discrepancy points (builtin Sobol, Joe–Kuo direction numbers, d ≤ 8)
# ---------------------------------------------------------------------------

#' Joe–Kuo direction-number tables for Sobol dims 2..8 (dim 1 = van der Corput).
#' @keywords internal
#' @noRd
.tt_lab_sobol_params <- function() {
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
.tt_lab_sobol_unit <- function(n, d, skip = 1L) {
  n <- as.integer(n)
  d <- as.integer(d)
  skip <- as.integer(skip)
  stopifnot(n >= 1L, d >= 1L, d <= 8L, skip >= 0L)
  maxbit <- 30L
  params <- .tt_lab_sobol_params()
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
  # Keep off exact 0/1 for log-λ safety if ever mapped without padding
  out <- pmin(pmax(out, .Machine$double.eps), 1 - .Machine$double.eps)
  out
}

#' Map unit cube Sobol to a hyper-rectangle.
#' @keywords internal
#' @noRd
.tt_lab_sobol_box <- function(n, lower, upper, skip = 1L) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  d <- length(lower)
  stopifnot(length(upper) == d, all(upper > lower))
  U <- .tt_lab_sobol_unit(n, d, skip = skip)
  sweep(U, 2L, upper - lower, `*`) +
    matrix(lower, nrow = n, ncol = d, byrow = TRUE)
}

#' Default Sobol budget by dimension.
#' @keywords internal
#' @noRd
.tt_lab_default_n_global <- function(d) {
  d <- as.integer(d)
  if (d <= 2L) 256L else if (d == 3L) 512L else if (d == 4L) 1024L else 256L * d
}

# ---------------------------------------------------------------------------
# Diversity selection among promising θ points
# ---------------------------------------------------------------------------

#' Greedy diversity keep: best Q, then next if far enough in θ-space.
#' @keywords internal
#' @noRd
.tt_lab_diverse_top <- function(theta, q, n_keep = 10L, min_dist = 0.75) {
  theta <- as.matrix(theta)
  q <- as.numeric(q)
  n_keep <- as.integer(n_keep)
  stopifnot(nrow(theta) == length(q), n_keep >= 1L)
  ok <- is.finite(q)
  if (!any(ok)) return(integer(0))
  ord <- order(q[ok])
  idx_all <- which(ok)[ord]
  selected <- integer(0)
  for (i in idx_all) {
    if (!length(selected)) {
      selected <- i
    } else {
      dists <- sqrt(rowSums((theta[selected, , drop = FALSE] -
                               matrix(theta[i, ], nrow = length(selected),
                                      ncol = ncol(theta), byrow = TRUE))^2))
      if (min(dists) >= min_dist) selected <- c(selected, i)
    }
    if (length(selected) >= n_keep) break
  }
  selected
}

# ---------------------------------------------------------------------------
# Anchors: cGCV, box centre, corners, isotropic line
# ---------------------------------------------------------------------------

#' Build a set of θ anchors (includes optional cGCV solution).
#' @keywords internal
#' @noRd
.tt_lab_theta_anchors <- function(d, lower, upper, theta_cgcv = NULL) {
  d <- as.integer(d)
  lower <- rep(as.numeric(lower), length.out = d)
  upper <- rep(as.numeric(upper), length.out = d)
  mid <- 0.5 * (lower + upper)
  rows <- list(centre = mid)
  # Representative corners (all-low, all-high, and one mixed pair if d>=2)
  rows$corner_lo <- lower
  rows$corner_hi <- upper
  if (d >= 2L) {
    mix <- mid
    mix[1L] <- lower[1L]
    mix[2L] <- upper[2L]
    rows$corner_mix <- mix
  }
  # Isotropic line: a few shared θ values
  iso_vals <- unique(c(mid[1L], seq(lower[1L], upper[1L], length.out = 5L)))
  for (i in seq_along(iso_vals)) {
    rows[[paste0("iso_", i)]] <- rep(iso_vals[i], d)
  }
  if (!is.null(theta_cgcv) && length(theta_cgcv) == d && all(is.finite(theta_cgcv))) {
    rows$cgcv <- as.numeric(theta_cgcv)
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Evaluator factory (deterministic SAA TT-gGCV)
# ---------------------------------------------------------------------------

#' Create a cached θ → TT-gGCV evaluator with fixed init cores + MC bank.
#' @keywords internal
#' @noRd
.tt_lab_make_theta_evaluator <- function(y, X, rank, common_init, probes,
                                         control, k, degree, penalty_order,
                                         epsilon_rel = 1e-3,
                                         scheme = "forward",
                                         lower, upper,
                                         cache = new.env(parent = emptyenv()),
                                         gdf_warm_start = TRUE,
                                         fit_backend = c("R", "Rcpp_fixed")) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  d <- ncol(X)
  lower <- rep(as.numeric(lower), length.out = d)
  upper <- rep(as.numeric(upper), length.out = d)
  fit_backend <- .tt_lab_match_fit_backend(fit_backend)
  force(list(y, X, rank, common_init, probes, control, k, degree,
             penalty_order, epsilon_rel, scheme, lower, upper, cache,
             gdf_warm_start, fit_backend))

  # Cost counters (shared across closures)
  cost <- new.env(parent = emptyenv())
  cost$n_theta_calls <- 0L          # including cache hits + barrier
  cost$n_theta_miss <- 0L           # full TT-gGCV evaluations
  cost$n_theta_hit <- 0L
  cost$n_als_fits_approx <- 0L      # base + M probes (forward) per miss
  cost$n_barrier <- 0L

  eval_one <- function(theta, M, mc_bank_id = "bank0",
                       init_cores = common_init,
                       control_use = control,
                       return_fit = TRUE,
                       use_cache = TRUE) {
    theta <- as.numeric(theta)
    cost$n_theta_calls <- cost$n_theta_calls + 1L
    if (length(theta) != d) stop("theta length must equal d.", call. = FALSE)
    # Soft barrier outside the box
    if (any(!is.finite(theta)) || any(theta < lower - 1e-12) ||
        any(theta > upper + 1e-12)) {
      cost$n_barrier <- cost$n_barrier + 1L
      return(list(
        global_gcv = Inf, gdf = NA_real_, gdf_mc_se = NA_real_,
        rss = NA_real_, lambda = 10^theta, theta = theta,
        valid = FALSE, boundary_violation = TRUE, fit = NULL,
        cache_hit = FALSE, M = as.integer(M)
      ))
    }
    theta_c <- pmin(pmax(theta, lower), upper)
    lambda <- 10^theta_c
    # High precision key: %.6f collapses nlminb finite-difference steps
    key <- paste(
      paste(sprintf("%.12f", theta_c), collapse = ","),
      as.integer(M), mc_bank_id, sep = "|"
    )
    if (isTRUE(use_cache) && exists(key, envir = cache, inherits = FALSE)) {
      cost$n_theta_hit <- cost$n_theta_hit + 1L
      out <- cache[[key]]
      out$cache_hit <- TRUE
      if (!isTRUE(return_fit)) out$fit <- NULL
      return(out)
    }
    cost$n_theta_miss <- cost$n_theta_miss + 1L
    M_use <- as.integer(M)
    # forward scheme: 1 base + M perturbed fits
    cost$n_als_fits_approx <- cost$n_als_fits_approx + 1L + M_use
    ev <- .tt_with_preserved_seed({
      tt_global_gcv(
        lambda = lambda,
        y = y,
        X = X,
        rank = rank,
        probes = probes,
        M = M_use,
        epsilon_rel = epsilon_rel,
        scheme = scheme,
        control = control_use,
        init = .tt_clone_cores(init_cores),
        fit_backend = fit_backend,
        k = k,
        degree = degree,
        penalty_order = penalty_order,
        on_nonconverged = "na",
        warm_start = gdf_warm_start
      )
    })
    out <- list(
      global_gcv = as.numeric(ev$global_gcv)[1L],
      gdf = as.numeric(ev$gdf)[1L],
      gdf_mc_se = as.numeric(ev$gdf_mc_se)[1L],
      rss = as.numeric(ev$rss)[1L],
      lambda = as.numeric(ev$lambda %||% lambda),
      theta = theta_c,
      valid = isTRUE(ev$valid),
      boundary_violation = FALSE,
      fit = if (isTRUE(return_fit)) ev$fit else NULL,
      fit_ok = ev$fit_ok,
      cache_hit = FALSE,
      M = M_use,
      invalid_reasons = ev$invalid_reasons %||% character(0)
    )
    # Store without bulky fit in cache to limit memory
    if (isTRUE(use_cache)) {
      cached <- out
      cached$fit <- NULL
      cache[[key]] <- cached
    }
    out
  }

  list(eval = eval_one, cache = cache, cost = cost, d = d,
       lower = lower, upper = upper)
}

#' Fit cGCV once (anchor); returns theta = log10(lambda).
#' @keywords internal
#' @noRd
.tt_lab_cgcv_anchor <- function(y, X, rank, common_init, control,
                                k, degree, penalty_order) {
  ctrl <- control
  ctrl$cgcv_update <- "outer_simultaneous"
  ctrl$warn_lambda_boundary <- FALSE
  fit <- tryCatch(
    .tt_with_preserved_seed({
      ttps(
        y, X,
        family = stats::gaussian(),
        rank = rank,
        k = k,
        degree = degree,
        penalty_order = penalty_order,
        lambda = "cGCV",
        optimizer = "ALS",
        backend = "R",
        init = .tt_clone_cores(common_init),
        control = ctrl
      )
    }),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, theta = NULL, lambda = NULL, fit = NULL,
                reason = conditionMessage(fit)))
  }
  lam <- as.numeric(fit$lambda)
  list(
    ok = TRUE,
    theta = log10(lam),
    lambda = lam,
    fit = fit,
    reason = NA_character_
  )
}

#' ALS multistart at fixed λ; pick basin by penalized objective, then GCV.
#' @keywords internal
#' @noRd
.tt_lab_reeval_multistart <- function(theta, y, X, rank, probes, M,
                                      n_starts, seed, control,
                                      k, degree, penalty_order,
                                      epsilon_rel, scheme,
                                      lower, upper,
                                      fit_backend = c("R", "Rcpp_fixed")) {
  d <- ncol(as.matrix(X))
  theta <- pmin(pmax(as.numeric(theta), lower), upper)
  lambda <- 10^theta
  ctrl <- .tt_lab_refit_control(control)
  fit_backend <- .tt_lab_match_fit_backend(fit_backend)
  fits <- vector("list", n_starts)
  objs <- rep(Inf, n_starts)
  for (b in seq_len(n_starts)) {
    init_b <- .tt_with_preserved_seed({
      tt_initialize(d = d, rank = rank, k = k, seed = as.integer(seed) + 1000L * b)
    })
    fit_b <- tryCatch(
      .tt_with_preserved_seed({
        .tt_lab_fit_fixed(
          y = y,
          X = X,
          lambda = lambda,
          rank = rank,
          control = ctrl,
          fit_backend = fit_backend,
          init = init_b,
          k = k,
          degree = degree,
          penalty_order = penalty_order
        )
      }),
      error = function(e) NULL
    )
    if (is.null(fit_b)) next
    obj <- tryCatch(tt_objective(fit_b, X)$value, error = function(e) Inf)
    fits[[b]] <- fit_b
    objs[b] <- as.numeric(obj)[1L]
  }
  if (!any(is.finite(objs))) {
    return(list(
      ok = FALSE, global_gcv = Inf, theta = theta, lambda = lambda,
      gdf = NA_real_, gdf_mc_se = NA_real_, fit = NULL,
      basin_index = NA_integer_, penalized_obj = NA_real_
    ))
  }
  b_star <- which.min(objs)
  fit_star <- fits[[b_star]]
  # GCV from this basin (warm GDF from its cores)
  gdf_res <- tt_global_gdf_mc(
    fit_base = fit_star,
    y = y,
    M = as.integer(M),
    epsilon_rel = epsilon_rel,
    scheme = scheme,
    probes = probes[, seq_len(min(as.integer(M), ncol(probes))), drop = FALSE],
    warm_start = TRUE,
    on_nonconverged = "na",
    control = ctrl
  )
  n <- length(y)
  rss <- sum((y - as.numeric(fitted(fit_star)))^2)
  gdf <- gdf_res$gdf
  denom <- (n - gdf)^2
  gcv <- if (is.finite(rss) && is.finite(denom) && denom >= 1e-12 &&
             is.finite(gdf) && gdf >= 0 && gdf < n) {
    n * rss / denom
  } else {
    Inf
  }
  list(
    ok = is.finite(gcv),
    global_gcv = as.numeric(gcv),
    theta = theta,
    lambda = lambda,
    gdf = as.numeric(gdf),
    gdf_mc_se = as.numeric(gdf_res$gdf_mc_se),
    rss = as.numeric(rss),
    fit = fit_star,
    basin_index = as.integer(b_star),
    penalized_obj = as.numeric(objs[b_star]),
    n_starts = as.integer(n_starts)
  )
}

# ---------------------------------------------------------------------------
# Main optimizer
# ---------------------------------------------------------------------------

#' Experimental joint global λ optimizer (Sobol + refine + SAA TT-gGCV).
#'
#' Optimizes \(\boldsymbol\theta=\log_{10}\boldsymbol\lambda\) jointly inside a
#' compact box. Not exported. Does not replace operational `cGCV`.
#'
#' @param y,X Gaussian response and covariates (preferred API).
#' @param formula,data Optional formula interface (`y ~ x1 + x2`); ignored if
#'   `y`/`X` supplied.
#' @param rank Fixed TT rank.
#' @param family Currently only [stats::gaussian()].
#' @param theta_lower,theta_upper Box bounds (scalar or length-d).
#' @param n_global Sobol sample size (default by `d` if `NULL`).
#' @param n_refine Number of local refinements (`nlminb` bounded).
#' @param n_diverse How many diverse Sobol elites to consider before refine.
#' @param M_search,M_final SAA probe counts for search vs final re-eval.
#' @param core_starts_search Core inits per Sobol/refine eval (default 1).
#' @param core_starts_final Core inits for rigorous re-eval (basin by ALS
#'   objective, then GCV).
#' @param min_dist Minimum θ-distance for diversity clustering.
#' @param boundary_tol Distance to box edge that flags `boundary = TRUE`.
#' @param seed Master seed (common cores + MC banks).
#' @param k,degree,penalty_order Basis / penalty.
#' @param control [tt_control()] for ALS.
#' @param include_cgcv_anchor If `TRUE`, evaluate and seed with cGCV.
#' @param extra_theta Optional matrix of additional θ starts (`n × d`), e.g. a
#'   prior solution used as an explicit anchor (`prior_solution_*`).
#' @param epsilon_rel,scheme Forwarded to MC-GDF.
#' @param verbose Print stage progress.
#' @return Unevaluated diagnostic list (see package lab docs).
#' @keywords internal
tt_global_lambda_optimize <- function(y = NULL,
                                      X = NULL,
                                      formula = NULL,
                                      data = NULL,
                                      rank,
                                      family = stats::gaussian(),
                                      theta_lower = -6,
                                      theta_upper = 6,
                                      n_global = NULL,
                                      n_refine = 5L,
                                      n_diverse = 10L,
                                      M_search = 25L,
                                      M_final = 200L,
                                      core_starts_search = 1L,
                                      core_starts_final = 5L,
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
                                      fit_backend = c("R", "Rcpp_fixed"),
                                      verbose = FALSE) {
  t_wall0 <- proc.time()[["elapsed"]]
  scheme <- match.arg(scheme)
  fit_backend <- .tt_lab_match_fit_backend(fit_backend)
  fam_key <- family_key(normalize_family(family))
  if (!identical(fam_key, "gaussian")) {
    stop("tt_global_lambda_optimize: Gaussian only in v0.", call. = FALSE)
  }

  if (is.null(y) || is.null(X)) {
    if (is.null(formula) || is.null(data)) {
      stop("Supply y+X or formula+data.", call. = FALSE)
    }
    mf <- stats::model.frame(formula, data = data)
    y <- stats::model.response(mf)
    X <- stats::model.matrix(formula, data = data)
    # drop intercept column if present
    if ("(Intercept)" %in% colnames(X)) {
      X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
    }
  }
  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y)
  d <- ncol(X)
  if (nrow(X) != n) stop("nrow(X) must equal length(y).", call. = FALSE)

  lower <- rep(as.numeric(theta_lower), length.out = d)
  upper <- rep(as.numeric(theta_upper), length.out = d)
  if (any(upper <= lower)) stop("theta_upper must exceed theta_lower.", call. = FALSE)

  if (is.null(n_global)) n_global <- .tt_lab_default_n_global(d)
  n_global <- as.integer(n_global)
  n_refine <- as.integer(n_refine)
  n_diverse <- as.integer(n_diverse)
  M_search <- as.integer(M_search)
  M_final <- as.integer(M_final)
  core_starts_search <- as.integer(core_starts_search)
  core_starts_final <- as.integer(core_starts_final)
  seed <- as.integer(seed)
  k <- as.integer(k)

  ctrl <- .tt_lab_refit_control(control)
  ctrl$seed <- seed

  # Common cores + MC banks (deterministic SAA)
  common_init <- .tt_with_preserved_seed({
    tt_initialize(d = d, rank = rank, k = k, seed = seed)
  })
  probes_search <- .tt_lab_rademacher_probes(n, M_search, probe_seed = seed + 11L)
  # Nested prefix: final bank extends search bank when M_final >= M_search
  M_bank <- max(M_search, M_final)
  probes_bank <- .tt_lab_rademacher_probes(n, M_bank, probe_seed = seed + 11L)
  probes_search <- probes_bank[, seq_len(M_search), drop = FALSE]
  probes_final <- probes_bank[, seq_len(M_final), drop = FALSE]
  probes_alt <- .tt_lab_rademacher_probes(n, M_final, probe_seed = seed + 911L)

  evaluator <- .tt_lab_make_theta_evaluator(
    y = y, X = X, rank = rank, common_init = common_init,
    probes = probes_search, control = ctrl, k = k, degree = degree,
    penalty_order = penalty_order, epsilon_rel = epsilon_rel, scheme = scheme,
    lower = lower, upper = upper, fit_backend = fit_backend
  )

  # --- Stage 0: cGCV anchor -------------------------------------------------
  cgcv <- list(ok = FALSE)
  if (isTRUE(include_cgcv_anchor)) {
    if (verbose) message("[global-lambda] cGCV anchor...")
    cgcv <- .tt_lab_cgcv_anchor(
      y, X, rank, common_init, ctrl, k, degree, penalty_order
    )
  }

  # --- Stage A/B: Sobol exploration + anchors -------------------------------
  use_sobol <- n_global > 0L
  if (verbose) {
    message(sprintf(
      "[global-lambda] exploration n_sobol=%d d=%d M=%d cgcv_anchor=%s",
      if (use_sobol) n_global else 0L, d, M_search, include_cgcv_anchor
    ))
  }
  if (use_sobol) {
    sobol_theta <- .tt_lab_sobol_box(n_global, lower, upper, skip = 1L)
    source_sobol <- rep("sobol", nrow(sobol_theta))
  } else {
    sobol_theta <- matrix(numeric(0), nrow = 0L, ncol = d)
    source_sobol <- character(0)
  }
  anchors <- .tt_lab_theta_anchors(
    d, lower, upper,
    theta_cgcv = if (isTRUE(cgcv$ok)) cgcv$theta else NULL
  )
  # Rename cgcv row for provenance clarity
  an_names <- rownames(anchors)
  an_names[an_names == "cgcv"] <- "cgcv_anchor"
  an_names[an_names == "centre"] <- "center"
  an_names[grepl("^iso_", an_names)] <- paste0(
    "isotropic_anchor_", gsub("^iso_", "", an_names[grepl("^iso_", an_names)])
  )
  an_names[an_names == "corner_lo"] <- "corner_lo"
  an_names[an_names == "corner_hi"] <- "corner_hi"
  an_names[an_names == "corner_mix"] <- "corner_mix"
  rownames(anchors) <- an_names

  if (!is.null(extra_theta)) {
    et <- as.matrix(extra_theta)
    if (ncol(et) != d && nrow(et) == d && ncol(et) == 1L) et <- t(et)
    if (ncol(et) != d) stop("`extra_theta` must have d columns.", call. = FALSE)
    rownames(et) <- paste0("prior_solution_", seq_len(nrow(et)))
    anchors <- rbind(anchors, et)
  }

  theta_all <- rbind(sobol_theta, anchors)
  source_lab <- c(source_sobol, rownames(anchors))

  if (!nrow(theta_all)) {
    stop("No exploration points: need n_global > 0 or anchors.", call. = FALSE)
  }

  sobol_rows <- vector("list", nrow(theta_all))
  for (i in seq_len(nrow(theta_all))) {
    # Search stage: single common init (multistart reserved for final)
    ev <- evaluator$eval(theta_all[i, ], M = M_search, mc_bank_id = "bank0")
    sobol_rows[[i]] <- data.frame(
      i = i,
      source = source_lab[i],
      matrix(theta_all[i, ], nrow = 1L,
             dimnames = list(NULL, paste0("theta", seq_len(d)))),
      global_gcv = ev$global_gcv,
      gdf = ev$gdf,
      rss = ev$rss,
      valid = ev$valid,
      cache_hit = ev$cache_hit,
      stringsAsFactors = FALSE
    )
  }
  sobol_results <- do.call(rbind, sobol_rows)
  q_vec <- sobol_results$global_gcv
  th_mat <- as.matrix(sobol_results[, paste0("theta", seq_len(d)), drop = FALSE])

  elite_idx <- .tt_lab_diverse_top(
    th_mat, q_vec, n_keep = max(n_diverse, 1L), min_dist = min_dist
  )
  if (!length(elite_idx)) {
    elite_idx <- which.min(q_vec)
  }
  n_ref_use <- max(0L, as.integer(n_refine))
  refine_starts <- if (n_ref_use > 0L) {
    elite_idx[seq_len(min(n_ref_use, length(elite_idx)))]
  } else {
    integer(0)
  }

  # --- Stage C: bounded local refinement (nlminb) ---------------------------
  if (verbose) message(sprintf("[global-lambda] refining %d starts...",
                               length(refine_starts)))
  refined_rows <- vector("list", length(refine_starts))
  for (j in seq_along(refine_starts)) {
    i0 <- refine_starts[j]
    th0 <- as.numeric(th_mat[i0, ])
    obj <- function(th) {
      # use_cache=FALSE so finite-difference steps are not collapsed
      ev <- evaluator$eval(th, M = M_search, mc_bank_id = "bank0",
                           return_fit = FALSE, use_cache = FALSE)
      q <- ev$global_gcv
      if (!is.finite(q)) 1e30 else q
    }
    opt <- tryCatch(
      stats::nlminb(
        start = th0,
        objective = obj,
        lower = lower,
        upper = upper,
        control = list(abs.tol = 1e-8, rel.tol = 1e-8, iter.max = 80L)
      ),
      error = function(e) list(par = th0, objective = obj(th0),
                               convergence = 99L, message = conditionMessage(e))
    )
    th_ref <- as.numeric(opt$par)
    ev_ref <- evaluator$eval(th_ref, M = M_search, mc_bank_id = "bank0")
    refined_rows[[j]] <- data.frame(
      start_i = i0,
      start_source = source_lab[i0],
      winner_source = paste0("refined_from_", source_lab[i0]),
      start_gcv = q_vec[i0],
      matrix(th_ref, nrow = 1L,
             dimnames = list(NULL, paste0("theta", seq_len(d)))),
      global_gcv = ev_ref$global_gcv,
      gdf = ev_ref$gdf,
      rss = ev_ref$rss,
      nlminb_conv = as.integer(opt$convergence %||% NA_integer_),
      improved = is.finite(ev_ref$global_gcv) && is.finite(q_vec[i0]) &&
        ev_ref$global_gcv < q_vec[i0] - 1e-12,
      stringsAsFactors = FALSE
    )
  }
  refined_results <- if (length(refined_rows)) {
    do.call(rbind, refined_rows)
  } else {
    data.frame(
      start_i = integer(0), start_source = character(0),
      winner_source = character(0), start_gcv = numeric(0),
      global_gcv = numeric(0), gdf = numeric(0), rss = numeric(0),
      nlminb_conv = integer(0), improved = logical(0),
      stringsAsFactors = FALSE
    )
  }

  # --- Stage D: rigorous re-evaluation --------------------------------------
  cand_theta <- list()
  cand_source <- character(0)
  # top refined
  if (nrow(refined_results)) {
    ord_r <- order(refined_results$global_gcv)
    take_r <- head(ord_r, 5L)
    for (ii in take_r) {
      cand_theta[[length(cand_theta) + 1L]] <- as.numeric(
        refined_results[ii, paste0("theta", seq_len(d))]
      )
      cand_source <- c(cand_source, refined_results$winner_source[ii])
    }
  }
  # best exploration point without refine
  i_best_raw <- which.min(q_vec)
  cand_theta[[length(cand_theta) + 1L]] <- as.numeric(th_mat[i_best_raw, ])
  cand_source <- c(cand_source, as.character(source_lab[i_best_raw]))
  # cGCV raw (even if already among anchors)
  if (isTRUE(cgcv$ok)) {
    cand_theta[[length(cand_theta) + 1L]] <- as.numeric(cgcv$theta)
    cand_source <- c(cand_source, "cgcv_anchor")
  }
  # unique by rounded θ (keep first provenance)
  key_th <- vapply(cand_theta, function(th) {
    paste(sprintf("%.5f", th), collapse = ",")
  }, character(1))
  keep <- !duplicated(key_th)
  cand_theta <- cand_theta[keep]
  cand_source <- cand_source[keep]

  # Snapshot search+refine cost before final multistart re-eval
  cost_search <- list(
    n_theta_calls = as.integer(evaluator$cost$n_theta_calls),
    n_theta_miss = as.integer(evaluator$cost$n_theta_miss),
    n_theta_hit = as.integer(evaluator$cost$n_theta_hit),
    n_als_fits_approx = as.integer(evaluator$cost$n_als_fits_approx),
    n_barrier = as.integer(evaluator$cost$n_barrier),
    n_explore_points = as.integer(nrow(theta_all)),
    n_refine_starts = as.integer(length(refine_starts))
  )

  if (verbose) {
    message(sprintf("[global-lambda] final re-eval %d candidates (M=%d, starts=%d)...",
                    length(cand_theta), M_final, core_starts_final))
  }

  ctrl_strict <- ctrl
  ctrl_strict$tol <- min(as.numeric(ctrl$tol %||% 1e-8), 1e-10)
  ctrl_strict$max_sweeps <- max(as.integer(ctrl$max_sweeps %||% 40L), 60L)

  n_starts_alt <- max(1L, min(3L, core_starts_final))
  # Per candidate: n_starts ALS + (1 base + M) GDF on winner; ×2 banks
  n_final_als_approx <- as.integer(length(cand_theta) * (
    (core_starts_final + 1L + M_final) + (n_starts_alt + 1L + M_final)
  ))
  n_final_theta_equiv <- as.integer(length(cand_theta)) # candidates re-scored

  final_rows <- vector("list", length(cand_theta))
  final_fits <- vector("list", length(cand_theta))
  for (j in seq_along(cand_theta)) {
    th <- cand_theta[[j]]
    # Primary bank
    re1 <- .tt_lab_reeval_multistart(
      theta = th, y = y, X = X, rank = rank,
      probes = probes_final, M = M_final,
      n_starts = core_starts_final, seed = seed,
      control = ctrl_strict, k = k, degree = degree,
      penalty_order = penalty_order, epsilon_rel = epsilon_rel,
      scheme = scheme, lower = lower, upper = upper,
      fit_backend = fit_backend
    )
    # Independent bank (ranking stability check)
    re2 <- .tt_lab_reeval_multistart(
      theta = th, y = y, X = X, rank = rank,
      probes = probes_alt, M = M_final,
      n_starts = n_starts_alt, seed = seed + 17L,
      control = ctrl_strict, k = k, degree = degree,
      penalty_order = penalty_order, epsilon_rel = epsilon_rel,
      scheme = scheme, lower = lower, upper = upper,
      fit_backend = fit_backend
    )
    final_fits[[j]] <- re1$fit
    final_rows[[j]] <- data.frame(
      cand = j,
      source = cand_source[j],
      matrix(th, nrow = 1L,
             dimnames = list(NULL, paste0("theta", seq_len(d)))),
      gcv_bank0 = re1$global_gcv,
      gcv_bank1 = re2$global_gcv,
      gdf_bank0 = re1$gdf,
      gdf_se_bank0 = re1$gdf_mc_se,
      penalized_obj = re1$penalized_obj,
      basin_index = re1$basin_index,
      ok = re1$ok,
      stringsAsFactors = FALSE
    )
  }
  final_tab <- do.call(rbind, final_rows)

  # Select by primary bank GCV (basin already chosen by ALS objective)
  j_best <- which.min(final_tab$gcv_bank0)
  j_alt <- which.min(final_tab$gcv_bank1)
  best_fit <- final_fits[[j_best]]
  theta_best <- as.numeric(final_tab[j_best, paste0("theta", seq_len(d))])
  theta_alt <- as.numeric(final_tab[j_alt, paste0("theta", seq_len(d))])
  gcv_best <- as.numeric(final_tab$gcv_bank0[j_best])
  gcv_alt_winner_bank0 <- as.numeric(final_tab$gcv_bank0[j_alt])
  gdf_best <- as.numeric(final_tab$gdf_bank0[j_best])
  gdf_se <- as.numeric(final_tab$gdf_se_bank0[j_best])
  winner_source <- as.character(final_tab$source[j_best])

  near_lo <- abs(theta_best - lower) <= boundary_tol
  near_hi <- abs(theta_best - upper) <= boundary_tol
  boundary <- any(near_lo | near_hi)

  # Did refinement beat best Sobol?
  best_sobol_gcv <- min(q_vec, na.rm = TRUE)
  best_refined_gcv <- if (nrow(refined_results)) {
    min(refined_results$global_gcv, na.rm = TRUE)
  } else {
    Inf
  }

  cost_total <- list(
    n_explore_points = cost_search$n_explore_points,
    n_theta_miss_search = cost_search$n_theta_miss,
    n_theta_calls_search = cost_search$n_theta_calls,
    n_als_fits_search = cost_search$n_als_fits_approx,
    n_final_candidates = n_final_theta_equiv,
    n_als_fits_final_approx = n_final_als_approx,
    n_als_fits_total_approx = cost_search$n_als_fits_approx + n_final_als_approx,
    note = paste(
      "Search ALS ≈ (1+M_search) per cache-miss TT-gGCV;",
      "final ≈ per candidate: starts + (1+M_final) GDF, ×2 banks."
    )
  )

  elapsed <- proc.time()[["elapsed"]] - t_wall0
  list(
    lambda = 10^theta_best,
    theta = theta_best,
    gcv = gcv_best,
    gdf = gdf_best,
    gdf_mc_se = gdf_se,
    fit = best_fit,
    winner_source = winner_source,
    sobol_results = sobol_results,
    refined_results = refined_results,
    final_candidates = final_tab,
    elite_index = elite_idx,
    cgcv = list(
      ok = isTRUE(cgcv$ok),
      lambda = cgcv$lambda,
      theta = cgcv$theta
    ),
    boundary = boundary,
    boundary_margins = list(near_lower = near_lo, near_upper = near_hi),
    box = list(lower = lower, upper = upper),
    cost = cost_total,
    diagnostics = list(
      best_sobol_gcv = best_sobol_gcv,
      best_refined_gcv = best_refined_gcv,
      refine_improved = is.finite(best_refined_gcv) &&
        best_refined_gcv < best_sobol_gcv - 1e-12,
      ranking_stable_alt_bank = identical(j_best, j_alt),
      alt_winner_theta = theta_alt,
      alt_winner_gcv_bank0 = gcv_alt_winner_bank0,
      n_global = n_global,
      M_search = M_search,
      M_final = M_final,
      core_starts_final = core_starts_final,
      seed = seed,
      include_cgcv_anchor = isTRUE(include_cgcv_anchor),
      winner_source = winner_source,
      criterion = "SAA_TT_gGCV",
      note = paste(
        "Experimental joint λ optimizer; not classical GCV;",
        "not a replacement for operational cGCV."
      )
    ),
    convergence = list(
      refine_nlminb = if (nrow(refined_results)) refined_results$nlminb_conv else integer(0),
      final_ok = isTRUE(final_tab$ok[j_best])
    ),
    elapsed = elapsed
  )
}
