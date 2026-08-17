# Experimental laboratory: global TT-gGCV via Monte Carlo GDF.
#
# NOT exported. NOT a validated public method. NOT classical exact GCV.
# Terminology (provisional):
#   - "global TT-gGCV"
#   - "Monte Carlo generalized degrees of freedom" (GDF)
#   - "global GDF of the converged fixed-rank TT fitter"
#
# GDF_TT != sum_k ed_k (conditional core EDFs). Rank stays fixed while
# evaluating the criterion; changing rank changes the map y |-> yhat.

#' Deep-copy a list of TT core arrays.
#' @keywords internal
#' @noRd
.tt_clone_cores <- function(cores) {
  lapply(cores, function(a) {
    out <- array(as.numeric(a), dim = dim(a))
    storage.mode(out) <- "double"
    out
  })
}

#' Elementwise equality of TT core lists (numeric, dims).
#' @keywords internal
#' @noRd
.tt_lab_cores_equal <- function(a, b, tol = 0) {
  if (length(a) != length(b)) return(FALSE)
  for (j in seq_along(a)) {
    if (!identical(dim(a[[j]]), dim(b[[j]]))) return(FALSE)
    if (isTRUE(tol > 0)) {
      if (max(abs(as.numeric(a[[j]]) - as.numeric(b[[j]]))) > tol) return(FALSE)
    } else if (!isTRUE(all.equal(as.numeric(a[[j]]), as.numeric(b[[j]]),
                                 tolerance = 0, check.attributes = FALSE))) {
      return(FALSE)
    }
  }
  TRUE
}

#' Resolve GDF perturbation init policy (P4B).
#'
#' `"probe_warm"`: each MC probe starts from a **clone** of `fit_base$cores`
#' (never chains probe→probe). `"cold"`: fresh random / NULL init per probe.
#' θ→θ `"continuation"` is **not** offered here (experimental elsewhere).
#'
#' @param gdf_init Preferred name; overrides `warm_start` when set.
#' @param warm_start Legacy: `TRUE`→`probe_warm`, `FALSE`→`cold`.
#' @keywords internal
#' @noRd
.tt_lab_match_gdf_init <- function(gdf_init = NULL, warm_start = NULL) {
  if (!is.null(gdf_init) && !identical(gdf_init, "")) {
    gdf_init <- as.character(gdf_init)[[1L]]
    return(match.arg(gdf_init, c("probe_warm", "cold")))
  }
  if (isFALSE(warm_start)) return("cold")
  "probe_warm"
}

#' Default stage-wise ALS / tol budgets for adaptive TT-gGCV fidelity (P4B).
#'
#' Sobol: cheap ranking · refine: locate valley · final: decision-quality.
#' Winner selection must use **final** fidelity only (same M bank + sweeps).
#'
#' @keywords internal
#' @noRd
.tt_lab_default_fidelity <- function() {
  list(
    sobol = list(max_sweeps = 12L, tol = 1e-8),
    refine = list(max_sweeps = 25L, tol = 1e-8),
    final = list(max_sweeps = 50L, tol = 1e-10)
  )
}

#' Apply stage fidelity overrides onto a lab control list.
#' @keywords internal
#' @noRd
.tt_lab_control_for_stage <- function(base_control,
                                      stage = c("sobol", "refine", "final"),
                                      fidelity = NULL,
                                      adaptive = TRUE) {
  stage <- match.arg(stage)
  ctrl <- .tt_lab_refit_control(base_control)
  if (!isTRUE(adaptive)) return(ctrl)
  fid <- fidelity %||% .tt_lab_default_fidelity()
  st <- fid[[stage]]
  if (is.null(st)) return(ctrl)
  if (!is.null(st$max_sweeps)) ctrl$max_sweeps <- as.integer(st$max_sweeps)
  if (!is.null(st$tol)) ctrl$tol <- as.numeric(st$tol)
  ctrl
}

#' Save / restore .Random.seed around an expression (package has no withr).
#' @keywords internal
#' @noRd
.tt_with_preserved_seed <- function(expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

#' Rademacher probe matrix (n x M) with seed isolation.
#' @keywords internal
#' @noRd
.tt_lab_rademacher_probes <- function(n, M, probe_seed = 1L) {
  n <- as.integer(n)
  M <- as.integer(M)
  stopifnot(n >= 1L, M >= 1L)
  .tt_with_preserved_seed({
    set.seed(as.integer(probe_seed))
    mat <- matrix(
      sample(c(-1, 1), size = n * M, replace = TRUE),
      nrow = n, ncol = M
    )
  })
  storage.mode(mat) <- "double"
  mat
}

#' Finite-difference step from relative scale of y.
#'
#' \eqn{\epsilon = \texttt{epsilon_rel} \|y\|_2 / \sqrt{n}}.
#' If \|y\|_2 is near zero, fall back to \texttt{epsilon_rel} (absolute).
#' @keywords internal
#' @noRd
.tt_lab_epsilon <- function(y, epsilon_rel = 1e-3) {
  y <- as.numeric(y)
  n <- length(y)
  epsilon_rel <- as.numeric(epsilon_rel)[[1L]]
  if (!is.finite(epsilon_rel) || epsilon_rel <= 0) {
    stop("`epsilon_rel` must be finite and > 0.", call. = FALSE)
  }
  scale <- sqrt(sum(y * y) / max(n, 1L))
  near_zero <- !is.finite(scale) || scale < .Machine$double.eps * 100
  eps <- if (near_zero) {
    epsilon_rel
  } else {
    epsilon_rel * scale
  }
  list(
    epsilon = as.numeric(eps),
    y_scale = as.numeric(scale),
    near_zero_y = isTRUE(near_zero),
    epsilon_rel = epsilon_rel
  )
}

#' Convergence diagnostics for lab refits.
#'
#' Fixed-λ Gaussian ALS in this package currently sets `converged = TRUE`
#' whenever the sweep loop finishes (even at `max_sweeps`). For the lab we:
#' - treat an explicit `converged = FALSE` as failure;
#' - record whether the sweep cap was hit and the last relative RSS change;
#' - optionally require a soft RSS tolerance (`require_rss_tol = TRUE`).
#'
#' Basin jumps / non-descent are reported in diagnostics, not hidden.
#' @keywords internal
#' @noRd
.tt_lab_fit_ok <- function(fit, control = NULL, require_rss_tol = FALSE,
                           soft_tol = 1e-4) {
  if (is.null(fit)) {
    return(list(
      ok = FALSE, reason = "NULL fit", n_sweeps = NA_integer_,
      hit_max_sweeps = NA, last_rel_rss = NA_real_
    ))
  }
  n_sweeps <- as.integer(fit$n_sweeps %||% NA_integer_)
  max_sw <- as.integer(
    (fit$control$max_sweeps %||% control$max_sweeps %||% NA_integer_)
  )
  hit_max <- is.finite(max_sw) && is.finite(n_sweeps) && n_sweeps >= max_sw
  last_rel <- NA_real_
  hist <- fit$history
  rss_hist <- NULL
  if (is.data.frame(hist) && nrow(hist) >= 2L && "rss" %in% names(hist)) {
    rss_hist <- as.numeric(hist$rss)
  } else if (is.list(hist) && length(hist) >= 2L) {
    rss_hist <- vapply(hist, function(h) as.numeric(h$rss %||% NA_real_), numeric(1))
  }
  if (!is.null(rss_hist) && length(rss_hist) >= 2L) {
    i1 <- length(rss_hist)
    i0 <- i1 - 1L
    if (all(is.finite(rss_hist[c(i0, i1)]))) {
      last_rel <- abs(rss_hist[i0] - rss_hist[i1]) / max(1, abs(rss_hist[i0]))
    }
  }

  if (!isTRUE(fit$converged)) {
    return(list(
      ok = FALSE, reason = "package_converged_FALSE",
      n_sweeps = n_sweeps, hit_max_sweeps = hit_max, last_rel_rss = last_rel
    ))
  }

  if (isTRUE(require_rss_tol)) {
    tol_use <- as.numeric(soft_tol)
    if (is.finite(last_rel) && last_rel >= tol_use) {
      return(list(
        ok = FALSE, reason = "rss_tol_not_met",
        n_sweeps = n_sweeps, hit_max_sweeps = hit_max, last_rel_rss = last_rel
      ))
    }
  }

  reason <- if (isTRUE(hit_max)) {
    "package_converged_hit_max_sweeps"
  } else if (is.finite(last_rel) && last_rel < soft_tol) {
    "package_converged_rss_stable"
  } else {
    "package_converged"
  }
  list(
    ok = TRUE, reason = reason, n_sweeps = n_sweeps,
    hit_max_sweeps = hit_max, last_rel_rss = last_rel
  )
}

#' Clone control with optional stricter ALS settings for perturbed refits.
#' @keywords internal
#' @noRd
.tt_lab_refit_control <- function(control, max_sweeps = NULL, tol = NULL) {
  if (!inherits(control, "tt_control")) {
    control <- do.call(tt_control, as.list(control))
  }
  ctrl <- control
  if (!is.null(max_sweeps)) ctrl$max_sweeps <- as.integer(max_sweeps)
  if (!is.null(tol)) ctrl$tol <- as.numeric(tol)
  ctrl$compute_edf <- FALSE
  ctrl$trace <- FALSE
  ctrl$monitor <- FALSE
  ctrl$warn_lambda_boundary <- FALSE
  ctrl
}

#' Resolve lab fixed-λ fit backend.
#'
#' @param fit_backend `"R"` (public [ttps()] ALS) or `"Rcpp_fixed"` (P3
#'   [tt_als_fit_fixed_global()] C++ path).
#' @keywords internal
#' @noRd
.tt_lab_match_fit_backend <- function(fit_backend = c("R", "Rcpp_fixed")) {
  fit_backend <- as.character(fit_backend)[[1L]]
  match.arg(fit_backend, c("R", "Rcpp_fixed"))
}

#' Convert P3 history data.frame to list-of-lists (ttps ALS style).
#' @keywords internal
#' @noRd
.tt_lab_history_list <- function(hist) {
  if (is.null(hist)) return(list())
  if (is.list(hist) && !is.data.frame(hist)) return(hist)
  if (!is.data.frame(hist) || nrow(hist) < 1L) return(list())
  lapply(seq_len(nrow(hist)), function(i) {
    list(
      sweep = as.integer(hist$sweep[i]),
      rss = as.numeric(hist$rss[i]),
      objective = as.numeric(hist$objective[i]),
      penalty = as.numeric(hist$penalty[i]),
      d_eta = as.numeric(hist$d_eta[i] %||% NA_real_),
      lambda = NULL
    )
  })
}

#' Wrap P3 fixed-λ ALS output as a minimal `"ttpspline"` for the GCV lab.
#' @keywords internal
#' @noRd
.tt_lab_wrap_fixed_fit <- function(raw, y, X, basis, knots, cyclic, ranks,
                                   lambda, control, k, degree, penalty_order,
                                   offset, weights, fit_backend,
                                   init_used = NULL) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y)
  d <- ncol(X)
  eta <- as.numeric(raw$eta)
  rss <- as.numeric(raw$rss)
  hist <- .tt_lab_history_list(raw$history)
  penalties <- tryCatch(
    tt_core_penalties_from_basis(ranks, basis, penalty_order),
    error = function(e) NULL
  )
  structure(
    list(
      call = match.call(),
      family = stats::gaussian(),
      family_key = "gaussian",
      y = y,
      X = X,
      d = d,
      n = n,
      k = as.integer(ncol(basis[[1L]])),
      degree = as.integer(degree),
      knots = knots,
      cyclic = cyclic,
      penalty_order = as.integer(penalty_order),
      penalty_mode = "global",
      cores = raw$cores,
      penalties = penalties,
      rank = ranks,
      rank_internal = ranks[-c(1L, length(ranks))],
      rank_max = max(ranks),
      lambda = as.numeric(lambda),
      lambda_method = "fixed",
      lambda_bounds = control$lambda_bounds %||% c(1e-4, 1e4),
      lambda_boundary = rep("interior", d),
      lambda_at_boundary = FALSE,
      intercept = as.numeric(raw$intercept)[1L],
      beta = numeric(0),
      linear = NULL,
      smooth = NULL,
      offset = normalize_offset(offset, n),
      null_space = "joint",
      null_space_info = NULL,
      weights = normalize_weights(weights, n),
      fitted.values = eta,
      linear.predictors = eta,
      residuals = y - eta,
      deviance = rss,
      edf = NA_real_,
      edf_tt = NA_real_,
      edf_margin = NULL,
      edf_margin_cond = NULL,
      edf_note = "lab fixed-λ fit (EDF not computed)",
      npar_tt = NA_integer_,
      npar_tt_intrinsic = NA_integer_,
      npar_dense = NA_integer_,
      compression_ratio = NA_real_,
      inference = NULL,
      ._inf = new.env(parent = emptyenv()),
      converged = isTRUE(raw$converged) || identical(raw$convergence_reason, "tol_rss") ||
        (is.finite(raw$n_sweeps) && raw$n_sweeps >= 1L),
      convergence = list(
        overall = isTRUE(raw$converged),
        pirls = NA,
        als = isTRUE(raw$converged),
        reason = as.character(raw$convergence_reason %||% NA_character_)
      ),
      optimizer = "ALS",
      optimizer_requested = "ALS",
      optimizer_used = "ALS",
      optimizer_reason = "lab fixed-λ dispatcher",
      n_sweeps = as.integer(raw$n_sweeps),
      n_pirls = NA_integer_,
      n_opt_iter = NA_integer_,
      n_outer = NA_integer_,
      n_criterion_evals = 0L,
      history = hist,
      q_descent = list(checked = FALSE),
      cgcv = NULL,
      backend = if (identical(fit_backend, "Rcpp_fixed")) "Rcpp_fixed" else "R",
      fit_backend = fit_backend,
      sparse_backend = control$sparse %||% "auto",
      timing = NA_real_,
      control = control,
      x_names = colnames(X),
      x_range = apply(X, 2L, range)
    ),
    class = "ttpspline"
  )
}

#' Single dispatcher for lab Gaussian fixed-λ ALS fits (P4A).
#'
#' All Sobol / GDF / refine / final fixed-λ evaluations should call this.
#' Does **not** change outer GCV algorithm — only the ALS backend.
#'
#' @param fit_backend `"R"` → [ttps()] ALS; `"Rcpp_fixed"` → P3 C++ fitter.
#' @keywords internal
#' @noRd
.tt_lab_fit_fixed <- function(y,
                              X,
                              lambda,
                              rank,
                              control,
                              fit_backend = c("R", "Rcpp_fixed"),
                              init = NULL,
                              k = 8L,
                              degree = 3L,
                              penalty_order = 2L,
                              knots = NULL,
                              cyclic = NULL,
                              period = NULL,
                              offset = NULL,
                              weights = NULL,
                              ...) {
  fit_backend <- .tt_lab_match_fit_backend(fit_backend)
  y <- as.numeric(y)
  X <- as.matrix(X)
  d <- ncol(X)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (length(lambda) != d) stop("`lambda` length must be 1 or d.", call. = FALSE)
  ctrl <- .tt_lab_refit_control(control)

  if (identical(fit_backend, "R")) {
    fit <- ttps(
      y = y,
      X = X,
      family = stats::gaussian(),
      rank = rank,
      k = k,
      degree = degree,
      penalty_order = penalty_order,
      lambda = lambda,
      optimizer = "ALS",
      backend = "R",
      init = init,
      control = ctrl,
      knots = knots,
      cyclic = cyclic,
      period = period,
      offset = offset,
      weights = weights,
      ...
    )
    fit$fit_backend <- "R"
    return(fit)
  }

  # Rcpp_fixed: P3 multi-sweep fitter (global P_k^full), wrap as ttpspline
  if (!exists("tt_als_fit_fixed_global_cpp", mode = "function")) {
    stop("fit_backend='Rcpp_fixed' requires compiled P3 fitter.", call. = FALSE)
  }
  bs <- build_marginal_bases(
    X, k = k, degree = degree, knots = knots,
    cyclic = cyclic, period = period
  )
  basis <- bs$basis
  ranks <- tt_rank(rank, d = d)
  if (is.null(init)) {
    init <- initialize_tt_cores(
      ncol(basis[[1L]]), ranks,
      seed = ctrl$seed %||% 1L,
      sd = ctrl$init_sd %||% 0.1
    )
  }
  raw <- tt_als_fit_fixed_global(
    y = y,
    cores = init,
    basis = basis,
    lambda = lambda,
    offset = offset,
    weights = weights,
    penalty_order = penalty_order,
    max_sweeps = as.integer(ctrl$max_sweeps %||% 50L),
    tol = as.numeric(ctrl$tol %||% 1e-8),
    backend = "Rcpp"
  )
  .tt_lab_wrap_fixed_fit(
    raw = raw, y = y, X = X, basis = basis, knots = bs$knots,
    cyclic = bs$cyclic, ranks = ranks, lambda = lambda, control = ctrl,
    k = k, degree = degree, penalty_order = penalty_order,
    offset = offset, weights = weights, fit_backend = "Rcpp_fixed"
  )
}

#' Default Gaussian fixed-λ refit via [.tt_lab_fit_fixed()].
#'
#' For `gdf_init = "probe_warm"`, always starts from a **fresh clone** of
#' `fit_base$cores` (or `init_cores` if supplied). Never reuses a previous
#' probe's cores.
#'
#' @keywords internal
#' @noRd
.tt_lab_refit_from_base <- function(y_new,
                                    fit_base,
                                    warm_start = TRUE,
                                    gdf_init = NULL,
                                    init_cores = NULL,
                                    control = NULL,
                                    fit_backend = NULL,
                                    backend = NULL,
                                    ...) {
  stopifnot(inherits(fit_base, "ttpspline"))
  if (!identical(fit_base$family_key %||% family_key(fit_base$family), "gaussian")) {
    stop("Experimental global GDF lab currently supports Gaussian only.",
         call. = FALSE)
  }
  if (!identical(fit_base$lambda_method, "fixed")) {
    stop(
      "Base fit must use numeric fixed lambda (lambda_method = 'fixed'). ",
      "Do not use lambda = 'cGCV' inside this laboratory.",
      call. = FALSE
    )
  }
  y_new <- as.numeric(y_new)
  if (length(y_new) != fit_base$n) {
    stop("Perturbed y length must match fit_base$n.", call. = FALSE)
  }
  ctrl <- .tt_lab_refit_control(control %||% fit_base$control)
  fb <- fit_backend %||% fit_base$fit_backend %||% "R"
  # Legacy alias: backend="Rcpp" from older call sites → Rcpp_fixed
  if (is.null(fit_backend) && !is.null(backend)) {
    if (identical(backend, "Rcpp") || identical(backend, "Rcpp_fixed")) {
      fb <- "Rcpp_fixed"
    } else if (identical(backend, "R")) {
      fb <- "R"
    }
  }
  fb <- .tt_lab_match_fit_backend(fb)
  gi <- .tt_lab_match_gdf_init(gdf_init = gdf_init, warm_start = warm_start)
  init <- NULL
  if (identical(gi, "probe_warm")) {
    src <- init_cores %||% fit_base$cores
    if (is.null(src)) {
      stop("probe_warm requires fit_base$cores or init_cores.", call. = FALSE)
    }
    init <- .tt_clone_cores(src)
  }
  .tt_lab_fit_fixed(
    y = y_new,
    X = fit_base$X,
    lambda = as.numeric(fit_base$lambda),
    rank = fit_base$rank,
    control = ctrl,
    fit_backend = fb,
    init = init,
    k = fit_base$k,
    degree = fit_base$degree %||% 3L,
    penalty_order = fit_base$penalty_order %||% 2L,
    knots = fit_base$knots,
    cyclic = fit_base$cyclic,
    offset = fit_base$offset,
    weights = fit_base$weights,
    ...
  )
}

#' Monte Carlo estimate of global GDF of a fixed-rank TT fitter.
#'
#' Experimental / internal. Perturbs the full converged map
#' \eqn{\hat y_{\boldsymbol\lambda}=\mathcal A_r(y;\boldsymbol\lambda)}
#' with shared Rademacher probes. Does **not** use conditional core EDFs.
#'
#' @param fit_base A Gaussian `"ttpspline"` fit with **fixed** numeric `lambda`.
#' @param y Response used for the base fit (defaults to `fit_base$y`).
#' @param refit_fun Function `(y_new, fit_base, ...)` returning a `"ttpspline"`
#'   (default: warm-started [ttps()] ALS at the same rank / lambda).
#' @param M Number of Monte Carlo probes.
#' @param epsilon_rel Relative FD step (see `.tt_lab_epsilon()`).
#' @param scheme `"forward"` or `"central"`.
#' @param probe_seed Seed for Rademacher draws (RNG state restored).
#' @param probes Optional `n x M` probe matrix (overrides `probe_seed` / `M`).
#' @param warm_start Legacy alias (`TRUE`→`probe_warm`, `FALSE`→`cold`).
#' @param gdf_init `"probe_warm"` (clone `fit_base$cores` per probe) or `"cold"`.
#'   Probes never chain; `fit_base$cores` must remain unchanged.
#' @param on_nonconverged `"error"`, `"warn"`, or `"na"` (record NA contribution).
#' @param control Optional [tt_control()] override for perturbed refits.
#' @param fit_backend Lab ALS backend for default refits (`"R"` / `"Rcpp_fixed"`).
#' @param ... Passed to `refit_fun`.
#' @return Diagnostic list with `gdf`, contributions, SE, probes, etc.
#' @keywords internal
tt_global_gdf_mc <- function(fit_base,
                             y = NULL,
                             refit_fun = NULL,
                             M = 5L,
                             epsilon_rel = 1e-3,
                             scheme = c("forward", "central"),
                             probe_seed = 1L,
                             probes = NULL,
                             warm_start = TRUE,
                             gdf_init = NULL,
                             on_nonconverged = c("error", "warn", "na"),
                             control = NULL,
                             fit_backend = NULL,
                             ...) {
  stopifnot(inherits(fit_base, "ttpspline"))
  scheme <- match.arg(scheme)
  on_nonconverged <- match.arg(on_nonconverged)
  gdf_init <- .tt_lab_match_gdf_init(gdf_init = gdf_init, warm_start = warm_start)
  warm_start <- identical(gdf_init, "probe_warm")
  y <- as.numeric(y %||% fit_base$y)
  n <- length(y)
  if (n != fit_base$n) stop("`y` length must match fit_base$n.", call. = FALSE)

  if (!identical(fit_base$lambda_method, "fixed")) {
    stop("fit_base must use fixed numeric lambda (not cGCV).", call. = FALSE)
  }
  fb <- fit_backend %||% fit_base$fit_backend %||% "R"
  fb <- .tt_lab_match_fit_backend(fb)

  if (is.null(probes)) {
    probes <- .tt_lab_rademacher_probes(n, M = as.integer(M), probe_seed = probe_seed)
  } else {
    probes <- as.matrix(probes)
    storage.mode(probes) <- "double"
    if (nrow(probes) != n) stop("`probes` must have n rows.", call. = FALSE)
    M <- ncol(probes)
  }
  if (M < 1L) stop("Need M >= 1 probes.", call. = FALSE)

  eps_info <- .tt_lab_epsilon(y, epsilon_rel = epsilon_rel)
  eps <- eps_info$epsilon

  base_fit_ok <- .tt_lab_fit_ok(fit_base, control = control %||% fit_base$control)
  yhat0 <- as.numeric(fitted(fit_base))
  base_rss <- sum((y - yhat0)^2)
  # Snapshot for isolation: probes must not mutate fit_base$cores
  base_cores_snap <- if (!is.null(fit_base$cores)) {
    .tt_clone_cores(fit_base$cores)
  } else {
    NULL
  }

  if (is.null(refit_fun)) {
    refit_fun <- function(y_new, fit_base, ...) {
      .tt_lab_refit_from_base(
        y_new, fit_base,
        gdf_init = gdf_init,
        init_cores = base_cores_snap,
        control = control,
        fit_backend = fb,
        ...
      )
    }
  }

  assert_base_untouched <- function() {
    if (is.null(base_cores_snap) || is.null(fit_base$cores)) return(invisible(TRUE))
    if (!.tt_lab_cores_equal(fit_base$cores, base_cores_snap)) {
      stop(
        "GDF probe mutated fit_base$cores (probe isolation violated). ",
        "Use gdf_init='probe_warm' with cloned inits only.",
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  .tt_with_preserved_seed({
  contrib <- rep(NA_real_, M)
  diag_rows <- vector("list", M)
  n_refits <- 0L

  handle_bad <- function(j, msg) {
    if (identical(on_nonconverged, "error")) {
      stop("Probe ", j, ": ", msg, call. = FALSE)
    }
    if (identical(on_nonconverged, "warn")) {
      warning("Probe ", j, ": ", msg, call. = FALSE)
    }
    # "na": leave contrib[j] as NA
  }

  for (j in seq_len(M)) {
    zj <- probes[, j]
    if (identical(scheme, "forward")) {
      y_p <- y + eps * zj
      t0 <- proc.time()[["elapsed"]]
      fit_p <- tryCatch(
        refit_fun(y_p, fit_base, ...),
        error = function(e) e
      )
      n_refits <- n_refits + 1L
      elapsed_p <- proc.time()[["elapsed"]] - t0
      assert_base_untouched()
      if (inherits(fit_p, "error")) {
        handle_bad(j, paste0("refit error: ", conditionMessage(fit_p)))
        diag_rows[[j]] <- data.frame(
          j = j, scheme = scheme, ok = FALSE, contribution = NA_real_,
          n_sweeps_p = NA_integer_, n_sweeps_m = NA_integer_,
          reason = conditionMessage(fit_p), elapsed = elapsed_p,
          lambda_ok = NA, rank_ok = NA,
          stringsAsFactors = FALSE
        )
        next
      }
      ok_p <- .tt_lab_fit_ok(fit_p, control = control %||% fit_base$control)
      lam_ok <- isTRUE(all.equal(as.numeric(fit_p$lambda), as.numeric(fit_base$lambda),
                                 tolerance = 0, check.attributes = FALSE))
      rank_ok <- identical(as.integer(fit_p$rank), as.integer(fit_base$rank))
      if (!isTRUE(ok_p$ok)) {
        handle_bad(j, paste0("perturbed refit not converged (", ok_p$reason, ")"))
        diag_rows[[j]] <- data.frame(
          j = j, scheme = scheme, ok = FALSE, contribution = NA_real_,
          n_sweeps_p = ok_p$n_sweeps, n_sweeps_m = NA_integer_,
          reason = ok_p$reason, elapsed = elapsed_p,
          lambda_ok = lam_ok, rank_ok = rank_ok,
          stringsAsFactors = FALSE
        )
        next
      }
      if (!lam_ok || !rank_ok) {
        handle_bad(j, "perturbed refit changed lambda or rank")
        diag_rows[[j]] <- data.frame(
          j = j, scheme = scheme, ok = FALSE, contribution = NA_real_,
          n_sweeps_p = ok_p$n_sweeps, n_sweeps_m = NA_integer_,
          reason = "lambda_or_rank_changed", elapsed = elapsed_p,
          lambda_ok = lam_ok, rank_ok = rank_ok,
          stringsAsFactors = FALSE
        )
        next
      }
      yhat_p <- as.numeric(fitted(fit_p))
      cj <- sum(zj * (yhat_p - yhat0)) / eps
      contrib[j] <- cj
      diag_rows[[j]] <- data.frame(
        j = j, scheme = scheme, ok = TRUE, contribution = cj,
        n_sweeps_p = ok_p$n_sweeps, n_sweeps_m = NA_integer_,
        reason = ok_p$reason, elapsed = elapsed_p,
        lambda_ok = TRUE, rank_ok = TRUE,
        stringsAsFactors = FALSE
      )
    } else {
      # central
      y_p <- y + eps * zj
      y_m <- y - eps * zj
      t0 <- proc.time()[["elapsed"]]
      fit_p <- tryCatch(
        refit_fun(y_p, fit_base, ...),
        error = function(e) e
      )
      assert_base_untouched()
      fit_m <- tryCatch(
        refit_fun(y_m, fit_base, ...),
        error = function(e) e
      )
      n_refits <- n_refits + 2L
      elapsed_p <- proc.time()[["elapsed"]] - t0
      assert_base_untouched()
      if (inherits(fit_p, "error") || inherits(fit_m, "error")) {
        msg <- if (inherits(fit_p, "error")) conditionMessage(fit_p) else conditionMessage(fit_m)
        handle_bad(j, paste0("refit error: ", msg))
        diag_rows[[j]] <- data.frame(
          j = j, scheme = scheme, ok = FALSE, contribution = NA_real_,
          n_sweeps_p = NA_integer_, n_sweeps_m = NA_integer_,
          reason = msg, elapsed = elapsed_p,
          lambda_ok = NA, rank_ok = NA,
          stringsAsFactors = FALSE
        )
        next
      }
      ok_p <- .tt_lab_fit_ok(fit_p, control = control %||% fit_base$control)
      ok_m <- .tt_lab_fit_ok(fit_m, control = control %||% fit_base$control)
      lam_ok <- isTRUE(all.equal(as.numeric(fit_p$lambda), as.numeric(fit_base$lambda),
                                 tolerance = 0, check.attributes = FALSE)) &&
        isTRUE(all.equal(as.numeric(fit_m$lambda), as.numeric(fit_base$lambda),
                         tolerance = 0, check.attributes = FALSE))
      rank_ok <- identical(as.integer(fit_p$rank), as.integer(fit_base$rank)) &&
        identical(as.integer(fit_m$rank), as.integer(fit_base$rank))
      if (!isTRUE(ok_p$ok) || !isTRUE(ok_m$ok)) {
        handle_bad(j, "central refit(s) not converged")
        diag_rows[[j]] <- data.frame(
          j = j, scheme = scheme, ok = FALSE, contribution = NA_real_,
          n_sweeps_p = ok_p$n_sweeps, n_sweeps_m = ok_m$n_sweeps,
          reason = paste(ok_p$reason, ok_m$reason, sep = ";"),
          elapsed = elapsed_p, lambda_ok = lam_ok, rank_ok = rank_ok,
          stringsAsFactors = FALSE
        )
        next
      }
      if (!lam_ok || !rank_ok) {
        handle_bad(j, "perturbed refit changed lambda or rank")
        diag_rows[[j]] <- data.frame(
          j = j, scheme = scheme, ok = FALSE, contribution = NA_real_,
          n_sweeps_p = ok_p$n_sweeps, n_sweeps_m = ok_m$n_sweeps,
          reason = "lambda_or_rank_changed", elapsed = elapsed_p,
          lambda_ok = lam_ok, rank_ok = rank_ok,
          stringsAsFactors = FALSE
        )
        next
      }
      yhat_p <- as.numeric(fitted(fit_p))
      yhat_m <- as.numeric(fitted(fit_m))
      cj <- sum(zj * (yhat_p - yhat_m)) / (2 * eps)
      contrib[j] <- cj
      diag_rows[[j]] <- data.frame(
        j = j, scheme = scheme, ok = TRUE, contribution = cj,
        n_sweeps_p = ok_p$n_sweeps, n_sweeps_m = ok_m$n_sweeps,
        reason = "ok", elapsed = elapsed_p,
        lambda_ok = TRUE, rank_ok = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }

  pert_diag <- do.call(rbind, diag_rows)
  ok_mask <- is.finite(contrib)
  n_ok <- sum(ok_mask)
  gdf <- if (n_ok > 0L) mean(contrib[ok_mask]) else NA_real_
  gdf_mc_se <- if (n_ok >= 2L) {
    stats::sd(contrib[ok_mask]) / sqrt(n_ok)
  } else {
    NA_real_
  }

  list(
    gdf = as.numeric(gdf),
    gdf_contributions = as.numeric(contrib),
    gdf_mc_se = as.numeric(gdf_mc_se),
    epsilon = eps,
    epsilon_info = eps_info,
    scheme = scheme,
    M = as.integer(M),
    M_ok = as.integer(n_ok),
    probes = probes,
    probe_seed = as.integer(probe_seed),
    base_rss = as.numeric(base_rss),
    base_fitted = yhat0,
    base_fit_ok = base_fit_ok,
    warm_start = isTRUE(warm_start),
    gdf_init = gdf_init,
    base_cores_unchanged = isTRUE(
      is.null(base_cores_snap) ||
        .tt_lab_cores_equal(fit_base$cores, base_cores_snap)
    ),
    n_refits = as.integer(n_refits),
    lambda = as.numeric(fit_base$lambda),
    rank = as.integer(fit_base$rank),
    # Explicit: never use / sum conditional EDF
    used_conditional_edf = FALSE,
    conditional_edf_note = "GDF_TT is MC divergence of full fitter; not sum_k ed_k",
    perturbation_diagnostics = pert_diag,
    on_nonconverged = on_nonconverged
  )
  }) # .tt_with_preserved_seed
}

#' Evaluate experimental global TT-gGCV at a fixed numeric lambda.
#'
#' @param lambda Numeric isotropic or length-d anisotropic smoothing.
#' @param y,X Gaussian response and covariates.
#' @param rank Fixed TT rank (scalar / chain); held fixed for GDF.
#' @param fit_fun Optional `(y, X, lambda, ...)` returning `"ttpspline"`.
#' @param probes Shared probe matrix for GDF (required for comparable surfaces).
#' @param M,epsilon_rel,scheme Passed to [tt_global_gdf_mc()].
#' @param control [tt_control()] for the base fit / refits.
#' @param init Optional warm-start cores for the base fit.
#' @param fit_backend Lab ALS backend: `"R"` ([ttps()]) or `"Rcpp_fixed"`
#'   (P3 C++ fixed-λ fitter). All fixed-λ fits in this call (base + GDF
#'   perturbations) use [.tt_lab_fit_fixed()].
#' @param backend Deprecated alias; ignored when `fit_backend` is set.
#' @param k,degree,penalty_order Basis / penalty settings.
#' @param on_nonconverged Forwarded to GDF estimator.
#' @param warm_start Legacy GDF init alias.
#' @param gdf_init `"probe_warm"` or `"cold"` for MC perturbations.
#' @param fidelity Optional stage label recorded in the result (`"sobol"`,
#'   `"refine"`, `"final"`, or custom).
#' @param ... Extra args to `fit_fun` / GDF.
#' @return Rich diagnostic list (criterion + GDF + validity flags).
#' @keywords internal
tt_global_gcv <- function(lambda,
                          y,
                          X,
                          rank,
                          probes,
                          fit_fun = NULL,
                          M = NULL,
                          epsilon_rel = 1e-3,
                          scheme = c("forward", "central"),
                          control = tt_control(max_sweeps = 40, tol = 1e-10,
                                               compute_edf = FALSE, seed = 1L),
                          init = NULL,
                          fit_backend = c("R", "Rcpp_fixed"),
                          backend = NULL,
                          k = 8L,
                          degree = 3L,
                          penalty_order = 2L,
                          on_nonconverged = c("error", "warn", "na"),
                          warm_start = TRUE,
                          gdf_init = NULL,
                          fidelity = NULL,
                          ...) {
  scheme <- match.arg(scheme)
  on_nonconverged <- match.arg(on_nonconverged)
  gdf_init <- .tt_lab_match_gdf_init(gdf_init = gdf_init, warm_start = warm_start)
  warm_start <- identical(gdf_init, "probe_warm")
  if (!is.null(backend) && missing(fit_backend)) {
    fit_backend <- if (identical(backend, "Rcpp") || identical(backend, "Rcpp_fixed")) {
      "Rcpp_fixed"
    } else {
      "R"
    }
  }
  fit_backend <- .tt_lab_match_fit_backend(fit_backend)
  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y)
  d <- ncol(X)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (length(lambda) != d) {
    stop("`lambda` must be scalar or length d.", call. = FALSE)
  }
  if (any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("`lambda` must be finite and strictly positive.", call. = FALSE)
  }
  probes <- as.matrix(probes)
  if (nrow(probes) != n) stop("`probes` must have n rows.", call. = FALSE)
  if (is.null(M)) M <- ncol(probes)
  M <- as.integer(M)
  if (M > ncol(probes)) stop("M cannot exceed ncol(probes).", call. = FALSE)
  probes_use <- probes[, seq_len(M), drop = FALSE]

  ctrl <- .tt_lab_refit_control(control)

  if (is.null(fit_fun)) {
    fit_fun <- function(y, X, lambda, ...) {
      .tt_lab_fit_fixed(
        y = y, X = X, lambda = lambda, rank = rank, control = ctrl,
        fit_backend = fit_backend, init = init, k = k, degree = degree,
        penalty_order = penalty_order, ...
      )
    }
  }

  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    fit_fun(y, X, lambda, ...),
    error = function(e) e
  )
  fit_time <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error")) {
    return(list(
      global_gcv = Inf,
      rss = NA_real_,
      gdf = NA_real_,
      gdf_mc_se = NA_real_,
      n = n,
      lambda = lambda,
      rank = rank,
      valid = FALSE,
      invalid_reasons = paste0("base_fit_error: ", conditionMessage(fit)),
      fit = NULL,
      fit_ok = list(ok = FALSE, reason = "base_fit_error", n_sweeps = NA_integer_),
      gdf_result = NULL,
      fit_time = fit_time,
      scheme = scheme,
      epsilon = NA_real_,
      M = M,
      fit_backend = fit_backend,
      criterion_name = "global_TT_gGCV",
      note = "Experimental; base fit failed"
    ))
  }

  fit_ok <- .tt_lab_fit_ok(fit, control = ctrl)
  yhat <- as.numeric(fitted(fit))
  rss <- sum((y - yhat)^2)

  gdf_res <- tt_global_gdf_mc(
    fit_base = fit,
    y = y,
    M = M,
    epsilon_rel = epsilon_rel,
    scheme = scheme,
    probes = probes_use,
    gdf_init = gdf_init,
    on_nonconverged = on_nonconverged,
    control = ctrl,
    fit_backend = fit_backend,
    ...
  )
  gdf <- gdf_res$gdf
  denom <- (n - gdf)^2
  reasons <- character(0)
  if (!isTRUE(fit_ok$ok)) reasons <- c(reasons, paste0("base_not_converged:", fit_ok$reason))
  if (!is.finite(gdf)) reasons <- c(reasons, "gdf_not_finite")
  if (is.finite(gdf) && gdf < 0) reasons <- c(reasons, "gdf_negative")
  if (is.finite(gdf) && gdf >= n) reasons <- c(reasons, "gdf_ge_n")
  if (!is.finite(denom) || denom < 1e-12) reasons <- c(reasons, "denom_too_small")
  if (isTRUE(gdf_res$M_ok) && gdf_res$M_ok < M) {
    reasons <- c(reasons, sprintf("partial_probes_ok_%d_of_%d", gdf_res$M_ok, M))
  }
  if (any(!gdf_res$perturbation_diagnostics$ok, na.rm = TRUE)) {
    reasons <- c(reasons, "some_perturbations_failed")
  }

  global_gcv <- if (length(reasons) && any(reasons %in% c(
    "gdf_not_finite", "gdf_negative", "gdf_ge_n", "denom_too_small"
  ))) {
    Inf
  } else if (is.finite(rss) && is.finite(denom) && denom >= 1e-12) {
    n * rss / denom
  } else {
    Inf
  }

  list(
    global_gcv = as.numeric(global_gcv),
    rss = as.numeric(rss),
    gdf = as.numeric(gdf),
    gdf_mc_se = gdf_res$gdf_mc_se,
    n = n,
    lambda = as.numeric(fit$lambda),
    rank = as.integer(fit$rank),
    valid = length(reasons) == 0L,
    invalid_reasons = reasons,
    fit = fit,
    fit_ok = fit_ok,
    gdf_result = gdf_res,
    fit_time = fit_time,
    scheme = scheme,
    epsilon = gdf_res$epsilon,
    M = M,
    n_sweeps = as.integer(fit_ok$n_sweeps %||% fit$n_sweeps %||% NA_integer_),
    converged = isTRUE(fit_ok$ok),
    fidelity = if (is.null(fidelity)) NA_character_ else as.character(fidelity)[[1L]],
    gdf_init = gdf_init,
    warm_start = warm_start,
    fit_backend = fit_backend,
    criterion_name = "global_TT_gGCV",
    note = "Experimental; not classical exact GCV; GDF_TT != sum ed_k"
  )
}

# ---------------------------------------------------------------------------
# Full tensor-product helpers (scattered design; small d,k only)
# ---------------------------------------------------------------------------

#' Row-wise Kronecker TP design from marginal bases (n x prod k_j).
#' @keywords internal
#' @noRd
.tt_lab_tp_design <- function(basis) {
  d <- length(basis)
  n <- nrow(basis[[1]])
  if (d == 1L) return(basis[[1]])
  X <- basis[[1]]
  for (j in 2:d) {
    Bj <- basis[[j]]
    X <- t(vapply(seq_len(n), function(i) {
      kronecker(Bj[i, ], X[i, ])
    }, numeric(ncol(X) * ncol(Bj))))
  }
  X
}

#' Anisotropic difference penalty on vec(Theta) (Kronecker sum).
#' @keywords internal
#' @noRd
.tt_lab_tp_penalty <- function(k_vec, lambda, penalty_order = 2L) {
  k_vec <- as.integer(k_vec)
  d <- length(k_vec)
  lambda <- rep(as.numeric(lambda), length.out = d)
  npar <- prod(k_vec)
  P <- matrix(0, npar, npar)
  for (m in seq_len(d)) {
    Dm <- difference_penalty(k_vec[m], penalty_order)
    blocks <- vector("list", d)
    for (j in seq_len(d)) {
      blocks[[j]] <- if (j == m) Dm else diag(k_vec[j])
    }
    K <- blocks[[1]]
    for (j in 2:d) K <- kronecker(blocks[[j]], K)
    P <- P + lambda[m] * K
  }
  P
}

#' Exact full-TP Gaussian EDF / RSS / fitted (intercept unpenalized).
#'
#' Uses stable solves: EDF = tr(M^{-1} Z'Z) with M = Z'Z + P_aug.
#' @keywords internal
#' @noRd
.tt_lab_full_tp_gaussian <- function(y, basis, lambda, penalty_order = 2L) {
  y <- as.numeric(y)
  n <- length(y)
  X <- .tt_lab_tp_design(basis)
  k_vec <- vapply(basis, ncol, integer(1))
  P <- .tt_lab_tp_penalty(k_vec, lambda, penalty_order = penalty_order)
  Z <- cbind(1, X)
  P_aug <- matrix(0, nrow(P) + 1L, ncol(P) + 1L)
  P_aug[-1, -1] <- P
  ZtZ <- crossprod(Z)
  M <- ZtZ + P_aug
  coef <- tryCatch(
    solve_spd(M, as.numeric(crossprod(Z, y))),
    error = function(e) solve_spd_ridge(M, as.numeric(crossprod(Z, y)))
  )
  fitted <- as.numeric(Z %*% coef)
  rss <- sum((y - fitted)^2)
  # edf = tr(M^{-1} Z'Z)
  Minv_ZtZ <- tryCatch(
    solve_spd(M, ZtZ),
    error = function(e) solve_spd_ridge(M, ZtZ)
  )
  edf <- sum(diag(Minv_ZtZ))
  denom <- (n - edf)^2
  gcv <- if (is.finite(denom) && denom >= 1e-12) n * rss / denom else Inf
  list(
    fitted = fitted,
    rss = rss,
    edf = as.numeric(edf),
    gcv = as.numeric(gcv),
    coef = coef,
    npar = ncol(Z),
    lambda = as.numeric(lambda),
    method = "full_TP_gaussian"
  )
}

# ---------------------------------------------------------------------------
# Phase 1: objective cache, CRN evaluators, designs, flatness, optim wrapper
# ---------------------------------------------------------------------------

#' Create an empty TT-gGCV evaluation cache (environment).
#' @keywords internal
#' @noRd
.tt_lab_new_cache <- function() {
  new.env(parent = emptyenv())
}

#' Cache key for a TT-gGCV evaluation.
#' @keywords internal
#' @noRd
.tt_lab_cache_key <- function(dataset_id, rank, theta, init_policy, mc_bank_id,
                              M = NULL) {
  th <- paste(sprintf("%.8f", as.numeric(theta)), collapse = ",")
  rlab <- paste(as.integer(rank), collapse = "-")
  paste(
    as.character(dataset_id),
    rlab,
    th,
    as.character(init_policy),
    as.character(mc_bank_id),
    if (is.null(M)) "Mna" else paste0("M", as.integer(M)),
    sep = "|"
  )
}

#' Local flatness: grid points within `tol_rel` of the minimum GCV.
#' @keywords internal
#' @noRd
.tt_lab_flatness_mask <- function(gcv, tol_rel = 0.01) {
  gcv <- as.numeric(gcv)
  ok <- is.finite(gcv)
  if (!any(ok)) return(rep(FALSE, length(gcv)))
  gmin <- min(gcv[ok])
  ok & (gcv <= gmin * (1 + as.numeric(tol_rel)))
}

#' Phase-1 Gaussian design factory (d = 2).
#'
#' Scenarios chosen so full-TP GCV is computable and at least one design
#' targets an interior lambda minimum (verified empirically in scripts).
#' @keywords internal
#' @noRd
.tt_lab_phase1_make_design <- function(scenario = c("smooth_smooth",
                                                    "smooth_rough",
                                                    "strong_aniso",
                                                    "weak_signal"),
                                       n = 200L,
                                       n_test = 200L,
                                       k = 8L,
                                       seed = 1L,
                                       degree = 3L,
                                       penalty_order = 2L) {
  scenario <- match.arg(scenario)
  seed <- as.integer(seed)
  n <- as.integer(n)
  n_test <- as.integer(n_test)
  k <- as.integer(k)

  gen_X <- function(nn, s) {
    .tt_with_preserved_seed({
      set.seed(s)
      matrix(runif(nn * 2L), nn, 2L)
    })
  }
  X <- gen_X(n, seed)
  X_test <- gen_X(n_test, seed + 10007L)
  colnames(X) <- colnames(X_test) <- c("x1", "x2")

  f_fun <- switch(
    scenario,
    # Low-frequency additive: tends to push lambda away from 0 with flexible k
    smooth_smooth = function(X) {
      sin(2 * pi * X[, 1]) + cos(2 * pi * X[, 2])
    },
    # Margin 2 more oscillatory
    smooth_rough = function(X) {
      0.8 * sin(2 * pi * X[, 1]) + 1.2 * sin(6 * pi * X[, 2])
    },
    # Strong anisotropy
    strong_aniso = function(X) {
      0.4 * X[, 1] + sin(8 * pi * X[, 2]) * cos(pi * X[, 1])
    },
    # Tiny signal / large noise → flatter GCV
    weak_signal = function(X) {
      0.15 * (sin(2 * pi * X[, 1]) + sin(2 * pi * X[, 2]))
    }
  )
  sigma <- switch(
    scenario,
    smooth_smooth = 0.35,
    smooth_rough = 0.40,
    strong_aniso = 0.45,
    weak_signal = 1.00
  )

  f <- f_fun(X)
  f_test <- f_fun(X_test)
  y <- .tt_with_preserved_seed({
    set.seed(seed + 17L)
    f + stats::rnorm(n, 0, sigma)
  })
  y_test <- .tt_with_preserved_seed({
    set.seed(seed + 19L)
    f_test + stats::rnorm(n_test, 0, sigma)
  })

  list(
    scenario = scenario,
    dataset_id = sprintf("%s_n%d_k%d_s%d", scenario, n, k, seed),
    n = n,
    n_test = n_test,
    d = 2L,
    k = k,
    degree = as.integer(degree),
    penalty_order = as.integer(penalty_order),
    seed = seed,
    sigma = as.numeric(sigma),
    X = X,
    y = y,
    f = f,
    X_test = X_test,
    y_test = y_test,
    f_test = f_test,
    truth_note = paste0("scenario=", scenario, "; sigma=", sigma)
  )
}

#' Phase-1-style Gaussian design factory for d = 3.
#'
#' @keywords internal
#' @noRd
.tt_lab_phase1_make_design_d3 <- function(scenario = c("smooth_smooth",
                                                       "strong_aniso"),
                                          n = 250L,
                                          n_test = 200L,
                                          k = 6L,
                                          seed = 1L,
                                          degree = 3L,
                                          penalty_order = 2L) {
  scenario <- match.arg(scenario)
  seed <- as.integer(seed)
  n <- as.integer(n)
  n_test <- as.integer(n_test)
  k <- as.integer(k)
  d <- 3L

  gen_X <- function(nn, s) {
    .tt_with_preserved_seed({
      set.seed(s)
      matrix(runif(nn * d), nn, d)
    })
  }
  X <- gen_X(n, seed)
  X_test <- gen_X(n_test, seed + 10007L)
  colnames(X) <- colnames(X_test) <- paste0("x", seq_len(d))

  f_fun <- switch(
    scenario,
    smooth_smooth = function(X) {
      sin(2 * pi * X[, 1]) + cos(2 * pi * X[, 2]) + 0.7 * sin(2 * pi * X[, 3])
    },
    # Margin 2 rough; 1 and 3 smoother → anisotropic λ pattern
    strong_aniso = function(X) {
      0.3 * X[, 1] +
        sin(8 * pi * X[, 2]) * cos(pi * X[, 1]) +
        0.5 * cos(2 * pi * X[, 3])
    }
  )
  sigma <- switch(
    scenario,
    smooth_smooth = 0.40,
    strong_aniso = 0.50
  )

  f <- f_fun(X)
  f_test <- f_fun(X_test)
  y <- .tt_with_preserved_seed({
    set.seed(seed + 17L)
    f + stats::rnorm(n, 0, sigma)
  })
  y_test <- .tt_with_preserved_seed({
    set.seed(seed + 19L)
    f_test + stats::rnorm(n_test, 0, sigma)
  })

  list(
    scenario = scenario,
    dataset_id = sprintf("d3_%s_n%d_k%d_s%d", scenario, n, k, seed),
    n = n,
    n_test = n_test,
    d = d,
    k = k,
    degree = as.integer(degree),
    penalty_order = as.integer(penalty_order),
    seed = seed,
    sigma = as.numeric(sigma),
    X = X,
    y = y,
    f = f,
    X_test = X_test,
    y_test = y_test,
    f_test = f_test,
    truth_note = paste0("d=3; scenario=", scenario, "; sigma=", sigma)
  )
}

#' Phase-1-style Gaussian design factory for d ≥ 4.
#'
#' Extends the d=2/d=3 lab scenarios without inventing new scientific cases:
#' - `smooth_smooth`: low-frequency additive across margins
#' - `strong_aniso`: margin 2 rough (high-frequency), others smooth
#'
#' @keywords internal
#' @noRd
.tt_lab_phase1_make_design_dn <- function(d,
                                          scenario = c("smooth_smooth",
                                                       "strong_aniso"),
                                          n = 250L,
                                          n_test = 200L,
                                          k = 5L,
                                          seed = 1L,
                                          degree = 3L,
                                          penalty_order = 2L) {
  scenario <- match.arg(scenario)
  d <- as.integer(d)[1L]
  if (!is.finite(d) || d < 4L) {
    stop("`d` must be an integer >= 4.", call. = FALSE)
  }
  seed <- as.integer(seed)
  n <- as.integer(n)
  n_test <- as.integer(n_test)
  k <- as.integer(k)

  gen_X <- function(nn, s) {
    .tt_with_preserved_seed({
      set.seed(s)
      matrix(runif(nn * d), nn, d)
    })
  }
  X <- gen_X(n, seed)
  X_test <- gen_X(n_test, seed + 10007L)
  colnames(X) <- colnames(X_test) <- paste0("x", seq_len(d))

  f_fun <- switch(
    scenario,
    smooth_smooth = function(X) {
      s <- 0
      for (j in seq_len(d)) {
        amp <- 1 / sqrt(j)
        s <- s + amp * sin(2 * pi * X[, j] + 0.3 * j)
      }
      s
    },
    # Margin 2 rough; others smoother → anisotropic λ (same role as d=2/3)
    strong_aniso = function(X) {
      # Match d=3 frequency (8π) so modest k can resolve the rough margin;
      # amplitude kept ≥1 so λ2 stays the identifiable small-λ direction.
      s <- 0.25 * sin(2 * pi * X[, 1]) +
        1.2 * sin(8 * pi * X[, 2]) * cos(pi * X[, 1])
      if (d >= 3L) s <- s + 0.5 * cos(2 * pi * X[, 3])
      if (d >= 4L) {
        for (j in 4:d) {
          s <- s + (0.35 / sqrt(j)) * sin(2 * pi * X[, j])
        }
      }
      s
    }
  )
  sigma <- switch(
    scenario,
    smooth_smooth = 0.45,
    strong_aniso = 0.55
  )

  f <- f_fun(X)
  f_test <- f_fun(X_test)
  y <- .tt_with_preserved_seed({
    set.seed(seed + 17L)
    f + stats::rnorm(n, 0, sigma)
  })
  y_test <- .tt_with_preserved_seed({
    set.seed(seed + 19L)
    f_test + stats::rnorm(n_test, 0, sigma)
  })

  list(
    scenario = scenario,
    dataset_id = sprintf("d%d_%s_n%d_k%d_s%d", d, scenario, n, k, seed),
    n = n,
    n_test = n_test,
    d = d,
    k = k,
    degree = as.integer(degree),
    penalty_order = as.integer(penalty_order),
    seed = seed,
    sigma = as.numeric(sigma),
    rough_margin = if (identical(scenario, "strong_aniso")) 2L else NA_integer_,
    X = X,
    y = y,
    f = f,
    X_test = X_test,
    y_test = y_test,
    f_test = f_test,
    truth_note = paste0("d=", d, "; scenario=", scenario, "; sigma=", sigma)
  )
}

#' Shared Phase-1 ALS control (lab defaults).
#' @keywords internal
#' @noRd
.tt_lab_phase1_control <- function(seed = 1L, max_sweeps = 40L, tol = 1e-8,
                                   margin_order = NULL) {
  tt_control(
    max_sweeps = as.integer(max_sweeps),
    tol = as.numeric(tol),
    compute_edf = FALSE,
    seed = as.integer(seed),
    backend = "R",
    warn_lambda_boundary = FALSE,
    cgcv_trace = TRUE,
    cgcv_margin_order = margin_order
  )
}

#' Evaluate TT-gGCV with CRN probes, optional cache, cold or continuation init.
#'
#' @param init_policy `"cold_common"` (clone fixed cores) or `"continuation"`
#'   (use `init_cores` from neighbour; still clones).
#' @keywords internal
#' @noRd
.tt_lab_eval_ggcv <- function(lambda,
                              design,
                              rank,
                              probes,
                              M,
                              init_cores,
                              init_policy = c("cold_common", "continuation"),
                              mc_bank_id = "bank0",
                              cache = NULL,
                              epsilon_rel = 1e-3,
                              scheme = "forward",
                              control = NULL,
                              gdf_warm_start = TRUE,
                              gdf_init = NULL,
                              fit_backend = c("R", "Rcpp_fixed"),
                              on_nonconverged = "na") {
  init_policy <- match.arg(init_policy)
  fit_backend <- .tt_lab_match_fit_backend(fit_backend)
  gdf_init <- .tt_lab_match_gdf_init(gdf_init = gdf_init, warm_start = gdf_warm_start)
  theta <- log10(as.numeric(lambda))
  key <- .tt_lab_cache_key(
    design$dataset_id, rank, theta, init_policy, mc_bank_id, M = M
  )
  if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
    out <- cache[[key]]
    out$cache_hit <- TRUE
    return(out)
  }

  ctrl <- control %||% .tt_lab_phase1_control(seed = design$seed)
  init_use <- if (is.null(init_cores)) NULL else .tt_clone_cores(init_cores)

  # Isolate RNG around ALS init side-effects
  ev <- .tt_with_preserved_seed({
    tt_global_gcv(
      lambda = lambda,
      y = design$y,
      X = design$X,
      rank = rank,
      probes = probes,
      M = M,
      epsilon_rel = epsilon_rel,
      scheme = scheme,
      control = ctrl,
      init = init_use,
      fit_backend = fit_backend,
      k = design$k,
      degree = design$degree,
      penalty_order = design$penalty_order,
      on_nonconverged = on_nonconverged,
      gdf_init = gdf_init
    )
  })

  # Predictive metrics on held-out test (same X_test)
  rmse_test <- NA_real_
  ise_true <- NA_real_
  if (!is.null(ev$fit) && inherits(ev$fit, "ttpspline")) {
    mu_te <- tryCatch(
      as.numeric(predict(ev$fit, newdata = design$X_test, type = "response")),
      error = function(e) NA_real_
    )
    if (length(mu_te) == design$n_test && all(is.finite(mu_te))) {
      rmse_test <- sqrt(mean((mu_te - design$y_test)^2))
      ise_true <- mean((mu_te - design$f_test)^2)
    }
  }

  # Coerce scalars so grid data.frame() never sees NULL/length-0 fields
  as_num1 <- function(x) {
    x <- tryCatch(as.numeric(x)[1L], error = function(e) NA_real_)
    if (length(x) != 1L || is.null(x)) NA_real_ else x
  }
  as_int1 <- function(x) {
    x <- tryCatch(as.integer(x)[1L], error = function(e) NA_integer_)
    if (length(x) != 1L || is.null(x)) NA_integer_ else x
  }
  out <- list(
    global_gcv = as_num1(ev$global_gcv),
    rss = as_num1(ev$rss),
    gdf = as_num1(ev$gdf),
    gdf_mc_se = as_num1(ev$gdf_mc_se),
    lambda = as.numeric(ev$lambda %||% lambda),
    theta = as.numeric(theta),
    rank = as_int1(ev$rank %||% rank),
    valid = isTRUE(ev$valid),
    invalid_reasons = ev$invalid_reasons %||% character(0),
    converged = isTRUE(ev$fit_ok$ok %||% ev$fit$converged),
    n_sweeps = as_int1(ev$fit$n_sweeps %||% ev$fit_ok$n_sweeps),
    fit_time = as_num1(ev$fit_time),
    init_policy = init_policy,
    mc_bank_id = mc_bank_id,
    M = as.integer(M)[1L],
    rmse_test = as_num1(rmse_test),
    ise_true = as_num1(ise_true),
    cores = if (!is.null(ev$fit) && inherits(ev$fit, "ttpspline")) {
      .tt_clone_cores(ev$fit$cores)
    } else {
      NULL
    },
    fit = ev$fit,
    cache_hit = FALSE,
    key = key
  )
  if (!is.null(cache)) cache[[key]] <- out
  out
}

#' Full-TP GCV on a design at fixed lambda (exact EDF).
#' @keywords internal
#' @noRd
.tt_lab_eval_full_gcv <- function(design, lambda) {
  bs <- build_marginal_bases(
    design$X, k = design$k, degree = design$degree,
    knots = NULL, cyclic = NULL, period = NULL
  )
  .tt_lab_full_tp_gaussian(
    design$y, bs$basis, lambda = lambda,
    penalty_order = design$penalty_order
  )
}

#' Evaluate a log10-lambda grid (cold_common or continuation serpentine).
#' @keywords internal
#' @noRd
.tt_lab_eval_lambda_grid <- function(design,
                                     rank,
                                     log10_grid,
                                     probes,
                                     M,
                                     init_policy = c("cold_common", "continuation"),
                                     traverse = c("row_asc", "row_desc", "serpentine"),
                                     common_init = NULL,
                                     mc_bank_id = "bank0",
                                     cache = NULL,
                                     control = NULL,
                                     epsilon_rel = 1e-3) {
  init_policy <- match.arg(init_policy)
  traverse <- match.arg(traverse)
  g1 <- as.numeric(log10_grid)
  g2 <- as.numeric(log10_grid)
  coords <- expand.grid(log10_l1 = g1, log10_l2 = g2, KEEP.OUT.ATTRS = FALSE)
  # Order rows for continuation
  if (identical(traverse, "row_asc")) {
    ord <- order(coords$log10_l1, coords$log10_l2)
  } else if (identical(traverse, "row_desc")) {
    ord <- order(-coords$log10_l1, coords$log10_l2)
  } else {
    # serpentine within each l1 row
    split_idx <- split(seq_len(nrow(coords)), coords$log10_l1)
    ord <- integer(0)
    flip <- FALSE
    for (nm in names(split_idx)) {
      ii <- split_idx[[nm]]
      ii <- ii[order(coords$log10_l2[ii])]
      if (flip) ii <- rev(ii)
      ord <- c(ord, ii)
      flip <- !flip
    }
  }
  coords <- coords[ord, , drop = FALSE]

  if (identical(init_policy, "cold_common") && is.null(common_init)) {
    common_init <- .tt_with_preserved_seed({
      tt_initialize(d = 2L, rank = rank, k = design$k, seed = design$seed)
    })
  }

  rows <- vector("list", nrow(coords))
  prev_cores <- NULL
  for (i in seq_len(nrow(coords))) {
    lam <- 10^c(coords$log10_l1[i], coords$log10_l2[i])
    init_cores <- if (identical(init_policy, "cold_common")) {
      common_init
    } else {
      prev_cores
    }
    t0 <- proc.time()[["elapsed"]]
    ev <- .tt_lab_eval_ggcv(
      lambda = lam, design = design, rank = rank, probes = probes, M = M,
      init_cores = init_cores, init_policy = init_policy,
      mc_bank_id = mc_bank_id, cache = cache, control = control,
      epsilon_rel = epsilon_rel
    )
    full <- .tt_lab_eval_full_gcv(design, lam)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (!is.null(ev$cores)) prev_cores <- ev$cores
    rows[[i]] <- data.frame(
      dataset_id = design$dataset_id,
      scenario = design$scenario,
      rank = as.integer(rank)[[1L]],
      init_policy = init_policy,
      traverse = traverse,
      log10_l1 = coords$log10_l1[i],
      log10_l2 = coords$log10_l2[i],
      lambda1 = lam[1],
      lambda2 = lam[2],
      rss_tt = as.numeric(ev$rss)[1],
      gdf = as.numeric(ev$gdf)[1],
      gdf_se = as.numeric(ev$gdf_mc_se %||% NA_real_)[1],
      global_gcv_tt = as.numeric(ev$global_gcv)[1],
      valid = isTRUE(ev$valid),
      converged = isTRUE(ev$converged),
      n_sweeps = as.integer(ev$n_sweeps %||% NA_integer_)[1],
      rss_full = as.numeric(full$rss)[1],
      edf_full = as.numeric(full$edf)[1],
      gcv_full = as.numeric(full$gcv)[1],
      rmse_test = as.numeric(ev$rmse_test)[1],
      ise_true = as.numeric(ev$ise_true)[1],
      elapsed_sec = as.numeric(elapsed)[1],
      cache_hit = isTRUE(ev$cache_hit),
      M = as.integer(M)[1],
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, rows)
  tab$flat_1pct <- .tt_lab_flatness_mask(tab$global_gcv_tt, 0.01)
  tab
}

#' Derivative-free joint optimization of TT-gGCV in log10-lambda.
#' @keywords internal
#' @noRd
.tt_lab_optimize_ggcv <- function(design,
                                  rank,
                                  probes,
                                  M_search,
                                  M_final,
                                  starts,
                                  box = c(-5, 5),
                                  common_init = NULL,
                                  mc_bank_id = "bank0",
                                  cache = NULL,
                                  control = NULL,
                                  epsilon_rel = 1e-3,
                                  method = "Nelder-Mead",
                                  optim_control = list(maxit = 40, reltol = 1e-3)) {
  box <- as.numeric(box)
  if (is.null(common_init)) {
    common_init <- .tt_with_preserved_seed({
      tt_initialize(d = 2L, rank = rank, k = design$k, seed = design$seed)
    })
  }
  # Objective with box projection
  obj_fun <- function(theta) {
    th <- pmin(pmax(as.numeric(theta), box[1]), box[2])
    lam <- 10^th
    ev <- .tt_lab_eval_ggcv(
      lambda = lam, design = design, rank = rank, probes = probes,
      M = M_search, init_cores = common_init, init_policy = "cold_common",
      mc_bank_id = mc_bank_id, cache = cache, control = control,
      epsilon_rel = epsilon_rel
    )
    val <- ev$global_gcv
    if (!is.finite(val) || !isTRUE(ev$valid)) return(1e12)
    val
  }

  starts <- lapply(starts, function(s) pmin(pmax(as.numeric(s), box[1]), box[2]))
  results <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    t0 <- proc.time()[["elapsed"]]
    op <- tryCatch(
      stats::optim(
        par = starts[[i]], fn = obj_fun, method = method,
        control = optim_control
      ),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    if (inherits(op, "error")) {
      results[[i]] <- data.frame(
        start_id = i, ok = FALSE, reason = conditionMessage(op),
        log10_l1 = starts[[i]][1], log10_l2 = starts[[i]][2],
        global_gcv_search = NA_real_, global_gcv_final = NA_real_,
        gdf = NA_real_, gdf_se = NA_real_, rss = NA_real_,
        rmse_test = NA_real_, ise_true = NA_real_,
        elapsed_sec = elapsed, stringsAsFactors = FALSE
      )
      next
    }
    th <- pmin(pmax(op$par, box[1]), box[2])
    # Final re-eval with larger M
    evf <- .tt_lab_eval_ggcv(
      lambda = 10^th, design = design, rank = rank, probes = probes,
      M = M_final, init_cores = common_init, init_policy = "cold_common",
      mc_bank_id = mc_bank_id, cache = cache, control = control,
      epsilon_rel = epsilon_rel
    )
    results[[i]] <- data.frame(
      start_id = i, ok = TRUE, reason = NA_character_,
      log10_l1 = th[1], log10_l2 = th[2],
      lambda1 = 10^th[1], lambda2 = 10^th[2],
      global_gcv_search = op$value,
      global_gcv_final = evf$global_gcv,
      gdf = evf$gdf, gdf_se = evf$gdf_mc_se, rss = evf$rss,
      rmse_test = evf$rmse_test, ise_true = evf$ise_true,
      converged_fit = evf$converged, n_sweeps = evf$n_sweeps,
      optim_counts = op$counts[["function"]] %||% NA_integer_,
      optim_convergence = op$convergence,
      elapsed_sec = elapsed,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, results)
}
