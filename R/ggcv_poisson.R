# Poisson / GLM TT-gGCV
#
# Two modes (tt_control$ggcv_glm_mode):
#   "working"     — pilot PIRLS → freeze (z,W) → Gaussian TT-gGCV on row-scaled
#                   WLS problem → refit Poisson at that λ.
#                   ponytail: ceiling = ignores W(y,λ); upgrade = algorithmic.
#   "algorithmic" — Hutchinson GDF of y |-> μ̂(y;λ) under fixed-λ PIRLS+ALS,
#                   GCV = n * D(y,μ̂) / (n - GDF)^2. Expensive.

#' Poisson TT-gGCV (working-response or algorithmic Hutchinson)
#'
#' @inheritParams tt_ggcv
#' @param family Must be Poisson (binomial reserved).
#' @param offset Optional offset (length n or scalar).
#' @param cyclic,period,knots Forwarded to [ttps()].
#' @param mode `"working"` or `"algorithmic"` (default from
#'   `control$ggcv_glm_mode`, else `"working"`).
#' @param anisotropic If `FALSE` (default for algorithmic when `d > 2`),
#'   search a shared isotropic λ; if `TRUE`, Sobol box in log10-λ^d.
#' @return Same shape as [tt_ggcv()]: `lambda`, `gcv`, `gdf`, `fit`, plus
#'   `ggcv_glm_mode` / diagnostics.
#' @export
tt_ggcv_poisson <- function(y,
                            X,
                            rank,
                            family = stats::poisson(),
                            offset = NULL,
                            cyclic = NULL,
                            period = NULL,
                            knots = NULL,
                            mode = NULL,
                            anisotropic = NULL,
                            k = 8L,
                            degree = 3L,
                            penalty_order = 2L,
                            control = tt_control(max_sweeps = 40L,
                                                 pirls_maxit = 25L,
                                                 compute_edf = FALSE,
                                                 seed = 1L),
                            n_global = NULL,
                            n_refine = NULL,
                            M_search = NULL,
                            M_final = NULL,
                            seed = NULL,
                            include_cgcv_anchor = NULL,
                            verbose = FALSE) {
  fam <- normalize_family(family)
  key <- family_key(fam)
  if (!identical(key, "poisson")) {
    stop("tt_ggcv_poisson() currently supports family = poisson() only.",
         call. = FALSE)
  }
  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y)
  d <- ncol(X)
  if (nrow(X) != n) stop("nrow(X) must equal length(y).", call. = FALSE)
  if (any(y < 0)) stop("poisson requires non-negative y.", call. = FALSE)

  mode <- mode %||% control$ggcv_glm_mode %||% "working"
  mode <- match.arg(mode, c("working", "algorithmic"))
  if (is.null(anisotropic)) {
    anisotropic <- isTRUE(control$ggcv_poisson_anisotropic %||% (d <= 2L))
  }
  anisotropic <- isTRUE(anisotropic)

  n_global <- as.integer(n_global %||% control$ggcv_n_global %||%
                           max(8L, 4L * d))
  n_refine <- as.integer(n_refine %||% control$ggcv_n_refine %||% 3L)
  M_search <- as.integer(M_search %||% control$ggcv_M_search %||% 8L)
  M_final <- as.integer(M_final %||% control$ggcv_M_final %||% 16L)
  seed <- as.integer(seed %||% control$seed %||% 1L)
  include_cgcv_anchor <- include_cgcv_anchor %||%
    isTRUE(control$ggcv_include_cgcv_anchor %||% TRUE)

  ctrl <- control
  ctrl$compute_edf <- isTRUE(control$compute_edf)
  ctrl$seed <- seed

  if (identical(mode, "working")) {
    out <- .tt_ggcv_poisson_working(
      y = y, X = X, rank = rank, fam = fam, offset = offset,
      cyclic = cyclic, period = period, knots = knots,
      k = k, degree = degree, penalty_order = penalty_order,
      control = ctrl, n_global = n_global, n_refine = n_refine,
      M_search = M_search, M_final = M_final, seed = seed,
      include_cgcv_anchor = include_cgcv_anchor, verbose = verbose
    )
  } else {
    out <- .tt_ggcv_poisson_algorithmic(
      y = y, X = X, rank = rank, fam = fam, offset = offset,
      cyclic = cyclic, period = period, knots = knots,
      k = k, degree = degree, penalty_order = penalty_order,
      control = ctrl, n_global = n_global, n_refine = n_refine,
      M_search = M_search, M_final = M_final, seed = seed,
      include_cgcv_anchor = include_cgcv_anchor,
      anisotropic = anisotropic, verbose = verbose
    )
  }
  out$method <- "gGCV"
  out$ggcv_glm_mode <- mode
  if (!is.null(out$fit) && inherits(out$fit, c("ttps", "ttpspline"))) {
    out$fit$lambda_method <- "gGCV"
    out$fit$ggcv <- list(
      gcv = out$gcv,
      gdf = out$gdf,
      gdf_mc_se = out$gdf_mc_se,
      theta = out$theta,
      mode = mode,
      anisotropic = anisotropic,
      working = out$working,
      boundary = out$boundary,
      elapsed = out$elapsed
    )
  }
  out
}

#' @keywords internal
#' @noRd
.ttps_dispatch_ggcv_poisson <- function(y, X, rank, k, degree, penalty_order,
                                        control, init = NULL, cl = NULL,
                                        offset = NULL, cyclic = NULL,
                                        period = NULL, knots = NULL) {
  opt <- tt_ggcv_poisson(
    y = y, X = X, rank = rank,
    family = stats::poisson(),
    offset = offset, cyclic = cyclic, period = period, knots = knots,
    k = k, degree = degree, penalty_order = penalty_order,
    control = control,
    verbose = isTRUE(control$trace)
  )
  fit <- opt$fit
  if (is.null(fit) || !inherits(fit, c("ttps", "ttpspline"))) {
    stop("lambda = \"gGCV\" (poisson) failed to produce a TT fit.", call. = FALSE)
  }
  if (!is.null(cl)) fit$call <- cl
  fit$lambda_method <- "gGCV"
  fit$ggcv <- opt$fit$ggcv %||% list(
    gcv = opt$gcv, gdf = opt$gdf, mode = opt$ggcv_glm_mode
  )
  fit
}

# ── working-response path ────────────────────────────────────────────────────

#' @keywords internal
#' @noRd
.tt_ggcv_poisson_working <- function(y, X, rank, fam, offset, cyclic, period,
                                     knots, k, degree, penalty_order, control,
                                     n_global, n_refine, M_search, M_final,
                                     seed, include_cgcv_anchor, verbose) {
  t0 <- proc.time()[["elapsed"]]
  d <- ncol(X)
  # Pilot: cGCV or fixed λ_start to get a working (z, W).
  if (isTRUE(include_cgcv_anchor)) {
    if (isTRUE(verbose)) message("poisson gGCV[working]: pilot cGCV …")
    pilot <- ttps(
      y = y, X = X, family = fam, offset = offset,
      rank = rank, k = k, degree = degree, penalty_order = penalty_order,
      lambda = "cGCV", cyclic = cyclic, period = period, knots = knots,
      control = control
    )
  } else {
    lam0 <- control$lambda_start %||% 1
    if (length(lam0) == 1L) lam0 <- rep(as.numeric(lam0), d)
    if (isTRUE(verbose)) message("poisson gGCV[working]: pilot fixed λ …")
    pilot <- ttps(
      y = y, X = X, family = fam, offset = offset,
      rank = rank, k = k, degree = degree, penalty_order = penalty_order,
      lambda = lam0, cyclic = cyclic, period = period, knots = knots,
      control = control
    )
  }
  eta <- as.numeric(stats::predict(pilot, type = "link"))
  work <- glm_working(fam, y, eta, control = control)
  sw <- sqrt(pmax(as.numeric(work$weight), 1e-12))
  y_g <- sw * as.numeric(work$z)
  X_g <- X * sw # row-scale: unweighted Gaussian ≡ WLS on (z,W)

  # Drop cyclic for the Gaussian proxy design (hour already in [0,1] scale);
  # final Poisson refit restores cyclic/period.
  if (isTRUE(verbose)) message("poisson gGCV[working]: Gaussian TT-gGCV on WLS …")
  ctrl_g <- control
  ctrl_g$ggcv_n_global <- n_global
  ctrl_g$ggcv_n_refine <- n_refine
  ctrl_g$ggcv_M_search <- M_search
  ctrl_g$ggcv_M_final <- M_final
  ctrl_g$ggcv_include_cgcv_anchor <- TRUE
  opt_g <- tt_ggcv(
    y = y_g, X = X_g, rank = rank, family = stats::gaussian(),
    k = k, degree = degree, penalty_order = penalty_order,
    n_global = n_global, n_refine = n_refine,
    M_search = M_search, M_final = M_final, seed = seed,
    control = ctrl_g, include_cgcv_anchor = TRUE,
    verbose = verbose
  )
  lam <- as.numeric(opt_g$lambda)
  if (length(lam) == 1L) lam <- rep(lam, d)

  if (isTRUE(verbose)) message("poisson gGCV[working]: final Poisson refit …")
  fit <- ttps(
    y = y, X = X, family = fam, offset = offset,
    rank = rank, k = k, degree = degree, penalty_order = penalty_order,
    lambda = lam, cyclic = cyclic, period = period, knots = knots,
    control = control
  )
  mu <- as.numeric(stats::fitted(fit))
  dev <- glm_deviance(fam, y, mu)
  list(
    lambda = lam,
    theta = log10(lam),
    gcv = opt_g$gcv,
    gdf = opt_g$gdf,
    gdf_mc_se = opt_g$gdf_mc_se,
    fit = fit,
    pilot = pilot,
    working = list(deviance = dev, gaussian_gcv = opt_g$gcv,
                   gaussian_gdf = opt_g$gdf, note = "GCV/GDF from WLS proxy"),
    boundary = opt_g$boundary,
    winner_source = opt_g$winner_source,
    elapsed = proc.time()[["elapsed"]] - t0
  )
}

# ── algorithmic Hutchinson path ──────────────────────────────────────────────

#' Fixed-λ Poisson TT fit (ALS/PIRLS).
#' @keywords internal
#' @noRd
.tt_pois_fit_fixed <- function(y, X, lambda, rank, fam, offset, cyclic, period,
                               knots, k, degree, penalty_order, control,
                               init = NULL) {
  d <- ncol(as.matrix(X))
  lam <- as.numeric(lambda)
  if (length(lam) == 1L) lam <- rep(lam, d)
  ttps(
    y = as.numeric(y), X = as.matrix(X), family = fam, offset = offset,
    rank = rank, k = k, degree = degree, penalty_order = penalty_order,
    lambda = lam, cyclic = cyclic, period = period, knots = knots,
    init = init, control = control
  )
}

#' Evaluate Poisson algorithmic TT-gGCV at one λ.
#' @keywords internal
#' @noRd
.tt_pois_ggcv_eval <- function(lambda, y, X, rank, fam, offset, cyclic, period,
                               knots, k, degree, penalty_order, control,
                               probes, scheme = "forward", epsilon_rel = 1e-3,
                               init = NULL) {
  n <- length(y)
  d <- ncol(X)
  lam <- as.numeric(lambda)
  if (length(lam) == 1L) lam <- rep(lam, d)
  fit <- tryCatch(
    .tt_pois_fit_fixed(
      y, X, lam, rank, fam, offset, cyclic, period, knots,
      k, degree, penalty_order, control, init = init
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(gcv = Inf, gdf = NA_real_, gdf_mc_se = NA_real_,
                deviance = NA_real_, fit = NULL, valid = FALSE,
                reason = conditionMessage(fit)))
  }
  mu0 <- as.numeric(stats::fitted(fit))
  dev <- glm_deviance(fam, y, mu0)
  # Hutchinson GDF of μ̂(y): cold-ish probe warm from cloned cores
  eps <- max(as.numeric(epsilon_rel) * stats::sd(y), 1e-3, mean(y) * 1e-4)
  if (!is.finite(eps) || eps <= 0) eps <- 1e-3
  M <- ncol(probes)
  contrib <- rep(NA_real_, M)
  base_cores <- if (!is.null(fit$cores)) .tt_clone_cores(fit$cores) else NULL
  for (j in seq_len(M)) {
    zj <- probes[, j]
    if (identical(scheme, "central")) {
      y_p <- pmax(y + eps * zj, 0)
      y_m <- pmax(y - eps * zj, 0)
      fit_p <- tryCatch(
        .tt_pois_fit_fixed(
          y_p, X, lam, rank, fam, offset, cyclic, period, knots,
          k, degree, penalty_order, control, init = base_cores
        ),
        error = function(e) NULL
      )
      fit_m <- tryCatch(
        .tt_pois_fit_fixed(
          y_m, X, lam, rank, fam, offset, cyclic, period, knots,
          k, degree, penalty_order, control, init = base_cores
        ),
        error = function(e) NULL
      )
      if (is.null(fit_p) || is.null(fit_m)) next
      mu_p <- as.numeric(stats::fitted(fit_p))
      mu_m <- as.numeric(stats::fitted(fit_m))
      contrib[j] <- sum(zj * (mu_p - mu_m) / (2 * eps))
    } else {
      y_p <- pmax(y + eps * zj, 0)
      fit_p <- tryCatch(
        .tt_pois_fit_fixed(
          y_p, X, lam, rank, fam, offset, cyclic, period, knots,
          k, degree, penalty_order, control, init = base_cores
        ),
        error = function(e) NULL
      )
      if (is.null(fit_p)) next
      mu_p <- as.numeric(stats::fitted(fit_p))
      contrib[j] <- sum(zj * (mu_p - mu0) / eps)
    }
  }
  ok <- is.finite(contrib)
  gdf <- if (any(ok)) mean(contrib[ok]) else NA_real_
  gdf_se <- if (sum(ok) >= 2L) stats::sd(contrib[ok]) / sqrt(sum(ok)) else NA_real_
  denom <- (n - gdf)^2
  gcv <- if (is.finite(gdf) && is.finite(denom) && denom > 1e-12 &&
              gdf > -1 && gdf < n) {
    n * dev / denom
  } else {
    Inf
  }
  list(
    gcv = gcv, gdf = gdf, gdf_mc_se = gdf_se, deviance = dev,
    fit = fit, valid = is.finite(gcv), epsilon = eps,
    M_ok = sum(ok), contrib = contrib
  )
}

#' @keywords internal
#' @noRd
.tt_ggcv_poisson_algorithmic <- function(y, X, rank, fam, offset, cyclic, period,
                                         knots, k, degree, penalty_order, control,
                                         n_global, n_refine, M_search, M_final,
                                         seed, include_cgcv_anchor,
                                         anisotropic, verbose) {
  t0 <- proc.time()[["elapsed"]]
  n <- length(y)
  d <- ncol(X)
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  lo <- log10(bounds[1L])
  hi <- log10(bounds[2L])

  # Shrink PIRLS/ALS for search evaluations
  ctrl_search <- control
  ctrl_search$pirls_maxit <- min(as.integer(control$pirls_maxit %||% 25L), 15L)
  ctrl_search$max_sweeps <- min(as.integer(control$max_sweeps %||% 40L), 12L)
  ctrl_search$compute_edf <- FALSE
  ctrl_final <- control
  ctrl_final$compute_edf <- isTRUE(control$compute_edf)

  probes_search <- .tt_lab_rademacher_probes(n, M_search, probe_seed = seed + 17L)
  probes_final <- .tt_lab_rademacher_probes(n, M_final, probe_seed = seed + 17L)

  eval_at <- function(lam, probes, ctrl, init = NULL) {
    .tt_pois_ggcv_eval(
      lambda = lam, y = y, X = X, rank = rank, fam = fam, offset = offset,
      cyclic = cyclic, period = period, knots = knots, k = k, degree = degree,
      penalty_order = penalty_order, control = ctrl, probes = probes,
      scheme = "forward", epsilon_rel = 1e-3, init = init
    )
  }

  # Optional cGCV anchor
  candidates <- list()
  if (isTRUE(include_cgcv_anchor)) {
    if (isTRUE(verbose)) message("poisson gGCV[algorithmic]: cGCV anchor …")
    anchor <- tryCatch(
      ttps(
        y = y, X = X, family = fam, offset = offset, rank = rank, k = k,
        degree = degree, penalty_order = penalty_order, lambda = "cGCV",
        cyclic = cyclic, period = period, knots = knots, control = control
      ),
      error = function(e) NULL
    )
    if (!is.null(anchor)) {
      ev <- eval_at(anchor$lambda, probes_search, ctrl_search)
      candidates[[length(candidates) + 1L]] <- list(
        lambda = as.numeric(anchor$lambda), gcv = ev$gcv, source = "cGCV_anchor",
        ev = ev
      )
    }
  }

  if (!isTRUE(anisotropic)) {
    # Isotropic 1D search on log10(λ)
    if (isTRUE(verbose)) message("poisson gGCV[algorithmic]: isotropic Brent …")
    grid <- seq(lo, hi, length.out = max(7L, n_global))
    for (th in grid) {
      lam <- rep(10^th, d)
      ev <- eval_at(lam, probes_search, ctrl_search)
      candidates[[length(candidates) + 1L]] <- list(
        lambda = lam, gcv = ev$gcv, source = "iso_grid", ev = ev
      )
    }
    best <- candidates[[which.min(vapply(candidates, function(z) z$gcv, 1))]]
    # Local refine
    if (n_refine > 0L && is.finite(best$gcv)) {
      th0 <- mean(log10(best$lambda))
      opt <- tryCatch(
        stats::optimize(
          function(th) {
            ev <- eval_at(rep(10^th, d), probes_search, ctrl_search)
            ev$gcv
          },
          interval = c(max(lo, th0 - 1), min(hi, th0 + 1))
        ),
        error = function(e) NULL
      )
      if (!is.null(opt)) {
        lam <- rep(10^opt$minimum, d)
        ev <- eval_at(lam, probes_final, ctrl_final)
        candidates[[length(candidates) + 1L]] <- list(
          lambda = lam, gcv = ev$gcv, source = "iso_refine", ev = ev
        )
      }
    }
  } else {
    if (isTRUE(verbose)) {
      message(sprintf("poisson gGCV[algorithmic]: Sobol d=%d n=%d …", d, n_global))
    }
    box <- .tt_lab_sobol_box(n_global, rep(lo, d), rep(hi, d), skip = 1L)
    for (i in seq_len(nrow(box))) {
      lam <- 10^as.numeric(box[i, ])
      ev <- eval_at(lam, probes_search, ctrl_search)
      candidates[[length(candidates) + 1L]] <- list(
        lambda = lam, gcv = ev$gcv, source = "sobol", ev = ev
      )
    }
  }

  ok_gcv <- vapply(candidates, function(z) is.finite(z$gcv), logical(1))
  if (!any(ok_gcv)) {
    stop("poisson algorithmic gGCV: no finite GCV evaluations.", call. = FALSE)
  }
  best <- candidates[[which.min(vapply(candidates, function(z) {
    if (is.finite(z$gcv)) z$gcv else Inf
  }, 1))]]

  # Final evaluation at winner with fuller probe bank
  if (isTRUE(verbose)) message("poisson gGCV[algorithmic]: final eval …")
  final <- eval_at(best$lambda, probes_final, ctrl_final)
  fit <- final$fit
  if (is.null(fit)) {
    fit <- .tt_pois_fit_fixed(
      y, X, best$lambda, rank, fam, offset, cyclic, period, knots,
      k, degree, penalty_order, ctrl_final
    )
  }

  list(
    lambda = as.numeric(best$lambda),
    theta = log10(as.numeric(best$lambda)),
    gcv = final$gcv,
    gdf = final$gdf,
    gdf_mc_se = final$gdf_mc_se,
    fit = fit,
    working = NULL,
    boundary = any(abs(log10(best$lambda) - lo) < 0.05) ||
      any(abs(log10(best$lambda) - hi) < 0.05),
    winner_source = best$source,
    anisotropic = anisotropic,
    search = data.frame(
      source = vapply(candidates, `[[`, "", "source"),
      gcv = vapply(candidates, function(z) z$gcv, 1),
      stringsAsFactors = FALSE
    ),
    elapsed = proc.time()[["elapsed"]] - t0
  )
}
