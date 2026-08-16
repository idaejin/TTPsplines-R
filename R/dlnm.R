#' Tensor-Train DLNM (distributed lag non-linear model)
#'
#' Fits
#' \deqn{\eta_t = \mathrm{offset}_t + \alpha + Z_t\beta + s(t)
#'   + \sum_{\ell=0}^{L} h(x^{(1)}_{t-\ell},\ldots,x^{(p)}_{t-\ell},\ell)}
#' where \(h\) is a tensor-product P-spline with coefficient tensor \(\Theta\)
#' in TT format. The dense cross-basis \(W\) is **never** formed: each ALS /
#' PIRLS step contracts the TT against lag-specific marginal bases and sums
#' over \(\ell\).
#'
#' Calendar confounding belongs in `smooth=` / `linear=` (e.g.
#' `s(year)+s(month,cc)+dow`), not inside the lag surface.
#'
#' @param y Response vector (length \(n\)).
#' @param x Named list of exposure vectors (each length \(n\)), e.g.
#'   `list(temp = d$temp, pm10 = d$pm10, o3 = d$o3)`.
#' @param lag Maximum lag \(L\) (non-negative integer). Histories use
#'   \(\ell=0,\ldots,L\).
#' @param family A family object (`poisson()`, `gaussian()`, …).
#' @param rank TT rank (scalar or compatible chain).
#' @param k Basis size for each **exposure** margin.
#' @param k_lag Basis size for the lag margin (defaults to `k`).
#' @param degree B-spline degree.
#' @param penalty_order Difference penalty order.
#' @param lambda Fixed nonnegative scalar / length-`d` vector, or `"cGCV"`
#'   (conditional GCV on each TT core with summed lag designs).
#' @param linear Optional parametric design (`dow`, …); see [ttps()].
#' @param smooth Optional additive smooths (e.g. year + month); see [ttps()].
#' @param lambda_smooth Default for smooth terms.
#' @param offset,weights Optional offset / observation weights.
#' @param control A [tt_control()] object.
#' @param ... Ignored.
#'
#' @return An object of class `c("ttps_dlnm", "ttpspline")` with extra
#'   components `dlnm` (lag histories, knots, names) for prediction.
#'
#' @examples
#' \dontrun{
#' data(chicagoNMMAPS, package = "dlnm")
#' d <- na.omit(chicagoNMMAPS[, c("death","temp","pm10","o3","year","month","dow")])
#' fit <- ttps_dlnm(
#'   d$death,
#'   list(temp = d$temp, pm10 = d$pm10, o3 = d$o3),
#'   lag = 14, family = poisson(), rank = 1, k = 6, k_lag = 5,
#'   lambda = 1,
#'   linear = model.matrix(~ 0 + dow, d),
#'   smooth = list(
#'     year = list(x = d$year, bs = "ps", k = 10, m = 2, target_edf = 6),
#'     month = list(x = d$month, bs = "cc", k = 12, m = 2,
#'                  period = c(0.5, 12.5), target_edf = 6)
#'   ),
#'   control = tt_control(max_sweeps = 8, pirls_maxit = 10, compute_edf = FALSE)
#' )
#' pr <- predict_dlnm(fit, var = "temp", at = seq(-10, 30, by = 2), cen = 21)
#' }
#' @export
ttps_dlnm <- function(y,
                      x,
                      lag = 14L,
                      family = stats::poisson(),
                      rank = 1L,
                      k = 8L,
                      k_lag = NULL,
                      degree = 3L,
                      penalty_order = 2L,
                      lambda = 1,
                      linear = NULL,
                      smooth = NULL,
                      lambda_smooth = "cGCV",
                      offset = NULL,
                      weights = NULL,
                      control = tt_control(),
                      ...) {
  cl <- match.call()
  fam <- normalize_family(family)
  key <- family_key(fam)
  y <- as.numeric(y)
  n0 <- length(y)
  if (!is.list(x) || length(x) < 1L) {
    stop("`x` must be a named list of exposure vectors.", call. = FALSE)
  }
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    stop("`x` must be a named list (e.g. list(temp=..., pm10=...)).",
         call. = FALSE)
  }
  for (nm in names(x)) {
    x[[nm]] <- as.numeric(x[[nm]])
    if (length(x[[nm]]) != n0) {
      stop("Exposure '", nm, "' length must equal length(y).", call. = FALSE)
    }
  }
  if (anyNA(y) || any(vapply(x, anyNA, logical(1)))) {
    stop("NA in y or exposures not supported; drop incomplete rows first.",
         call. = FALSE)
  }
  lag <- as.integer(lag)
  if (length(lag) != 1L || !is.finite(lag) || lag < 0L) {
    stop("`lag` must be a nonnegative integer.", call. = FALSE)
  }
  if (is.null(k_lag)) k_lag <- k
  k <- as.integer(k)
  k_lag <- as.integer(k_lag)
  degree <- as.integer(degree)
  penalty_order <- as.integer(penalty_order)

  if (!inherits(control, "tt_control")) {
    control <- do.call(tt_control, as.list(control))
  }

  dlnm <- build_tt_dlnm(
    x = x, lag = lag, k = k, k_lag = k_lag, degree = degree
  )
  ok <- dlnm$ok
  n <- sum(ok)
  y_ok <- y[ok]
  offset_ok <- normalize_offset(offset, n0)[ok]
  weights_ok <- normalize_weights(weights, n0)[ok]
  linear_ok <- if (is.null(linear)) NULL else normalize_linear(linear, n0)[ok, , drop = FALSE]
  if (!is.null(linear_ok)) attr(linear_ok, "ttps_linear_ok") <- TRUE
  smooth_ok <- NULL
  if (!is.null(smooth)) {
    # Build on full data then subset rows to `ok`
    smooth_full <- normalize_smooth(smooth, n0, lambda_smooth = lambda_smooth)
    smooth_ok <- subset_smooth(smooth_full, which(ok))
  }

  p_exp <- length(x)
  d_tt <- p_exp + 1L
  ranks <- tt_rank(rank, d = d_tt)
  # Margin sizes: exposures k, lag k_lag
  p_vec <- c(rep(k, p_exp), k_lag)
  lambda_spec <- parse_lambda_spec(lambda, d = d_tt, control = control)

  basis_lags <- dlnm$basis_lags
  cyclic <- rep(FALSE, d_tt)
  attr(basis_lags, "cyclic") <- cyclic
  # Attach per-margin sizes for penalties
  attr(basis_lags, "p_vec") <- p_vec

  t0 <- proc.time()[["elapsed"]]
  if (identical(key, "gaussian")) {
    raw <- tt_dlnm_als_fit(
      y_ok, basis_lags, ranks, lambda_spec, control, penalty_order,
      p_vec = p_vec, offset = offset_ok, weights = weights_ok,
      linear = linear_ok, smooth = smooth_ok
    )
  } else if (key %in% c("poisson", "bernoulli")) {
    raw <- tt_dlnm_pirls_fit(
      y_ok, basis_lags, fam, ranks, lambda_spec, control, penalty_order,
      p_vec = p_vec, offset = offset_ok, weights = weights_ok,
      linear = linear_ok, smooth = smooth_ok
    )
  } else {
    stop("ttps_dlnm supports gaussian / poisson / bernoulli in v0.", call. = FALSE)
  }
  elapsed <- proc.time()[["elapsed"]] - t0

  f_tt <- tt_dlnm_contraction(raw$cores, basis_lags)
  eta <- as.numeric(
    offset_ok + raw$intercept +
      tt_linear_contrib(raw$linear %||% linear_ok, raw$beta) +
      tt_smooth_contrib(raw$smooth %||% smooth_ok) +
      f_tt
  )
  mu <- invlink_eta(fam, eta)

  # Working sample = complete lag windows (length n). Also keep full-length
  # copies with NA burn-in for aligned time-series plots.
  eta_full <- rep(NA_real_, n0)
  mu_full <- rep(NA_real_, n0)
  eta_full[ok] <- eta
  mu_full[ok] <- mu

  npar_tt <- sum(vapply(raw$cores, length, integer(1)))
  npar_dense <- as.numeric(prod(p_vec))
  npar_intr <- npar_tt - tt_gauge_dim(ranks)

  out <- structure(
    list(
      call = cl,
      family = fam,
      family_key = key,
      y = y_ok,
      y_full = y,
      n = n,
      n_full = n0,
      d = d_tt,
      p_exp = p_exp,
      k = k,
      k_lag = k_lag,
      p_vec = p_vec,
      degree = degree,
      penalty_order = penalty_order,
      rank = ranks,
      cores = raw$cores,
      intercept = raw$intercept,
      beta = raw$beta,
      linear = raw$linear %||% linear_ok,
      smooth = raw$smooth %||% smooth_ok,
      lambda = raw$lambda,
      lambda_method = raw$method_lambda %||% lambda_spec$method,
      fitted.values = mu,
      fitted.values_full = mu_full,
      linear.predictors = eta,
      linear.predictors_full = eta_full,
      residuals = y_ok - mu,
      deviance = raw$deviance %||% sum(fam$dev.resids(y_ok, mu, rep(1, n))),
      edf = raw$edf %||% NA_real_,
      edf_note = if (!is.finite(raw$edf %||% NA_real_)) {
        "joint EDF not yet computed for TT-DLNM"
      } else {
        NULL
      },
      npar_tt = npar_tt,
      npar_tt_intrinsic = npar_intr,
      npar_dense = npar_dense,
      compression_ratio = npar_dense / max(npar_tt, 1),
      optimizer = raw$optimizer %||% "DLNM-ALS",
      optimizer_used = raw$optimizer %||% "DLNM-ALS",
      backend = "R",
      control = control,
      converged = isTRUE(raw$converged %||% TRUE),
      n_sweeps = raw$n_sweeps %||% NA_integer_,
      n_iter = raw$n_iter %||% raw$n_pirls %||% NA_integer_,
      n_pirls = raw$n_pirls %||% raw$n_iter %||% NA_integer_,
      n_eval = raw$n_eval %||% NA_integer_,
      time = elapsed,
      timing = elapsed,
      dlnm = dlnm,
      offset = offset_ok,
      weights = weights_ok,
      cyclic = cyclic,
      # No scattered X: prediction uses dlnm histories / predict_dlnm()
      X = NULL
    ),
    class = c("ttps_dlnm", "ttpspline")
  )
  out
}

# ---- construction -----------------------------------------------------------

#' Lag matrix (columns = lags) without tsModel dependency.
#' @keywords internal
#' @noRd
.tt_lag_matrix <- function(x, lags) {
  x <- as.numeric(x)
  n <- length(x)
  lags <- as.integer(lags)
  out <- matrix(NA_real_, n, length(lags))
  colnames(out) <- paste0("lag", lags)
  for (j in seq_along(lags)) {
    ell <- lags[[j]]
    if (ell == 0L) {
      out[, j] <- x
    } else if (ell < n) {
      out[(ell + 1L):n, j] <- x[seq_len(n - ell)]
    }
  }
  out
}

#' Build TT-DLNM lag histories, knots, and per-lag basis lists.
#' @keywords internal
#' @noRd
build_tt_dlnm <- function(x, lag, k, k_lag, degree = 3L) {
  nms <- names(x)
  p <- length(x)
  lag <- as.integer(lag)
  lags <- seq.int(0L, lag)
  n <- length(x[[1L]])
  Q <- lapply(x, function(v) .tt_lag_matrix(v, lags))
  names(Q) <- nms
  ok <- Reduce(`&`, lapply(Q, function(M) apply(M, 1L, function(r) all(is.finite(r)))))
  if (!any(ok)) stop("No complete lag windows; reduce `lag` or check data.", call. = FALSE)

  # Knots from observed (unlagged) exposure range + lag grid
  knots <- vector("list", p + 1L)
  for (j in seq_len(p)) {
    knots[[j]] <- structure(
      make_knots(x[[j]], k = k, degree = degree),
      cyclic = FALSE
    )
  }
  knots[[p + 1L]] <- structure(
    make_knots(as.numeric(lags), k = k_lag, degree = degree),
    cyclic = FALSE
  )
  names(knots) <- c(nms, "lag")

  # Lag basis on 0:L (shared)
  C_lag <- bspline_basis(as.numeric(lags), knots[[p + 1L]], degree = degree)
  storage.mode(C_lag) <- "double"

  # Per-lag basis lists on complete cases only
  idx <- which(ok)
  n_ok <- length(idx)
  basis_lags <- vector("list", length(lags))
  for (li in seq_along(lags)) {
    Bj <- vector("list", p + 1L)
    for (j in seq_len(p)) {
      Bj[[j]] <- bspline_basis(Q[[j]][idx, li], knots[[j]], degree = degree)
      storage.mode(Bj[[j]]) <- "double"
    }
    # Lag margin: each row is C_lag[li, ]
    Bj[[p + 1L]] <- matrix(C_lag[li, ], nrow = n_ok, ncol = ncol(C_lag),
                           byrow = TRUE)
    storage.mode(Bj[[p + 1L]]) <- "double"
    names(Bj) <- c(nms, "lag")
    basis_lags[[li]] <- Bj
  }

  list(
    names = nms,
    lag = lag,
    lags = lags,
    Q = Q,
    ok = ok,
    knots = knots,
    degree = as.integer(degree),
    k = as.integer(k),
    k_lag = as.integer(k_lag),
    C_lag = C_lag,
    basis_lags = basis_lags,
    n_ok = n_ok
  )
}

#' Sum of TT contractions over lags: f_t = sum_ℓ h(x_{t-ℓ}, ℓ).
#' @keywords internal
#' @noRd
tt_dlnm_contraction <- function(cores, basis_lags) {
  f <- 0
  for (ell in seq_along(basis_lags)) {
    f <- f + tt_contraction(cores, basis_lags[[ell]])
  }
  as.numeric(f)
}

#' Summed conditional design for core k across lags.
#' @keywords internal
#' @noRd
tt_dlnm_design_core <- function(cores, basis_lags, k) {
  X <- NULL
  for (ell in seq_along(basis_lags)) {
    b <- basis_lags[[ell]]
    L <- left_interfaces(cores, b)
    R <- right_interfaces(cores, b)
    Xe <- tt_design_core(L[[k]], R[[k]], b[[k]])
    if (is.null(X)) X <- Xe else X <- X + Xe
  }
  X
}

#' Per-margin core penalties with possibly unequal basis sizes.
#' @keywords internal
#' @noRd
tt_core_penalties_pvec <- function(ranks, p_vec, penalty_order = 2L,
                                   cyclic = NULL) {
  d <- length(ranks) - 1L
  p_vec <- as.integer(p_vec)
  if (length(p_vec) != d) {
    stop("p_vec must have length d = ", d, ".", call. = FALSE)
  }
  cyclic <- normalize_cyclic(cyclic, d)
  lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p_vec[k], ranks[k + 1L], penalty_order,
                 cyclic = cyclic[k])
  })
}

# ---- ALS / PIRLS ------------------------------------------------------------

#' One ALS sweep for TT-DLNM (summed lag designs), fixed or cGCV λ.
#' @keywords internal
#' @noRd
.tt_dlnm_als_sweep <- function(yc, cores, basis_lags, lambda, method, ranks,
                               control, penalties, p_vec, weight = NULL) {
  d <- length(cores)
  w <- if (is.null(weight)) rep(1, length(yc)) else as.numeric(weight)
  n_eval <- 0L
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  tol <- control$tol_lambda %||% 1e-3
  use_spec <- isTRUE(control$use_spectral_gcv %||% TRUE)
  for (k in seq_len(d)) {
    X <- tt_dlnm_design_core(cores, basis_lags, k)
    Pk <- penalties[[k]]
    if (identical(method, "cGCV")) {
      ws <- make_core_workspace(
        zc = yc, X = X, P = Pk, lambda0 = lambda[k],
        bounds = bounds, tol = tol, weight = w,
        use_spectral = use_spec, P0 = NULL
      )
      upd <- update_lambda_cgcv(ws)
      n_eval <- n_eval + upd$n_eval
      lambda[k] <- as.numeric(upd$lambda)
      cores[[k]] <- array(upd$g, dim = dim(cores[[k]]))
    } else {
      sw <- sqrt(pmax(w, 0))
      Xw <- X * sw
      yw <- yc * sw
      Gram <- crossprod(Xw)
      b <- as.numeric(crossprod(Xw, yw))
      M <- Gram + as.numeric(lambda[k]) * Pk
      g <- solve_spd_ridge(M, b)
      cores[[k]] <- array(g, dim = dim(cores[[k]]))
    }
  }
  list(cores = cores, lambda = lambda, n_eval = n_eval, penalties = penalties)
}

#' Gaussian ALS for TT-DLNM.
#' @keywords internal
#' @noRd
tt_dlnm_als_fit <- function(y, basis_lags, ranks, lambda_spec, control,
                            penalty_order = 2L, p_vec,
                            init_cores = NULL, offset = NULL, weights = NULL,
                            linear = NULL, smooth = NULL) {
  method <- lambda_spec$method
  lambda <- as.numeric(lambda_spec$values %||% lambda_spec$lambda0)
  n <- length(y)
  offset <- normalize_offset(offset, n)
  w <- normalize_weights(weights, n)
  linear <- normalize_linear(linear, n)
  smooth <- normalize_smooth(smooth, n)
  d <- length(ranks) - 1L
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p_vec, ranks, seed = control$seed,
                                 sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- tt_core_penalties_pvec(ranks, p_vec, penalty_order)

  ab0 <- tt_update_intercept_beta(y, offset, f = 0, linear = linear, weights = w)
  intercept <- ab0$intercept
  beta <- ab0$beta
  if (!is.null(smooth)) {
    add <- tt_refresh_additive(
      y, offset, f_tt = 0, linear = linear, smooth = smooth,
      weights = w, control = control
    )
    intercept <- add$intercept
    beta <- add$beta
    smooth <- add$smooth
  }

  n_eval <- 0L
  n_sweeps <- 0L
  for (sw in seq_len(control$max_sweeps %||% 20L)) {
    n_sweeps <- sw
    f <- tt_dlnm_contraction(cores, basis_lags)
    add <- tt_refresh_additive(
      y, offset, f_tt = f, linear = linear, smooth = smooth,
      weights = w, control = control
    )
    intercept <- add$intercept
    beta <- add$beta
    smooth <- add$smooth
    yc <- y - offset - intercept - tt_linear_contrib(linear, beta) -
      tt_smooth_contrib(smooth)
    step <- .tt_dlnm_als_sweep(
      yc, cores, basis_lags, lambda, method, ranks, control, penalties, p_vec,
      weight = w
    )
    cores <- step$cores
    lambda <- step$lambda
    n_eval <- n_eval + step$n_eval
  }
  f <- tt_dlnm_contraction(cores, basis_lags)
  mu <- offset + intercept + tt_linear_contrib(linear, beta) +
    tt_smooth_contrib(smooth) + f
  list(
    cores = cores,
    intercept = intercept,
    beta = beta,
    linear = linear,
    smooth = smooth,
    lambda = lambda,
    method_lambda = method,
    deviance = sum(w * (y - mu)^2),
    n_sweeps = n_sweeps,
    n_eval = n_eval,
    converged = TRUE,
    optimizer = "DLNM-ALS",
    edf = NA_real_
  )
}

#' PIRLS-ALS for Poisson/Bernoulli TT-DLNM.
#' @keywords internal
#' @noRd
tt_dlnm_pirls_fit <- function(y, basis_lags, family, ranks, lambda_spec, control,
                              penalty_order = 2L, p_vec,
                              init_cores = NULL, offset = NULL, weights = NULL,
                              linear = NULL, smooth = NULL) {
  method <- lambda_spec$method
  lambda <- as.numeric(lambda_spec$values %||% lambda_spec$lambda0)
  n <- length(y)
  fam <- normalize_family(family)
  key <- family_key(fam)
  offset <- normalize_offset(offset, n)
  w_obs <- normalize_weights(weights, n)
  linear <- normalize_linear(linear, n)
  smooth <- normalize_smooth(smooth, n)
  d <- length(ranks) - 1L

  # Link-scale intercept (NOT Gaussian OLS on counts)
  intercept <- init_intercept(fam, y, offset = offset, weights = w_obs)
  beta <- if (is.null(linear)) {
    numeric(0)
  } else {
    b <- rep(0, ncol(linear))
    names(b) <- colnames(linear)
    b
  }

  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p_vec, ranks, seed = control$seed,
                                 sd = control$init_sd)
    for (k in seq_len(d)) cores[[k]] <- cores[[k]] * 0.05
  } else {
    cores <- init_cores
  }
  penalties <- tt_core_penalties_pvec(ranks, p_vec, penalty_order)

  if (!is.null(smooth)) {
    # Soft start: smooths against working residual on link scale
    eta0 <- as.numeric(
      offset + intercept + tt_linear_contrib(linear, beta) +
        tt_dlnm_contraction(cores, basis_lags)
    )
    work0 <- glm_working(fam, y, eta0, control = control)
    add <- tt_refresh_additive(
      work0$z, offset, f_tt = tt_dlnm_contraction(cores, basis_lags),
      linear = linear, smooth = smooth,
      weights = work0$weight * w_obs, control = control
    )
    intercept <- add$intercept
    beta <- add$beta
    smooth <- add$smooth
  }

  eta <- as.numeric(
    offset + intercept + tt_linear_contrib(linear, beta) +
      tt_smooth_contrib(smooth) + tt_dlnm_contraction(cores, basis_lags)
  )
  n_eval <- 0L
  n_iter <- 0L
  maxit <- as.integer(control$pirls_maxit %||% 25L)
  als_sweeps <- as.integer(control$als_sweeps_per_pirls %||% 1L)
  tol <- as.numeric(control$tol %||% 1e-6)
  delta <- Inf

  for (it in seq_len(maxit)) {
    n_iter <- it
    work <- glm_working(fam, y, eta, control = control)
    ww <- work$weight * w_obs
    z <- work$z

    for (sw in seq_len(max(1L, als_sweeps))) {
      f <- tt_dlnm_contraction(cores, basis_lags)
      add <- tt_refresh_additive(
        z, offset, f_tt = f, linear = linear, smooth = smooth,
        weights = ww, control = control
      )
      intercept <- add$intercept
      beta <- add$beta
      smooth <- add$smooth
      yc <- z - offset - intercept - tt_linear_contrib(linear, beta) -
        tt_smooth_contrib(smooth)
      step <- .tt_dlnm_als_sweep(
        yc, cores, basis_lags, lambda, method, ranks, control, penalties,
        p_vec, weight = ww
      )
      cores <- step$cores
      lambda <- step$lambda
      n_eval <- n_eval + step$n_eval
    }
    eta_new <- as.numeric(
      offset + intercept + tt_linear_contrib(linear, beta) +
        tt_smooth_contrib(smooth) + tt_dlnm_contraction(cores, basis_lags)
    )
    # Keep eta in a sane range for Poisson/Bernoulli
    if (identical(key, "poisson")) {
      eta_new <- pmin(pmax(eta_new, -20), 20)
    }
    delta <- max(abs(eta_new - eta)) / (1 + max(abs(eta)))
    eta <- eta_new
    if (is.finite(delta) && delta < tol) break
  }
  mu <- invlink_eta(fam, eta)
  dev <- glm_deviance(fam, y, mu, weights = w_obs)
  list(
    cores = cores,
    intercept = intercept,
    beta = beta,
    linear = linear,
    smooth = smooth,
    lambda = lambda,
    method_lambda = method,
    deviance = dev,
    n_sweeps = n_iter * als_sweeps,
    n_iter = n_iter,
    n_pirls = n_iter,
    n_eval = n_eval,
    converged = is.finite(delta) && delta < tol,
    optimizer = "DLNM-PIRLS-ALS",
    edf = NA_real_
  )
}

# ---- Prediction (Gasparrini-style) ------------------------------------------

#' Evaluate TT-DLNM surface h at exposure grid × lag basis row.
#' @keywords internal
#' @noRd
.tt_dlnm_eval_h <- function(fit, Xgrid, lag_basis_row) {
  # Xgrid: m × p_exp; lag_basis_row: length k_lag (or m × k_lag)
  dlnm <- fit$dlnm
  p <- fit$p_exp
  m <- nrow(Xgrid)
  Bj <- vector("list", p + 1L)
  for (j in seq_len(p)) {
    Bj[[j]] <- bspline_basis(Xgrid[, j], dlnm$knots[[j]], degree = dlnm$degree)
  }
  if (is.null(dim(lag_basis_row))) {
    Bj[[p + 1L]] <- matrix(as.numeric(lag_basis_row), m, length(lag_basis_row),
                           byrow = TRUE)
  } else {
    Bj[[p + 1L]] <- as.matrix(lag_basis_row)
  }
  as.numeric(tt_contraction(fit$cores, Bj))
}

#' Predict overall / lag-slice cumulative effects from a [ttps_dlnm()] fit.
#'
#' @param fit A `"ttps_dlnm"` object.
#' @param var Name of exposure to vary (must be in `fit$dlnm$names`).
#' @param at Grid of values for `var`.
#' @param cen Named vector / list of centering values for all exposures
#'   (defaults to training medians of unlagged series). `cen[[var]]` is the
#'   reference for RR.
#' @param type `"overall"` (sum over lags; Gasparrini overall) or `"slice"`
#'   (single lag).
#' @param lag_at Lag value for `type = "slice"` (in `0:L`).
#' @param se Currently ignored (Level-1 SE for DLNM not yet wired).
#' @return A data.frame with `x`, `logRR`, `RR` (and `lag` for slices).
#' @export
predict_dlnm <- function(fit, var, at, cen = NULL,
                         type = c("overall", "slice"),
                         lag_at = 0L,
                         se = FALSE) {
  stopifnot(inherits(fit, "ttps_dlnm"))
  type <- match.arg(type)
  dlnm <- fit$dlnm
  nms <- dlnm$names
  if (!var %in% nms) {
    stop("`var` must be one of: ", paste(nms, collapse = ", "), call. = FALSE)
  }
  at <- as.numeric(at)
  # centering
  if (is.null(cen)) {
    cen <- vapply(nms, function(nm) {
      median(dlnm$Q[[nm]][dlnm$ok, 1L], na.rm = TRUE)
    }, numeric(1))
  } else {
    cen <- unlist(cen)
    if (is.null(names(cen))) {
      stop("`cen` must be a named vector matching exposure names.", call. = FALSE)
    }
    miss <- setdiff(nms, names(cen))
    if (length(miss)) {
      for (nm in miss) {
        cen[nm] <- median(dlnm$Q[[nm]][dlnm$ok, 1L], na.rm = TRUE)
      }
    }
  }
  m <- length(at)
  Xg <- matrix(as.numeric(cen[nms]), m, length(nms), byrow = TRUE)
  colnames(Xg) <- nms
  Xg[, var] <- at
  X0 <- matrix(as.numeric(cen[nms]), 1L, length(nms), byrow = TRUE)
  colnames(X0) <- nms

  if (identical(type, "overall")) {
    # sum_ℓ C(ℓ)  ->  length-k_lag vector
    clag <- as.numeric(colSums(dlnm$C_lag))
    h <- .tt_dlnm_eval_h(fit, Xg, clag)
    h0 <- .tt_dlnm_eval_h(fit, X0, clag)
    logRR <- as.numeric(h - h0)
    data.frame(
      var = var, x = at, cen = unname(cen[[var]]),
      type = "overall",
      logRR = logRR, RR = exp(logRR),
      row.names = NULL
    )
  } else {
    lags <- dlnm$lags
    lag_at <- as.integer(lag_at)
    if (!lag_at %in% lags) {
      stop("`lag_at` must be in 0:L (L = ", dlnm$lag, ").", call. = FALSE)
    }
    li <- match(lag_at, lags)
    clag <- as.numeric(dlnm$C_lag[li, ])
    h <- .tt_dlnm_eval_h(fit, Xg, clag)
    h0 <- .tt_dlnm_eval_h(fit, X0, clag)
    logRR <- as.numeric(h - h0)
    data.frame(
      var = var, x = at, cen = unname(cen[[var]]),
      type = "slice", lag = lag_at,
      logRR = logRR, RR = exp(logRR),
      row.names = NULL
    )
  }
}

#' @export
print.ttps_dlnm <- function(x, ...) {
  cat("Tensor-Train DLNM fit\n")
  cat(sprintf("  Exposures: %s\n", paste(x$dlnm$names, collapse = ", ")))
  cat(sprintf("  Lag: 0:%d  (k=%d, k_lag=%d, rank=%s)\n",
              x$dlnm$lag, x$k, x$k_lag, paste(x$rank, collapse = "-")))
  cat(sprintf("  Family: %s   deviance: %.4g\n",
              x$family$family, x$deviance))
  cat(sprintf("  Lambda (%s): %s\n",
              x$lambda_method,
              paste(sprintf("%.3g", x$lambda), collapse = ", ")))
  if (!is.null(x$smooth)) {
    cat(sprintf("  Smooths: %s\n", paste(names(x$smooth), collapse = ", ")))
  }
  invisible(x)
}
