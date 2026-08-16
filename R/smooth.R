# Additive 1D P-spline / cyclic smooths (mgcv-like bs = "ps" / "cc").
# Estimated jointly with TT via backfitting; lambda_sm via conditional cGCV.

#' Contribution of additive smooths: sum_j B_j gamma_j.
#' @keywords internal
#' @noRd
tt_smooth_contrib <- function(smooth) {
  if (is.null(smooth) || length(smooth) == 0L) return(0)
  out <- 0
  for (sm in smooth) {
    if (is.null(sm$B) || is.null(sm$gamma) || length(sm$gamma) == 0L) next
    out <- out + as.numeric(sm$B %*% sm$gamma)
  }
  out
}

#' Penalty value 0.5 sum_j lambda_j gamma_j' S_j gamma_j.
#' @keywords internal
#' @noRd
tt_smooth_penalty_value <- function(smooth) {
  if (is.null(smooth) || length(smooth) == 0L) return(0)
  pen <- 0
  for (sm in smooth) {
    if (is.null(sm$gamma) || is.null(sm$S) || !is.finite(sm$lambda)) next
    g <- as.numeric(sm$gamma)
    pen <- pen + 0.5 * as.numeric(sm$lambda) * as.numeric(crossprod(g, sm$S %*% g))
  }
  pen
}

#' Build one additive smooth term (basis B, penalty S, init gamma).
#'
#' @param x Covariate vector (length n).
#' @param bs `"ps"` (open P-spline) or `"cc"` (cyclic / circular).
#' @param k Basis size.
#' @param m Difference penalty order (alias of `penalty_order`).
#' @param degree B-spline degree (default 3).
#' @param period For `bs = "cc"`: `c(xl, xr)`; default from data / unit interval.
#' @param lambda Numeric (fixed), `"cGCV"`, or ignored when `target_edf` is set.
#' @param target_edf Optional target effective df: choose \(\lambda\) so that
#'   \(\mathrm{edf}(\lambda)\approx\) `target_edf` (monotone root on log-\(\lambda\)).
#' @param name Term label for summary / predict.
#' @keywords internal
#' @noRd
build_smooth_term <- function(x,
                              bs = c("ps", "cc"),
                              k = 10L,
                              m = 2L,
                              degree = 3L,
                              period = NULL,
                              lambda = "cGCV",
                              target_edf = NULL,
                              name = "s") {
  bs <- match.arg(bs)
  x <- as.numeric(x)
  if (anyNA(x)) stop("smooth covariate '", name, "' contains NA.", call. = FALSE)
  n <- length(x)
  k <- as.integer(k)
  m <- as.integer(m)
  degree <- as.integer(degree)
  if (k <= degree) {
    stop("smooth '", name, "': need k > degree (", degree, ").", call. = FALSE)
  }
  if (m < 1L || m >= k) {
    stop("smooth '", name, "': penalty order m must satisfy 1 <= m < k.",
         call. = FALSE)
  }

  if (identical(bs, "cc")) {
    pr <- cyclic_period_range(x, period)
    knots <- make_cyclic_knots(pr[1], pr[2], k = k, degree = degree)
    B <- cyclic_bspline_basis(x, knots = knots, degree = degree, k = k,
                              xl = pr[1], xr = pr[2])
    S <- circular_difference_penalty(k, m)
    period_out <- pr
  } else {
    knots <- structure(make_knots(x, k = k, degree = degree), cyclic = FALSE)
    B <- bspline_basis(x, knots, degree = degree)
    if (ncol(B) != k) {
      stop("smooth '", name, "': open basis ncol=", ncol(B), " != k=", k, ".",
           call. = FALSE)
    }
    S <- difference_penalty(k, m)
    period_out <- NULL
  }
  storage.mode(B) <- "double"
  storage.mode(S) <- "double"

  tgt <- target_edf
  if (!is.null(tgt)) {
    tgt <- as.numeric(tgt)
    if (length(tgt) != 1L || !is.finite(tgt)) {
      stop("smooth '", name, "': target_edf must be a finite scalar.",
           call. = FALSE)
    }
    # Difference penalty null space has dimension m; edf is roughly in (m, k]
    if (tgt <= m || tgt > k) {
      stop(
        "smooth '", name, "': target_edf must satisfy m < target_edf <= k ",
        "(got target_edf=", tgt, ", m=", m, ", k=", k, ").",
        call. = FALSE
      )
    }
    if (is.character(lambda) && identical(lambda, "cGCV")) {
      # target_edf overrides automatic cGCV
    } else if (is.numeric(lambda)) {
      warning(
        "smooth '", name, "': target_edf set; ignoring numeric lambda=",
        lambda[1], " for selection (still used only as a start if needed).",
        call. = FALSE
      )
    }
    lam_method <- "target_edf"
    lam0 <- 1
  } else if (is.character(lambda) && identical(lambda, "cGCV")) {
    lam_method <- "cGCV"
    lam0 <- 1
    tgt <- NA_real_
  } else {
    lam_method <- "fixed"
    lam0 <- as.numeric(lambda)
    if (length(lam0) != 1L || !is.finite(lam0) || lam0 < 0) {
      stop(
        "smooth '", name,
        "': lambda must be a nonnegative scalar, \"cGCV\", or use target_edf.",
        call. = FALSE
      )
    }
    tgt <- NA_real_
  }

  list(
    name = as.character(name),
    bs = bs,
    k = k,
    m = m,
    degree = degree,
    x = x,
    B = B,
    S = S,
    knots = knots,
    period = period_out,
    lambda = lam0,
    lambda_method = lam_method,
    target_edf = tgt,
    gamma = rep(0, k),
    edf = NA_real_
  )
}

#' Normalize `smooth=` argument into a list of built terms.
#'
#' Accepts:
#' - `NULL`
#' - named list: `list(time = list(x = ..., bs = "ps", k = 20, m = 2), ...)`
#' - unnamed list of the same element form (names from `name=` or `s1`, `s2`, ...)
#'
#' @keywords internal
#' @noRd
normalize_smooth <- function(smooth, n, lambda_smooth = "cGCV") {
  n <- as.integer(n)
  if (is.null(smooth)) return(NULL)
  if (!is.list(smooth)) {
    stop("`smooth` must be NULL or a list of smooth term specs.", call. = FALSE)
  }
  if (length(smooth) == 0L) return(NULL)

  # Already built (idempotent)
  if (isTRUE(attr(smooth, "ttps_smooth_ok"))) {
    for (sm in smooth) {
      if (is.null(sm$B) || nrow(sm$B) != n) {
        stop("built `smooth` terms must have B with ", n, " rows.", call. = FALSE)
      }
    }
    return(smooth)
  }

  nms <- names(smooth)
  out <- vector("list", length(smooth))
  for (i in seq_along(smooth)) {
    spec <- smooth[[i]]
    if (is.null(spec) || !is.list(spec)) {
      stop("`smooth[[", i, "]]` must be a list with at least `x`.", call. = FALSE)
    }
    if (is.null(spec$x)) {
      stop("`smooth[[", i, "]]` needs element `x`.", call. = FALSE)
    }
    nm <- if (!is.null(nms) && nzchar(nms[i])) {
      nms[i]
    } else if (!is.null(spec$name)) {
      as.character(spec$name)
    } else {
      paste0("s", i)
    }
    if (length(spec$x) != n) {
      stop("smooth '", nm, "': length(x) must equal n = ", n, ".", call. = FALSE)
    }
    m_use <- spec$m %||% spec$penalty_order %||% 2L
    lam_use <- spec$lambda %||% lambda_smooth
    out[[i]] <- build_smooth_term(
      x = spec$x,
      bs = spec$bs %||% "ps",
      k = spec$k %||% 10L,
      m = m_use,
      degree = spec$degree %||% 3L,
      period = spec$period,
      lambda = lam_use,
      target_edf = spec$target_edf,
      name = nm
    )
  }
  names(out) <- vapply(out, `[[`, character(1), "name")
  if (anyDuplicated(names(out))) {
    stop("`smooth` term names must be unique.", call. = FALSE)
  }
  attr(out, "ttps_smooth_ok") <- TRUE
  out
}

#' Subset built smooth terms to observation index (CV folds).
#' Keeps knots / S / lambda from the full-data build.
#' @keywords internal
#' @noRd
subset_smooth <- function(smooth, idx) {
  if (is.null(smooth)) return(NULL)
  idx <- as.integer(idx)
  out <- lapply(smooth, function(sm) {
    sm$B <- sm$B[idx, , drop = FALSE]
    sm$x <- sm$x[idx]
    sm$gamma <- sm$gamma %||% rep(0, sm$k)
    sm
  })
  names(out) <- names(smooth)
  attr(out, "ttps_smooth_ok") <- TRUE
  out
}

#' Evaluate fitted smooths at new covariate values.
#'
#' `newdata_smooth` is a named list of numeric vectors (or `list(x = ...)`)
#' with the same names as `object$smooth`.
#' @keywords internal
#' @noRd
eval_smooth_newdata <- function(smooth_fit, newdata_smooth) {
  if (is.null(smooth_fit) || length(smooth_fit) == 0L) return(0)
  if (is.null(newdata_smooth)) {
    stop(
      "This fit has additive `smooth` terms; pass matching `smooth=` ",
      "for newdata (named list of covariate vectors).",
      call. = FALSE
    )
  }
  if (!is.list(newdata_smooth)) {
    stop("`smooth` for predict must be a named list.", call. = FALSE)
  }
  nms <- names(smooth_fit)
  if (is.null(names(newdata_smooth)) ||
      !all(nms %in% names(newdata_smooth))) {
    stop(
      "`smooth` for predict must include names: ",
      paste(nms, collapse = ", "),
      call. = FALSE
    )
  }
  n <- NULL
  contrib <- 0
  for (nm in nms) {
    sm <- smooth_fit[[nm]]
    raw <- newdata_smooth[[nm]]
    xnew <- if (is.list(raw) && !is.null(raw$x)) {
      as.numeric(raw$x)
    } else {
      as.numeric(raw)
    }
    if (is.null(n)) n <- length(xnew)
    if (length(xnew) != n) {
      stop("smooth '", nm, "' newdata length mismatch.", call. = FALSE)
    }
    if (identical(sm$bs, "cc")) {
      pr <- sm$period
      Bnew <- cyclic_bspline_basis(
        xnew, knots = sm$knots, degree = sm$degree, k = sm$k,
        xl = pr[1], xr = pr[2]
      )
    } else {
      Bnew <- bspline_basis(xnew, sm$knots, degree = sm$degree)
    }
    contrib <- contrib + as.numeric(Bnew %*% sm$gamma)
  }
  contrib
}

#' Choose \(\lambda\) so \(\mathrm{edf}(\lambda)\approx\) `target` (monotone).
#'
#' Uses `uniroot` on \(\log\lambda\); expands the upper bound if needed.
#' As \(\lambda\to 0\), edf \(\to\) nearly \(k\); as \(\lambda\to\infty\),
#' edf \(\to m\) (difference-penalty null space).
#'
#' @keywords internal
#' @noRd
.lambda_for_target_edf <- function(Gram, S, target, bounds = c(1e-4, 1e4),
                                   tol = 1e-4) {
  target <- as.numeric(target)
  lo <- max(as.numeric(bounds[[1]]), 1e-12)
  hi <- max(as.numeric(bounds[[2]]), lo * 10)
  ed_at <- function(lam) .ed_S(Gram, S, lam, P0 = NULL)

  ed_lo <- ed_at(lo)
  if (target >= ed_lo - 1e-6) {
    return(list(lambda = lo, edf = ed_lo))
  }

  ed_hi <- ed_at(hi)
  expand <- 0L
  while (ed_hi > target + 1e-6 && hi < 1e14 && expand < 24L) {
    hi <- hi * 10
    ed_hi <- ed_at(hi)
    expand <- expand + 1L
  }
  if (ed_hi > target + 1e-6) {
    return(list(lambda = hi, edf = ed_hi))
  }

  root <- tryCatch(
    stats::uniroot(
      function(ll) ed_at(exp(ll)) - target,
      interval = log(c(lo, hi)),
      tol = tol
    ),
    error = function(e) NULL
  )
  if (is.null(root)) {
    ll_grid <- seq(log(lo), log(hi), length.out = 81L)
    eds <- vapply(ll_grid, function(ll) ed_at(exp(ll)), numeric(1))
    j <- which.min(abs(eds - target))
    lam <- exp(ll_grid[[j]])
  } else {
    lam <- exp(root$root)
  }
  list(lambda = as.numeric(lam), edf = ed_at(lam))
}

#' Update one smooth term (penalized WLS + cGCV / target_edf / fixed).
#' @keywords internal
#' @noRd
tt_update_smooth_term <- function(sm, residual, weights, control) {
  w <- as.numeric(weights)
  r <- as.numeric(residual)
  B <- sm$B
  S <- sm$S
  bounds <- control$lambda_bounds %||% c(1e-4, 1e4)
  tol <- control$tol_lambda %||% 1e-3
  use_spec <- isTRUE(control$use_spectral_gcv %||% TRUE)

  sw <- sqrt(pmax(w, 0))
  Bw <- B * sw
  yw <- r * sw
  Gram <- crossprod(Bw)
  b <- as.numeric(crossprod(Bw, yw))

  if (identical(sm$lambda_method, "cGCV")) {
    ws <- make_core_workspace(
      zc = r, X = B, P = S, lambda0 = sm$lambda,
      bounds = bounds, tol = tol, weight = w,
      use_spectral = use_spec, P0 = NULL
    )
    upd <- update_lambda_cgcv(ws)
    sm$lambda <- as.numeric(upd$lambda)
    sm$gamma <- as.numeric(upd$g)
    sm$edf <- as.numeric(upd$ed)
  } else if (identical(sm$lambda_method, "target_edf")) {
    hit <- .lambda_for_target_edf(
      Gram, S, target = sm$target_edf, bounds = bounds, tol = max(tol, 1e-4)
    )
    sm$lambda <- hit$lambda
    M <- Gram + as.numeric(sm$lambda) * S
    sm$gamma <- as.numeric(solve_spd_ridge(M, b))
    sm$edf <- hit$edf
  } else {
    M <- Gram + as.numeric(sm$lambda) * S
    sm$gamma <- as.numeric(solve_spd_ridge(M, b))
    sm$edf <- .ed_S(Gram, S, sm$lambda, P0 = NULL)
  }
  names(sm$gamma) <- paste0(sm$name, ".", seq_along(sm$gamma))
  sm
}

#' Backfit all additive smooths given fixed TT + parametric parts.
#'
#' `f_fixed` = TT + intercept is NOT included; pass
#' `base = offset + intercept + linear*beta + f_tt` and we form
#' residuals against other smooths internally.
#'
#' @param target Working response (y or PIRLS z).
#' @param base offset + intercept + linear contrib + TT (no smooths).
#' @keywords internal
#' @noRd
tt_update_smooths <- function(target, base, smooth, weights, control,
                              passes = 1L) {
  if (is.null(smooth) || length(smooth) == 0L) return(smooth)
  target <- as.numeric(target)
  base <- as.numeric(base)
  w <- normalize_weights(weights, length(target))
  passes <- max(1L, as.integer(passes))
  for (pass in seq_len(passes)) {
    for (j in seq_along(smooth)) {
      others <- 0
      if (length(smooth) > 1L) {
        for (ell in seq_along(smooth)) {
          if (ell == j) next
          others <- others + as.numeric(smooth[[ell]]$B %*% smooth[[ell]]$gamma)
        }
      }
      resid_j <- target - base - others
      smooth[[j]] <- tt_update_smooth_term(smooth[[j]], resid_j, w, control)
    }
  }
  attr(smooth, "ttps_smooth_ok") <- TRUE
  smooth
}

#' Refresh intercept/beta and smooths given current TT surface.
#' @keywords internal
#' @noRd
tt_refresh_additive <- function(target, offset, f_tt, linear, smooth,
                                weights, control, smooth_passes = 1L) {
  offset <- normalize_offset(offset, length(target))
  f_tt <- as.numeric(f_tt)
  # Smooths first (with current intercept absorbed into residual via
  # provisional intercept update after), then parametric.
  # Start: update smooths treating current parametric as fixed via base
  # that excludes smooth; use a provisional intercept from residual mean.
  ab0 <- tt_update_intercept_beta(
    target, offset, f_tt + tt_smooth_contrib(smooth),
    linear = linear, weights = weights
  )
  intercept <- ab0$intercept
  beta <- ab0$beta
  base <- offset + intercept + tt_linear_contrib(linear, beta) + f_tt
  smooth <- tt_update_smooths(
    target, base, smooth, weights, control, passes = smooth_passes
  )
  ab <- tt_update_intercept_beta(
    target, offset, f_tt + tt_smooth_contrib(smooth),
    linear = linear, weights = weights
  )
  list(
    intercept = ab$intercept,
    beta = ab$beta,
    smooth = smooth
  )
}

#' Compact smooth summary table (mgcv-like).
#' @keywords internal
#' @noRd
.tt_smooth_summary_table <- function(smooth) {
  if (is.null(smooth) || length(smooth) == 0L) return(NULL)
  rows <- lapply(smooth, function(sm) {
    data.frame(
      term = paste0("s(", sm$name, ")"),
      bs = sm$bs,
      k = sm$k,
      m = sm$m,
      edf = if (is.finite(sm$edf)) sm$edf else NA_real_,
      target_edf = if (is.finite(sm$target_edf %||% NA_real_)) {
        as.numeric(sm$target_edf)
      } else {
        NA_real_
      },
      lambda = sm$lambda,
      method = sm$lambda_method,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- tab$term
  tab$term <- NULL
  tab
}
