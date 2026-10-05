# Exact effective degrees of freedom of the fixed-lambda TT estimator.
#
# The penalized objective minimised by ALS / PIRLS is
#   Q(theta; y) = sum_i w_i l(y_i, eta_i) + pen(theta),
# l = (y - eta)^2 / 2 (Gaussian) or mu - y log mu (Poisson, log link), pen the
# global penalty. At a stationary point theta* the implicit function theorem
# gives d theta* / d y = H^+ J' W0, with
#   H  = d^2 Q / d theta^2 = J'WJ + sum_i c_i d^2 eta_i / d theta^2 + d^2 pen,
#   J  = d eta / d theta (stacked core designs plus an intercept column),
#   c_i = w_i (mu_i - y_i),  W = diag(w_i dmu_i/deta_i),  W0 = diag(w_i).
# Hence EDF = tr(d mu / d y) = tr(H^+ J'WJ).
#
# H is assembled analytically. eta is linear in each core, so the curvature
# term has only cross-core blocks; block (k, l), k < l, is
#   sum_i c_i [L_k(i) (x) B_k(i)] M_kl(i) [B_l(i) (x) R_l(i)],
# M_kl(i) the product of the evaluated cores strictly between k and l. The
# penalty is quadratic in each core, so central differences of its gradient
# are exact up to rounding (and cost no data pass).
#
# TT gauge directions (G_n A, -A G_{n+1} on every bond) leave Q unchanged.
# They are projected out analytically. An eigenvalue cutoff is not enough:
# away from an exact stationary point their curvature is of order
# |grad| |theta| and can be negative.
#
# The trace is the derivative of the converged map. ALS stopped at a
# tolerance is not stationary at low rank (6 x 5 x 4 Gaussian, rank 2,
# lambda = 0.03: |grad| 1e-2 after 3000 sweeps, EDF 21.09 against 22.55 at
# the stationary point), so the cores are first polished by a gauge-projected
# trust-region Newton iteration. Checks against finite-difference
# traces of the converged map: 22.5504 / 22.554 (lambda 0.03) and
# 12.4978 / 12.4978 (lambda 1); saturating rank reproduces the dense EDF.

#' Gradient, Hessian and J'WJ of the penalized TT objective (internal).
#' @keywords internal
#' @noRd
.tt_exact_parts <- function(cores, intercept, basis, y, key, w, off, lambda,
                            penalty_order, cyclic, hessian = TRUE) {
  d <- length(cores)
  n <- length(y)
  sz <- vapply(cores, length, integer(1))
  pos <- c(0L, cumsum(sz))
  np <- pos[d + 1L] + 1L
  eta <- off + intercept + tt_contraction(cores, basis)
  if (identical(key, "gaussian")) {
    mu <- eta
    dmu <- rep(1, n)
    nll <- 0.5 * sum(w * (y - eta)^2)
  } else {
    mu <- exp(eta)
    dmu <- mu
    nll <- sum(w * (mu - y * eta))
  }
  cvec <- w * (mu - y)
  pen <- .tt_penalty_value_grad_full(cores, lambda, penalty_order = penalty_order,
                                     cyclic = cyclic)
  value <- nll + pen$value
  L <- left_interfaces(cores, basis)
  R <- right_interfaces(cores, basis)
  if (!hessian) {
    g <- unlist(lapply(seq_len(d), function(k) {
      as.numeric(crossprod(tt_design_core(L[[k]], R[[k]], basis[[k]]), cvec))
    }), use.names = FALSE)
    grad <- c(g + .tt_pack_cores(pen$grads), sum(cvec))
    return(list(value = value, grad = grad, eta = eta, mu = mu))
  }
  J <- cbind(do.call(cbind, lapply(seq_len(d), function(k) {
    tt_design_core(L[[k]], R[[k]], basis[[k]])
  })), 1)
  grad <- as.numeric(crossprod(J, cvec))
  grad[-np] <- grad[-np] + .tt_pack_cores(pen$grads)
  G <- crossprod(J, (w * dmu) * J)
  H <- G
  # curvature of the multilinear map: cross-core blocks only
  if (d >= 2L) {
    for (k in seq_len(d - 1L)) {
      rkl <- dim(cores[[k]])[1L]
      Kk <- dim(cores[[k]])[2L]
      rk <- dim(cores[[k]])[3L]
      # U[, a + rkl (j - 1)] = L_k[, a] B_k[, j]
      U <- L[[k]][, rep(seq_len(rkl), Kk), drop = FALSE] *
        basis[[k]][, rep(seq_len(Kk), each = rkl), drop = FALSE]
      # M[[b]]: rows e_b' prod_{k < s < l} G_s[B_s(i)], advanced with l
      M <- lapply(seq_len(rk), function(b) {
        m <- matrix(0, n, rk)
        m[, b] <- 1
        m
      })
      for (l in (k + 1L):d) {
        if (l > k + 1L) {
          M <- lapply(M, function(m) {
            contract_left_step(m, cores[[l - 1L]], basis[[l - 1L]])
          })
        }
        rll <- dim(cores[[l]])[1L]
        Kl <- dim(cores[[l]])[2L]
        rl <- dim(cores[[l]])[3L]
        # V[, m + Kl (e - 1)] = B_l[, m] R_l[, e]
        V <- basis[[l]][, rep(seq_len(Kl), rl), drop = FALSE] *
          R[[l]][, rep(seq_len(rl), each = Kl), drop = FALSE]
        blk <- matrix(0, sz[k], sz[l])
        col_base <- (seq_len(Kl * rl) - 1L) * rll  # core-l index of (1, m, e) - 1
        for (b in seq_len(rk)) {
          rows <- (b - 1L) * rkl * Kk + seq_len(rkl * Kk)
          for (g in seq_len(rll)) {
            blk[rows, col_base + g] <- crossprod(U, (cvec * M[[b]][, g]) * V)
          }
        }
        ik <- pos[k] + seq_len(sz[k])
        il <- pos[l] + seq_len(sz[l])
        H[ik, il] <- H[ik, il] + blk
        H[il, ik] <- H[il, ik] + t(blk)
      }
    }
  }
  # penalty Hessian: central differences of its gradient (exact: the penalty
  # is quadratic in each core)
  th <- .tt_pack_cores(cores)
  bscale <- rep(vapply(cores, function(a) sqrt(mean(a^2)), numeric(1)), sz)
  pgrad <- function(t) {
    .tt_pack_cores(.tt_penalty_value_grad_full(
      .tt_unpack_cores(t, cores), lambda, penalty_order = penalty_order,
      cyclic = cyclic
    )$grads)
  }
  Hp <- matrix(0, np - 1L, np - 1L)
  for (j in seq_len(np - 1L)) {
    h <- max(bscale[j], 1e-8)
    e <- numeric(np - 1L)
    e[j] <- h
    Hp[, j] <- (pgrad(th + e) - pgrad(th - e)) / (2 * h)
  }
  H[-np, -np] <- H[-np, -np] + (Hp + t(Hp)) / 2
  list(value = value, grad = grad, H = (H + t(H)) / 2, G = G, eta = eta, mu = mu)
}

#' Orthonormal basis of the complement of the TT gauge directions (internal).
#' Columns span the packed parameter space (cores, then intercept) minus the
#' tangent of G_n -> G_n A, G_{n+1} -> A^{-1} G_{n+1} on every bond.
#' @keywords internal
#' @noRd
.tt_gauge_complement <- function(cores) {
  d <- length(cores)
  sz <- vapply(cores, length, integer(1))
  pos <- c(0L, cumsum(sz))
  np <- pos[d + 1L] + 1L
  cols <- list()
  if (d >= 2L) {
    for (nb in seq_len(d - 1L)) {
      Gn <- cores[[nb]]
      Gm <- cores[[nb + 1L]]
      r <- dim(Gn)[3L]
      for (a in seq_len(r)) for (b in seq_len(r)) {
        dn <- array(0, dim(Gn))
        dn[, , b] <- Gn[, , a]                 # G_n E_ab
        dm <- array(0, dim(Gm))
        dm[a, , ] <- -Gm[b, , ]                # -E_ab G_{n+1}
        v <- numeric(np)
        v[pos[nb] + seq_len(sz[nb])] <- as.vector(dn)
        v[pos[nb + 1L] + seq_len(sz[nb + 1L])] <- as.vector(dm)
        cols[[length(cols) + 1L]] <- v
      }
    }
  }
  if (!length(cols)) return(diag(np))
  qz <- qr(do.call(cbind, cols))
  Q <- qr.Q(qz, complete = TRUE)
  Q[, -seq_len(qz$rank), drop = FALSE]
}

#' Exact EDF of a TT fit, with trust-region Newton polishing (internal engine).
#'
#' @param cores,intercept Fitted TT cores and intercept (starting point).
#' @param basis List of marginal bases evaluated at the observations
#'   (attribute `"cyclic"` as from `eval_marginal_bases()`).
#' @param y Response; `fam` family object; `offset`, `weights` vectors.
#' @param lambda Smoothing vector; `penalty_order` difference order.
#' @param polish Polish the cores to a stationary point before the trace.
#' @param maxit Trust-region Newton iterations (polishing).
#' @param tol Stop when the Newton decrement is at most `tol * max(1, |Q|)`.
#' @param rel_cut Eigenvalue cutoff (relative to the largest |eigenvalue|) of
#'   the gauge-projected Hessian; only exactly flat directions fall below it.
#' @return List with `edf`, `converged` (Newton decrement below `tol` and
#'   positive curvature off the gauge), `n_hessian`,
#'   `grad_norm` (projected gradient at the returned point), `decrement`,
#'   `n_newton`, `n_negative` (negative curvature directions: a saddle),
#'   `n_flat`, `npar`, `cores`, `intercept`, `eta`, `mu`, `time_s`.
#' @keywords internal
#' @noRd
.tt_exact_edf_core <- function(cores, intercept, basis, y, fam, offset = NULL,
                               weights = NULL, lambda, penalty_order = 2L,
                               polish = TRUE, maxit = 12L, tol = 1e-12,
                               rel_cut = 1e-10, trace = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  y <- as.numeric(y)
  n <- length(y)
  key <- family_key(fam)
  if (!key %in% c("gaussian", "poisson")) {
    stop("Exact EDF supports gaussian() and poisson() only.", call. = FALSE)
  }
  link <- fam$link %||% if (identical(key, "poisson")) "log" else "identity"
  if (!identical(link, if (identical(key, "poisson")) "log" else "identity")) {
    stop("Exact EDF supports gaussian(identity) and poisson(log) only.",
         call. = FALSE)
  }
  off <- normalize_offset(offset, n)
  w <- normalize_weights(weights, n)
  lam <- as.numeric(lambda)
  cyclic <- attr(basis, "cyclic")
  # gauge balancing: equal Frobenius norms (function and penalty unchanged)
  nr <- vapply(cores, function(a) sqrt(sum(a^2)), numeric(1))
  if (any(!is.finite(nr)) || any(nr <= 0)) {
    stop("Exact EDF needs finite, non-zero cores.", call. = FALSE)
  }
  cores <- Map(function(a, cc) a * cc, cores, exp(mean(log(nr))) / nr)
  a0 <- as.numeric(intercept)
  np <- sum(vapply(cores, length, integer(1))) + 1L
  parts_at <- function(cr, a, hessian = TRUE) {
    .tt_exact_parts(cr, a, basis, y, key, w, off, lam, penalty_order, cyclic,
                    hessian = hessian)
  }
  project <- function(cr, p) {
    Q <- .tt_gauge_complement(cr)
    Hq <- crossprod(Q, p$H %*% Q)
    list(Q = Q, ev = eigen((Hq + t(Hq)) / 2, symmetric = TRUE),
         gq = as.numeric(crossprod(Q, p$grad)))
  }
  # Trust-region Newton (Levenberg) on the eigendecomposition of the
  # gauge-projected Hessian: s(m) = -sum_i u_i / (l_i + m) v_i with
  # m > -min(l_i), m adapted by the ratio of actual to predicted decrease.
  # Retrying a step costs one objective value. ALS can stop near a saddle
  # (Gaussian Dette d = 4, rank 4: negative curvature at the ALS point), where
  # a line-searched saddle-free step wanders along nearly flat directions;
  # near a minimum m -> 0 and the iteration is Newton (quadratic). Converged
  # = Newton decrement below tol at a point with positive curvature in every
  # non-gauge direction.
  p <- parts_at(cores, a0)
  pr <- project(cores, p)
  n_hess <- 1L
  it <- 0L
  dec <- NA_real_
  converged <- FALSE
  mu_tr <- NA_real_
  repeat {
    ev <- pr$ev
    l <- ev$values
    lmax <- max(abs(l))
    u <- as.numeric(crossprod(ev$vectors, pr$gq))
    # flat: exact null directions beyond the gauge (e.g. intercept against a
    # constant inside the tensor, which the penalty does not see)
    flat <- abs(l) <= rel_cut * lmax
    pos <- all(l[!flat] > 0)
    dec <- if (pos) sum(u[!flat]^2 / l[!flat]) else Inf
    if (pos && dec <= tol * max(1, abs(p$value))) {
      converged <- TRUE
      break
    }
    if (isTRUE(trace)) {
      message(sprintf("TR it %2d  Q %.10g  dec %.3e  |g| %.3e  min ev %.3e  n_neg %d  mu %.2e",
                      it, p$value, dec, sqrt(sum(pr$gq^2)), min(l), sum(l < 0),
                      mu_tr))
    }
    if (!polish || it >= maxit) break
    m_min <- max(0, -min(l)) + rel_cut * lmax
    if (!is.finite(mu_tr)) mu_tr <- if (pos) 0 else 1e-3 * lmax
    m <- m_min + mu_tr
    th <- c(.tt_pack_cores(cores), a0)
    ok <- FALSE
    for (try in seq_len(30L)) {
      z <- u / (l + m)
      pred <- sum(u * z) - 0.5 * sum(l * z^2)   # predicted decrease (> 0)
      th1 <- th - as.numeric(pr$Q %*% (ev$vectors %*% z))
      cr1 <- .tt_unpack_cores(th1[-np], cores)
      v1 <- parts_at(cr1, th1[np], hessian = FALSE)$value
      rho <- if (is.finite(v1) && pred > 0) (p$value - v1) / pred else -Inf
      if (rho > 1e-4) {
        ok <- TRUE
        break
      }
      mu_tr <- max(4 * mu_tr, 1e-6 * lmax)
      m <- m_min + mu_tr
    }
    if (!ok) break
    mu_tr <- if (rho > 0.75) mu_tr / 4 else if (rho < 0.25) 4 * mu_tr else mu_tr
    if (mu_tr < 1e-10 * lmax) mu_tr <- 0
    it <- it + 1L
    cores <- cr1
    a0 <- th1[np]
    p <- parts_at(cores, a0)
    pr <- project(cores, p)
    n_hess <- n_hess + 1L
  }
  ev <- pr$ev
  keep <- abs(ev$values) > rel_cut * max(abs(ev$values))
  V <- ev$vectors[, keep, drop = FALSE]
  Gq <- crossprod(pr$Q, p$G %*% pr$Q)
  edf <- sum(colSums(V * (Gq %*% V)) / ev$values[keep])
  list(
    edf = edf,
    converged = converged,
    grad_norm = sqrt(sum(pr$gq^2)),
    decrement = dec,
    n_newton = it,
    n_hessian = n_hess,
    n_negative = sum(ev$values[keep] < 0),
    n_flat = sum(!keep),
    npar = np,
    cores = cores,
    intercept = a0,
    eta = p$eta,
    mu = p$mu,
    time_s = proc.time()[["elapsed"]] - t0
  )
}

#' A ttps() fit moved to the polished cores of .tt_exact_edf_core() (internal).
#' @keywords internal
#' @noRd
.tt_fit_polished <- function(fit, e, y) {
  fit$cores <- e$cores
  fit$intercept <- e$intercept
  fit$fitted.values <- e$mu
  fit$linear.predictors <- e$eta
  fit$residuals <- as.numeric(y) - e$mu
  fit$deviance <- glm_deviance(fit$family, as.numeric(y), e$mu,
                               weights = fit$weights)
  fit$polished <- list(converged = e$converged, n_newton = e$n_newton,
                       grad_norm = e$grad_norm, n_negative = e$n_negative)
  fit
}

#' Exact effective degrees of freedom of a fitted TT P-spline
#'
#' Computes \eqn{\mathrm{tr}(d\hat\mu/dy)} of the fixed-\eqn{\lambda}
#' estimator by implicit differentiation with the full Hessian of the
#' penalized objective (see Details). Unlike the block-diagonal T-EDF stored
#' in `fit$edf`, it includes the curvature of the multilinear TT map and the
#' cross-core terms of the global penalty, and it reproduces the exact dense
#' EDF at saturating rank.
#'
#' @details The trace is the derivative of the converged estimator, so the
#'   cores are first polished to a stationary point of the penalized
#'   objective by a trust-region Newton iteration (`polish = TRUE`); ALS
#'   stopped at a tolerance can be far from stationary at low rank, or near
#'   a saddle. The Hessian is
#'   analytic (data term and multilinear curvature) plus exact central
#'   differences of the penalty gradient; TT gauge directions are projected
#'   out. Cost per Newton step: one stacked design (`n x npar`), one `J'WJ`
#'   and one dense eigendecomposition of size about `npar`, so it is meant for
#'   hundreds to a few thousand TT parameters. Gaussian (identity) and
#'   Poisson (log) families; no `linear=` / `smooth=` terms.
#' @param fit A [ttps()] fit (scattered or `array = TRUE`).
#' @param y,X Response and covariates of the fit (scattered fits; taken from
#'   the fit when stored). Ignored for array fits, which use `fit$Y` /
#'   `fit$axes` when available, or `Y` / `axes`.
#' @param Y,axes Array response and grid axes for array fits.
#' @param polish Polish the cores to a stationary point first (default
#'   `TRUE`). With `FALSE` the trace is taken at the fitted cores.
#' @param maxit,tol Trust-region Newton iterations and decrement tolerance
#'   (relative to the objective) of the polishing.
#' @return A list with `edf`, `converged`, `grad_norm` (projected gradient
#'   norm at the returned point), `decrement`, `n_newton`, `n_negative`
#'   (negative-curvature directions: the fit is a saddle), `n_flat`, `npar`,
#'   the polished `cores`, `intercept`, `eta`, `mu`, and `time_s`.
#' @seealso [tt_gdf()], [tt_gdf_array()] for Monte Carlo estimates.
#' @export
tt_edf_exact <- function(fit, y = NULL, X = NULL, Y = NULL, axes = NULL,
                         polish = TRUE, maxit = 12L, tol = 1e-12) {
  if (!is.null(fit$linear) || length(fit$smooth %||% list())) {
    stop("tt_edf_exact() does not support linear= / smooth= terms.", call. = FALSE)
  }
  if (isTRUE(fit$array %||% FALSE) || !is.null(Y)) {
    Y <- Y %||% fit$Y
    axes <- axes %||% fit$axes
    if (is.null(Y) || is.null(axes)) {
      stop("Array fit: supply Y and axes.", call. = FALSE)
    }
    X <- as.matrix(expand.grid(axes))
    y <- as.numeric(Y)
  } else {
    y <- y %||% fit$y
    X <- X %||% fit$X
    if (is.null(y) || is.null(X)) stop("Supply y and X.", call. = FALSE)
  }
  basis <- eval_marginal_bases(as.matrix(X), fit$knots, fit$degree,
                               cyclic = fit$cyclic)
  .tt_exact_edf_core(
    fit$cores, fit$intercept, basis, y, fit$family,
    offset = fit$offset, weights = fit$weights, lambda = fit$lambda,
    penalty_order = fit$penalty_order %||% 2L,
    polish = polish, maxit = maxit, tol = tol
  )
}
