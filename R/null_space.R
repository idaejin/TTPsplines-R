# Discrete P-spline null-space helpers (TODO-SC-NULL / SC §nullspace).
#
# Modes for ttps(null_space=):
#   joint    — single TT on full Theta (default)
#   profiled — Gaussian NSP: profile β0 exactly; ALS on Q0 y and Q0 Z_k

#' Basis for N(Δ^q) on coefficient vectors of length `k`.
#' @keywords internal
#' @noRd
tt_difference_null_basis <- function(k, order = 2L) {
  k <- as.integer(k)
  order <- as.integer(order)
  if (order < 1L) stop("penalty_order must be >= 1 for null-space split.", call. = FALSE)
  if (k <= order) {
    stop("Need k > penalty_order for a nontrivial penalized complement.",
         call. = FALSE)
  }
  D <- diff(diag(k), differences = order)
  sv <- svd(D, nu = 0L, nv = k)
  q <- order
  U <- sv$v[, (k - q + 1L):k, drop = FALSE]
  qrU <- qr(U)
  qr.Q(qrU)[, seq_len(q), drop = FALSE]
}

#' Row-wise tensor-product null design (n x q^d).
#' @keywords internal
#' @noRd
tt_null_space_design <- function(basis, U) {
  d <- length(basis)
  n <- nrow(basis[[1]])
  q <- ncol(U)
  BU <- lapply(seq_len(d), function(j) {
    as.matrix(basis[[j]] %*% U)
  })
  npar <- as.integer(q^d)
  Z <- matrix(1, n, npar)
  for (i in seq_len(n)) {
    v <- 1
    for (j in seq_len(d)) {
      v <- kronecker(BU[[j]][i, ], v)
    }
    Z[i, ] <- v
  }
  colnames(Z) <- paste0("null", seq_len(npar))
  Z
}

#' Build X0 = cbind(1, null tensor design) and check size / cyclic.
#' @keywords internal
#' @noRd
tt_null_space_X0 <- function(basis, penalty_order, max_npar = 4096L,
                             mode_label = "null_space") {
  cyclic <- attr(basis, "cyclic")
  if (!is.null(cyclic) && any(as.logical(cyclic))) {
    stop(
      sprintf("%s is not supported with cyclic margins yet.", mode_label),
      call. = FALSE
    )
  }
  k <- ncol(basis[[1]])
  d <- length(basis)
  U <- tt_difference_null_basis(k, order = penalty_order)
  q <- ncol(U)
  npar <- as.integer(q^d)
  if (npar > as.integer(max_npar)) {
    stop(
      sprintf(
        "%s needs q^d = %d parameters (q=%d, d=%d) > max_npar=%d.",
        mode_label, npar, q, d, as.integer(max_npar)
      ),
      call. = FALSE
    )
  }
  Z <- tt_null_space_design(basis, U)
  X0 <- cbind(`(Intercept)` = 1, Z)
  list(U = U, q = q, d = d, k = k, npar_tensor = npar, Z = Z, X0 = X0)
}

#' Empirical projector onto col(X0): Q0 = I - P0 (optional weights W).
#' @keywords internal
#' @noRd
tt_null_space_projector <- function(X0, weights = NULL) {
  X0 <- as.matrix(X0)
  n <- nrow(X0)
  w <- normalize_weights(weights, n)
  sw <- sqrt(w)
  Xw <- X0 * sw
  qrX <- qr(Xw, tol = 1e-10)
  rnk <- qrX$rank
  apply_vec <- function(v) {
    v <- as.numeric(v)
    beta <- qr.coef(qrX, sw * v)
    beta[!is.finite(beta)] <- 0
    as.numeric(v - as.numeric(X0 %*% beta))
  }
  apply_mat <- function(M) {
    M <- as.matrix(M)
    beta <- qr.coef(qrX, sw * M)
    beta[!is.finite(beta)] <- 0
    M - X0 %*% beta
  }
  list(
    X0 = X0,
    weights = w,
    qr = qrX,
    rank = as.integer(rnk),
    apply_vec = apply_vec,
    apply_mat = apply_mat
  )
}

#' Recover profiled null coefficients given TT surface.
#' @keywords internal
#' @noRd
tt_profile_null_coef <- function(y, offset, mu_tt, X0, weights = NULL) {
  n <- length(y)
  off <- normalize_offset(offset, n)
  w <- normalize_weights(weights, n)
  r <- as.numeric(y - off - mu_tt)
  fit <- stats::lm.wfit(x = as.matrix(X0), y = r, w = w)
  coef <- as.numeric(fit$coefficients)
  coef[!is.finite(coef)] <- 0
  coef
}

#' Predict null-space contribution at new bases.
#' @keywords internal
#' @noRd
tt_null_space_eta <- function(null_info, basis) {
  if (is.null(null_info)) return(rep(0, nrow(basis[[1]])))
  if (!identical(null_info$method %||% "", "profiled")) {
    return(rep(0, nrow(basis[[1]])))
  }
  Z <- tt_null_space_design(basis, null_info$U)
  Xns <- cbind(1, Z)
  coef <- null_info$coef
  if (length(coef) != ncol(Xns)) {
    stop("null_space coef length mismatch in prediction.", call. = FALSE)
  }
  as.numeric(Xns %*% coef)
}

#' Joint linearized EDF of the Q0-residualized TT map.
#' @keywords internal
#' @noRd
tt_joint_edf_profiled <- function(cores, basis, penalties, lambda, proj,
                                  weight = NULL, max_npar = 2500L) {
  tt_joint_edf_parts(
    cores, basis, penalties, lambda,
    weight = weight, max_npar = max_npar, null_proj = proj
  )$edf
}

#' Gaussian ALS with profiled (orthogonalized) discrete null space.
#'
#' Minimizes |Q0 (y - offset - μ_TT)|² + J_λ(μ_TT) over TT cores, then
#' recovers β0 = argmin |y - offset - X0 β0 - μ_TT|.
#' @keywords internal
#' @noRd
tt_als_fit_profiled_null <- function(y, basis, ranks, lambda_spec, control,
                                     penalty_order = 2L, init_cores = NULL,
                                     offset = NULL, weights = NULL,
                                     max_npar = 4096L) {
  method <- lambda_spec$method
  d <- length(basis)
  p <- ncol(basis[[1]])
  n <- length(y)
  offset <- normalize_offset(offset, n)
  w <- normalize_weights(weights, n)

  built <- tt_null_space_X0(
    basis, penalty_order, max_npar = max_npar,
    mode_label = "null_space = 'profiled'"
  )
  X0 <- built$X0
  proj <- tt_null_space_projector(X0, weights = w)
  # Working response in the orthogonal complement of col(X0)
  y_work <- proj$apply_vec(y - offset)

  lambda <- as.numeric(lambda_spec$values %||% lambda_spec$lambda0)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- tt_core_penalties_from_basis(ranks, basis, penalty_order)
  cyclic <- attr(basis, "cyclic")
  margin_order <- .cgcv_margin_order(control$cgcv_margin_order, d)
  rho <- control$cgcv_damping %||% 1
  delta <- control$cgcv_max_log10_step %||% Inf
  bounds <- control$lambda_bounds
  do_trace <- identical(method, "cGCV") && isTRUE(control$cgcv_trace %||% TRUE)

  # Force material design path so Q0 can residualize Z_k
  control_loc <- control
  control_loc$gram_method <- "legacy"

  n_eval <- 0L
  n_sweeps <- 0L
  history <- list()
  cgcv_trace <- list()
  prev_lam <- lambda
  prev_eta_tt <- NULL
  t0 <- proc.time()[["elapsed"]]
  intercept <- 0
  beta <- numeric(0)

  for (sw in seq_len(control$max_sweeps)) {
    use_cache <- isTRUE(control$design_interface_cache %||% TRUE)
    use_ltr <- use_cache && .tt_is_ltr_order(margin_order, d)
    use_rtl <- use_cache && .tt_is_rtl_order(margin_order, d)
    if (use_ltr) {
      R_all <- .tt_design_prepare_right(cores, basis)
      L_cur <- matrix(1, nrow(basis[[1]]), 1)
    } else if (use_rtl) {
      L_all <- .tt_design_prepare_left(cores, basis)
      R_cur <- matrix(1, nrow(basis[[1]]), 1)
    }
    for (k in margin_order) {
      if (use_ltr) {
        Left <- L_cur
        Right <- R_all[[k]]
      } else if (use_rtl) {
        Left <- L_all[[k]]
        Right <- R_cur
      } else {
        Left <- NULL
        Right <- NULL
      }
      built_ws <- .cgcv_core_workspace(
        cores, k, lambda, basis, y_work, ranks, control_loc,
        weight = w, penalty_order = penalty_order,
        use_spectral = identical(method, "cGCV") && isTRUE(control$use_spectral_gcv),
        compute_op_norms = do_trace,
        Left = Left,
        Right = Right,
        null_proj = proj
      )
      Pk <- built_ws$P_own
      penalties[[k]] <- Pk
      ws <- built_ws$workspace
      lambda_old_k <- lambda[k]
      upd <- update_lambda(method, ws)
      n_eval <- n_eval + upd$n_eval

      if (identical(method, "cGCV") && (rho < 1 - 1e-15 || is.finite(delta))) {
        step <- .cgcv_damped_trust_update(
          lambda_old = lambda_old_k,
          lambda_tilde = upd$lambda,
          rho = rho,
          max_log10_step = delta,
          bounds = bounds
        )
        lam_new <- step$lambda_new[[1L]]
        fit_clip <- .cgcv_eval_at(ws, lam_new)
        g_use <- fit_clip$g
        ed_use <- fit_clip$ed
        gcv_use <- fit_clip$value
        tilde <- upd$lambda
      } else {
        lam_new <- upd$lambda
        g_use <- upd$g
        ed_use <- upd$ed %||% NA_real_
        gcv_use <- upd$value %||% NA_real_
        tilde <- upd$lambda
      }
      cores[[k]] <- array(g_use, c(ranks[k], p, ranks[k + 1L]))
      lambda[k] <- lam_new

      if (use_ltr && k < d) {
        L_cur <- .tt_design_left_absorb(L_cur, cores[[k]], basis[[k]])
      } else if (use_rtl && k > 1L) {
        R_cur <- .tt_design_right_absorb(R_cur, cores[[k]], basis[[k]])
      }

      if (do_trace) {
        cgcv_trace[[length(cgcv_trace) + 1L]] <- data.frame(
          sweep = sw, margin = k,
          lambda_old = lambda_old_k, lambda_tilde = tilde, lambda_new = lam_new,
          ed = ed_use, gcv = gcv_use, mode = "profiled_null",
          stringsAsFactors = FALSE
        )
      }
    }

    mu_tt <- as.numeric(tt_contraction(cores, basis))
    # Profiled RSS on the original scale
    coef0 <- tt_profile_null_coef(y, offset, mu_tt, X0, weights = w)
    eta <- as.numeric(offset + as.numeric(X0 %*% coef0) + mu_tt)
    rss <- sum(w * (y - eta)^2)
    pen_val <- tt_global_penalty_value(
      cores, lambda, penalty_order = penalty_order, cyclic = cyclic
    )
    obj <- 0.5 * rss + pen_val
    d_eta <- if (sw == 1L) NA_real_ else sqrt(mean((mu_tt - prev_eta_tt)^2))
    history[[sw]] <- list(
      sweep = sw, rss = rss, objective = obj, penalty = pen_val,
      lambda = lambda, d_eta = d_eta
    )
    prev_eta_tt <- mu_tt
    n_sweeps <- sw
    if (control$trace) {
      cat(sprintf(
        "  ALS-profiled sweep %2d | obj=%.6g | RSS=%.6g | lambda=%s\n",
        sw, obj, rss, paste(sprintf("%.3g", lambda), collapse = ",")
      ))
    }
    if (identical(method, "cGCV")) {
      dlog <- max(abs(log(lambda) - log(pmax(prev_lam, 1e-12))))
      if (dlog < control$tol_lambda && sw > 1L) break
      prev_lam <- lambda
    } else if (sw > 2L) {
      prev_rss <- history[[sw - 1L]]$rss
      if (abs(prev_rss - rss) / max(1, abs(prev_rss)) < control$tol) break
    }
  }

  mu_tt <- as.numeric(tt_contraction(cores, basis))
  coef0 <- tt_profile_null_coef(y, offset, mu_tt, X0, weights = w)
  eta_null <- as.numeric(X0 %*% coef0)
  eta <- as.numeric(offset + eta_null + mu_tt)
  mu_tt_perp <- proj$apply_vec(mu_tt)
  orth <- as.numeric(crossprod(X0, w * mu_tt_perp))
  # In-sample NSP identity: yhat = offset + P0(y-offset) + Q0 mu_tt
  y_c <- y - offset
  eta_id <- as.numeric(offset + y_c - proj$apply_vec(y_c) + mu_tt_perp)
  cgcv_df <- if (length(cgcv_trace)) do.call(rbind, cgcv_trace) else NULL

  null_info <- list(
    U = built$U,
    q = built$q,
    d = built$d,
    k = built$k,
    coef = coef0,
    npar = length(coef0),
    npar_tensor = built$npar_tensor,
    eta = eta_null,
    design_rank = proj$rank,
    method = "profiled",
    X0 = X0,
    projector = proj,
    mu_tt = mu_tt,
    mu_tt_perp = mu_tt_perp,
    ortho_X0_mu_tt_perp = orth,
    # legacy alias used in early drafts
    ortho_X0_mu_tt = orth,
    identity_max_abs = max(abs(eta - eta_id))
  )

  list(
    cores = cores,
    intercept = 0,
    beta = numeric(0),
    linear = NULL,
    smooth = NULL,
    lambda = lambda,
    ranks = ranks,
    eta = eta,
    mu = eta,
    deviance = sum(w * (y - eta)^2),
    n_sweeps = n_sweeps,
    n_pirls = NA_integer_,
    n_criterion_evals = n_eval,
    history = history,
    penalties = penalties,
    penalty_mode = "global",
    cgcv = list(
      update = if (identical(method, "cGCV")) "sequential_profiled" else NA_character_,
      trace = cgcv_df
    ),
    q_descent = list(checked = FALSE),
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = TRUE,
    method_lambda = method,
    optimizer = "ALS",
    backend = "R",
    null_space_info = null_info
  )
}
