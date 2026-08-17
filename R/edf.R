# Joint linearized EDF for TT P-splines at convergence.
# edf = tr(H), H = (J'J + P_λ)^{-1} J'J, J = [X_1 | … | X_d] stacked designs.
# See lab docs/COORDINATE_AIC_BIC_LAMBDA.md.
#
# fit$edf_margin[k] = tr(H_kk): diagonal block of H for core-k parameters.
#   These SUM to the joint EDF (parameter-block partition).
# fit$edf_margin_cond[k]: conditional ALS EDF for core k (cGCV diagnostic).
#   These do NOT sum to the joint EDF.

#' Build stacked Jacobian of the TT map w.r.t. packed cores.
#' @keywords internal
tt_stacked_jacobian <- function(cores, basis, weight = NULL) {
  d <- length(cores)
  L <- left_interfaces(cores, basis)
  R <- right_interfaces(cores, basis)
  blocks <- vector("list", d)
  for (k in seq_len(d)) {
    Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
    if (!is.null(weight)) {
      sw <- sqrt(pmax(as.numeric(weight), 0))
      Xk <- Xk * sw
    }
    blocks[[k]] <- Xk
  }
  do.call(cbind, blocks)
}

#' Block-diagonal TT penalty P_λ = blkdiag(λ_1 P_1, …, λ_d P_d).
#' @keywords internal
tt_block_penalty <- function(penalties, lambda) {
  d <- length(penalties)
  stopifnot(length(lambda) == d)
  blocks <- lapply(seq_len(d), function(k) as.numeric(lambda[k]) * penalties[[k]])
  m <- sum(vapply(blocks, nrow, integer(1)))
  P <- matrix(0, m, m)
  off <- 0L
  for (k in seq_len(d)) {
    mk <- nrow(blocks[[k]])
    idx <- (off + 1L):(off + mk)
    P[idx, idx] <- blocks[[k]]
    off <- off + mk
  }
  P
}

#' Packed core parameter counts (column blocks of J / H).
#' @keywords internal
#' @noRd
tt_core_block_sizes <- function(cores) {
  vapply(cores, length, integer(1))
}

#' Influence matrix H = (J'J + P)^{-1} J'J in packed TT coordinates.
#' @keywords internal
#' @noRd
tt_influence_H <- function(J, P) {
  xtx <- crossprod(J)
  ridge <- ridge_scale(xtx, multiplier = 1e-9)
  m <- nrow(xtx)
  system <- xtx + P + ridge * diag(m)
  solve_spd(system, xtx)
}

#' Diagonal-block traces of a square matrix.
#' @keywords internal
#' @noRd
tt_block_diag_traces <- function(H, sizes, names = NULL) {
  sizes <- as.integer(sizes)
  d <- length(sizes)
  out <- rep(NA_real_, d)
  if (!is.null(names) && length(names) == d) {
    names(out) <- names
  } else {
    names(out) <- paste0("margin", seq_len(d))
  }
  if (is.null(H) || !is.matrix(H) || nrow(H) != sum(sizes)) return(out)
  off <- 0L
  for (k in seq_len(d)) {
    mk <- sizes[[k]]
    if (mk <= 0L) {
      out[[k]] <- 0
      next
    }
    idx <- (off + 1L):(off + mk)
    out[[k]] <- sum(diag(H[idx, idx, drop = FALSE]))
    off <- off + mk
  }
  out
}

#' Joint EDF and additive block-margin EDFs from one influence matrix.
#'
#' @param null_proj Optional projector (e.g. profiled \(Q_0\)) applied to rows
#'   of \(J\) before forming \(H\).
#' @return List with `edf`, `edf_margin` (block traces), `ok`.
#' @keywords internal
#' @noRd
tt_joint_edf_parts <- function(cores, basis, penalties, lambda,
                               weight = NULL, max_npar = 2500L,
                               null_proj = NULL, names = NULL) {
  npar <- sum(vapply(cores, length, integer(1)))
  empty <- list(
    edf = NA_real_,
    edf_margin = {
      d <- length(cores)
      v <- rep(NA_real_, d)
      if (!is.null(names) && length(names) == d) names(v) <- names
      else names(v) <- paste0("margin", seq_len(d))
      v
    },
    ok = FALSE
  )
  if (npar > as.integer(max_npar)) return(empty)
  n <- nrow(basis[[1]])
  if (as.numeric(n) * as.numeric(npar) > 5e7) return(empty)

  J <- tryCatch(
    tt_stacked_jacobian(cores, basis, weight = weight),
    error = function(e) NULL
  )
  if (is.null(J)) return(empty)
  if (!is.null(null_proj)) {
    J <- null_proj$apply_mat(J)
  }
  P <- tryCatch(
    tt_block_penalty(penalties, lambda),
    error = function(e) NULL
  )
  if (is.null(P) || nrow(P) != ncol(J)) return(empty)

  parts <- tryCatch({
    H <- tt_influence_H(J, P)
    sizes <- tt_core_block_sizes(cores)
    list(
      edf = sum(diag(H)),
      edf_margin = tt_block_diag_traces(H, sizes, names = names),
      ok = TRUE
    )
  }, error = function(e) empty)

  if (!isTRUE(parts$ok)) return(empty)
  if (!is.finite(parts$edf) || parts$edf < 0) {
    parts$edf <- NA_real_
    parts$ok <- FALSE
  }
  parts
}

#' Joint effective degrees of freedom (linearized TT map).
#'
#' @keywords internal
tt_joint_edf <- function(cores, basis, penalties, lambda,
                         weight = NULL, max_npar = 2500L) {
  tt_joint_edf_parts(
    cores, basis, penalties, lambda,
    weight = weight, max_npar = max_npar
  )$edf
}

#' Additive per-margin EDF = diagonal blocks of the joint influence \(H\).
#'
#' \(\mathrm{edf}_k=\operatorname{tr}(H_{kk})\) with
#' \(H=(J^\top J+P_\lambda)^{-1}J^\top J\). These sum to the joint EDF.
#'
#' @keywords internal
#' @noRd
tt_margin_edf <- function(cores, basis, penalties, lambda,
                          weight = NULL, names = NULL, max_npar = 2500L,
                          null_proj = NULL) {
  tt_joint_edf_parts(
    cores, basis, penalties, lambda,
    weight = weight, max_npar = max_npar,
    null_proj = null_proj, names = names
  )$edf_margin
}

#' Conditional EDF for each TT margin (ALS core design; cGCV diagnostic).
#'
#' For each margin, the conditional EDF is the trace of the core-level
#' smoother at the fitted interfaces (same quantity used inside conditional
#' cGCV updates). These do **not** sum to the joint EDF.
#'
#' @keywords internal
#' @noRd
tt_margin_edf_cond <- function(cores, basis, penalties, lambda,
                               weight = NULL, names = NULL) {
  d <- length(cores)
  stopifnot(length(penalties) == d, length(lambda) == d)
  out <- rep(NA_real_, d)
  if (!is.null(names) && length(names) == d) {
    names(out) <- names
  } else {
    names(out) <- paste0("margin", seq_len(d))
  }

  L <- tryCatch(left_interfaces(cores, basis), error = function(e) NULL)
  R <- tryCatch(right_interfaces(cores, basis), error = function(e) NULL)
  if (is.null(L) || is.null(R)) return(out)

  sw <- NULL
  if (!is.null(weight)) {
    sw <- sqrt(pmax(as.numeric(weight), 0))
  }

  for (k in seq_len(d)) {
    out[[k]] <- tryCatch({
      Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
      if (!is.null(sw)) Xk <- Xk * sw
      Pk <- as.numeric(lambda[[k]]) * penalties[[k]]
      if (exists("effective_df_cpp", mode = "function")) {
        as.numeric(effective_df_cpp(Xk, Pk))
      } else {
        .effective_df_r(Xk, Pk)
      }
    }, error = function(e) NA_real_)
  }
  out
}

.effective_df_r <- function(jacobian, penalty) {
  xtx <- crossprod(jacobian)
  ridge <- ridge_scale(xtx, multiplier = 1e-9)
  m <- nrow(xtx)
  system <- xtx + penalty + ridge * diag(m)
  infl <- solve_spd(system, xtx)
  sum(diag(infl))
}

#' Working weights at a fitted object for GLM joint EDF.
#' @keywords internal
tt_edf_weights <- function(fit_raw, family_key, y, weights = NULL) {
  w_obs <- normalize_weights(weights, length(y))
  if (identical(family_key, "gaussian")) {
    if (all(w_obs == 1)) return(NULL)
    return(w_obs)
  }
  fam <- normalize_family(if (identical(family_key, "bernoulli")) "binomial" else family_key)
  work <- glm_working(fam, y, fit_raw$eta)
  work$weight * w_obs
}

#' Extract joint and per-margin EDF from a TT fit
#'
#' Returns the joint linearized EDF (`fit$edf`) and two per-margin summaries:
#'
#' \describe{
#'   \item{`margin` / `fit$edf_margin`}{Diagonal blocks
#'     \(\operatorname{tr}(H_{kk})\) of the joint influence
#'     \(H=(J^\top J+P_\lambda)^{-1}J^\top J\). These **sum** to the joint EDF
#'     (parameter-block partition of TT core coordinates — not univariate
#'     GAM-style margin EDFs).}
#'   \item{`margin_cond` / `fit$edf_margin_cond`}{Conditional ALS core EDFs
#'     used inside cGCV. Diagnostic only; \(\sum_k\mathrm{ed}_k\) need not
#'     equal `edf`.}
#' }
#'
#' If `compute_edf` was `FALSE` at fit time, recomputes when `cores`,
#' `penalties`, and bases are available (rebuilds bases from `fit$X`).
#'
#' @param object A `"ttpspline"` / `"ttps"` fit from [ttps()].
#' @param recompute If `TRUE`, always recompute from cores (default `FALSE`
#'   uses stored slots when present).
#' @param max_npar Memory guard for joint EDF (default from
#'   `object$control$edf_max_npar` or `2500`).
#'
#' @return A list of class `"tt_edf"` with `joint`, `margin`, `margin_cond`,
#'   `sum_margin`, `note`, and `additive` (`TRUE` for block margins).
#'
#' @examples
#' data(ishigami)
#' X <- as.matrix(ishigami[1:150, c("x1", "x2", "x3")])
#' fit <- ttps(
#'   ishigami$y[1:150], X, rank = 2, k = 6, lambda = 1,
#'   control = tt_control(max_sweeps = 5, compute_edf = TRUE)
#' )
#' tt_edf(fit)
#' # sum(fit$edf_margin) ~= fit$edf
#' # fit$edf_margin_cond  # cGCV-style conditional traces
#'
#' @export
tt_edf <- function(object, recompute = FALSE, max_npar = NULL) {
  if (!inherits(object, "ttpspline") && !inherits(object, "ttps")) {
    stop("tt_edf() expects a ttps() / ttpspline fit.", call. = FALSE)
  }
  max_npar <- max_npar %||% object$control$edf_max_npar %||% 2500L
  max_npar <- as.integer(max_npar)

  joint <- object$edf
  margin <- object$edf_margin
  margin_cond <- object$edf_margin_cond
  note <- object$edf_note %||% "stored"

  need <- isTRUE(recompute) ||
    !is.finite(joint) ||
    is.null(margin) ||
    !length(margin) ||
    !all(is.finite(margin))

  if (need) {
    if (is.null(object$cores) || is.null(object$penalties) ||
        is.null(object$X) || is.null(object$knots)) {
      stop("Cannot (re)compute EDF: fit lacks cores/penalties/X/knots.",
           call. = FALSE)
    }
    bs <- build_marginal_bases(
      object$X, k = object$k, degree = object$degree, knots = object$knots,
      cyclic = object$cyclic %||% FALSE
    )
    basis <- bs$basis
    w <- NULL
    if (!identical(object$family_key %||% "gaussian", "gaussian") ||
        !is.null(object$weights)) {
      w <- tt_edf_weights(
        list(eta = object$linear.predictors, mu = object$fitted.values),
        object$family_key %||% "gaussian",
        object$y,
        weights = object$weights
      )
    }
    mnames <- colnames(object$X)
    null_proj <- NULL
    if (identical(object$null_space, "profiled") &&
        !is.null(object$null_space_info$projector)) {
      null_proj <- object$null_space_info$projector
    }
    parts <- tt_joint_edf_parts(
      object$cores, basis, object$penalties, as.numeric(object$lambda),
      weight = w, max_npar = max_npar, null_proj = null_proj, names = mnames
    )
    joint <- parts$edf
    if (identical(object$null_space, "profiled") &&
        !is.null(object$null_space_info$design_rank) &&
        is.finite(joint)) {
      joint <- as.integer(object$null_space_info$design_rank) + joint
    }
    margin <- parts$edf_margin
    margin_cond <- tt_margin_edf_cond(
      object$cores, basis, object$penalties, as.numeric(object$lambda),
      weight = w, names = mnames
    )
    note <- if (isTRUE(recompute)) {
      "recomputed joint + block margin EDF (+ conditional diagnostic)"
    } else {
      "computed on demand (was missing / incomplete on fit)"
    }
  }

  sum_m <- if (length(margin) && all(is.finite(margin))) sum(margin) else NA_real_
  # For profiled fits, block margins sum to TT-perp EDF, not GDF_total
  additive_ok <- is.finite(joint) && is.finite(sum_m)
  if (identical(object$null_space, "profiled") &&
      !is.null(object$edf_tt) && is.finite(object$edf_tt)) {
    additive_ok <- is.finite(sum_m) &&
      abs(sum_m - object$edf_tt) <= 1e-4 * max(1, abs(object$edf_tt)) + 1e-6
  } else if (additive_ok) {
    additive_ok <- abs(sum_m - joint) <= 1e-4 * max(1, abs(joint)) + 1e-6
  }

  structure(
    list(
      joint = as.numeric(joint),
      margin = margin,
      margin_cond = margin_cond,
      sum_margin = sum_m,
      note = note,
      additive = isTRUE(additive_ok),
      warning = if (identical(object$null_space, "profiled")) {
        "edf_margin blocks sum to GDF_TT_perp (edf_tt), not GDF_total = rank(X0)+edf_tt."
      } else {
        NULL
      }
    ),
    class = "tt_edf"
  )
}

#' @export
print.tt_edf <- function(x, digits = 3, ...) {
  cat("TT effective degrees of freedom\n")
  cat(sprintf("  Joint (linearized): %s\n",
              if (is.finite(x$joint)) format(round(x$joint, digits), nsmall = digits)
              else "NA"))
  cat("  Per-margin (block tr(H_kk); additive partition of TT map):\n")
  print(round(x$margin, digits))
  if (is.finite(x$sum_margin)) {
    cat(sprintf(
      "  sum(margin) = %s%s\n",
      format(round(x$sum_margin, digits), nsmall = digits),
      if (isTRUE(x$additive)) "  [= joint TT EDF]" else ""
    ))
  }
  if (!is.null(x$margin_cond) && length(x$margin_cond) &&
      any(is.finite(x$margin_cond))) {
    cat("  Per-margin conditional (ALS / cGCV diagnostic; not additive):\n")
    print(round(x$margin_cond, digits))
  }
  if (!is.null(x$note) && nzchar(x$note)) {
    cat("  Note:", x$note, "\n")
  }
  if (!is.null(x$warning) && nzchar(x$warning)) {
    cat("  Warning:", x$warning, "\n")
  }
  invisible(x)
}
