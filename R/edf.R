# Joint linearized EDF for TT P-splines at convergence.
# edf = tr[(J'J + P_λ)^{-1} J'J], J = [X_1 | … | X_d] stacked conditional designs.
# See lab docs/COORDINATE_AIC_BIC_LAMBDA.md. Gauge directions contribute little;
# do not confuse with npar_TT.

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
  # base::bdiag needs Matrix; build densely for modest npar
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

#' Joint effective degrees of freedom (linearized TT map).
#'
#' @param cores TT cores at convergence.
#' @param basis Marginal B-spline bases.
#' @param penalties List of core penalty matrices.
#' @param lambda Length-`d` smoothing vector.
#' @param weight Optional PIRLS weights (sqrt-weights applied to rows of J).
#' @param max_npar Skip if packed TT parameter count exceeds this (memory).
#' @return Numeric EDF, or `NA_real_` if skipped / failed.
#' @keywords internal
tt_joint_edf <- function(cores, basis, penalties, lambda,
                         weight = NULL, max_npar = 2500L) {
  npar <- sum(vapply(cores, length, integer(1)))
  if (npar > as.integer(max_npar)) return(NA_real_)
  n <- nrow(basis[[1]])
  # Rough memory guard for dense J (n × npar doubles)
  if (as.numeric(n) * as.numeric(npar) > 5e7) return(NA_real_)

  J <- tryCatch(
    tt_stacked_jacobian(cores, basis, weight = weight),
    error = function(e) NULL
  )
  if (is.null(J)) return(NA_real_)
  P <- tryCatch(
    tt_block_penalty(penalties, lambda),
    error = function(e) NULL
  )
  if (is.null(P) || nrow(P) != ncol(J)) return(NA_real_)

  ed <- tryCatch({
    if (exists("effective_df_cpp", mode = "function")) {
      as.numeric(effective_df_cpp(J, P))
    } else {
      .effective_df_r(J, P)
    }
  }, error = function(e) NA_real_)

  if (!is.finite(ed) || ed < 0) return(NA_real_)
  ed
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
