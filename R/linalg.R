# Internal linear algebra helpers (no explicit inverses).

solve_spd <- function(a, b) {
  chol_a <- tryCatch(chol(a), error = function(e) NULL)
  if (is.null(chol_a)) {
    return(solve(a, b))
  }
  backsolve(chol_a, forwardsolve(t(chol_a), b))
}

# Minimum-norm solution of M g = b for a symmetric positive semi-definite M,
# by one symmetric eigendecomposition: g = sum of v_i (v_i' b) / l_i over the
# eigenvalues l_i above rel_tol * max eigenvalue. Used for every global-mode
# TT core system (update_lambda_fixed()). At saturating rank under heavy
# smoothing the other cores' interfaces are rank deficient and M is singular
# along gauge directions (core changes that leave Theta unchanged; rcond
# down to 1e-21). There a Cholesky -> LU -> ridge cascade succeeds or fails
# by roundoff and leaves a gauge component of roundoff over a near-zero
# pivot, so the fixed-lambda ALS map jumps under tiny changes of y and
# finite-difference GDFs blow up. One formula for every system instead: on
# a well-conditioned M (every eigenvalue kept) it is the exact solve, and
# gauge directions get no component. A ridge (M + delta I) also removes the
# jumps but biases every conditional solve and slows ALS at saturating rank;
# a smooth cutoff of the eigenvalues makes the fixed-budget map steep
# wherever an eigenvalue sits inside the cutoff band.
solve_psd_pinv <- function(M, b, rel_tol = 1e-12) {
  if (!all(is.finite(M)) || !all(is.finite(b))) {
    stop("Core system has non-finite entries.", call. = FALSE)
  }
  e <- eigen(0.5 * (M + t(M)), symmetric = TRUE)
  lmax <- e$values[1L]
  if (!(lmax > 0)) return(matrix(0, nrow(M), NCOL(b)))
  keep <- e$values > rel_tol * lmax
  V <- e$vectors[, keep, drop = FALSE]
  V %*% (crossprod(V, b) / e$values[keep])
}

solve_spd_ridge <- function(A, b, base_ridge = NULL) {
  m <- nrow(A)
  if (is.null(base_ridge)) base_ridge <- ridge_scale(A, multiplier = 1e-6)
  for (fac in c(1, 10, 1e2, 1e3, 1e4, 1e5)) {
    out <- tryCatch(
      solve_spd(A + (fac * base_ridge) * diag(m), b),
      error = function(e) NULL
    )
    if (!is.null(out) && all(is.finite(out))) return(out)
  }
  qr.solve(A + 1e-3 * mean(diag(A)) * diag(m), b)
}

ridge_scale <- function(xtx, multiplier = 1e-7) {
  scale <- mean(diag(xtx))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  multiplier * scale
}
