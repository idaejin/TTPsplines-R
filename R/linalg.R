# Internal linear algebra helpers (no explicit inverses).

solve_spd <- function(a, b) {
  chol_a <- tryCatch(chol(a), error = function(e) NULL)
  if (is.null(chol_a)) {
    return(solve(a, b))
  }
  backsolve(chol_a, forwardsolve(t(chol_a), b))
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
