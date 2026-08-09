#' Difference penalty D'D for P-splines.
#' @keywords internal
difference_penalty <- function(p, order = 2) {
  crossprod(diff(diag(p), differences = order))
}

#' Core-wise Kronecker penalty for vec(G_k) with order (a, j, b).
#' @keywords internal
core_penalty <- function(rl, p, rr, penalty_order = 2) {
  DtD <- difference_penalty(p, penalty_order)
  kronecker(diag(rr), kronecker(DtD, diag(rl)))
}

#' Block-diagonal joint penalty from per-core penalties.
#' @keywords internal
block_penalty <- function(penalties, lambda) {
  sizes <- vapply(penalties, nrow, integer(1))
  m <- sum(sizes)
  out <- matrix(0, m, m)
  offset <- 0L
  for (k in seq_along(penalties)) {
    idx <- offset + seq_len(sizes[k])
    out[idx, idx] <- lambda[k] * penalties[[k]]
    offset <- offset + sizes[k]
  }
  out
}
