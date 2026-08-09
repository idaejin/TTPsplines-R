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

#' Penalty value and per-core gradients: 0.5 sum_k λ_k g'P g.
#' @keywords internal
#' @noRd
.tt_penalty_value_grad <- function(cores, penalties, lambda) {
  val <- 0
  grads <- vector("list", length(cores))
  for (k in seq_along(cores)) {
    g <- as.numeric(cores[[k]])
    Pk <- penalties[[k]]
    Pg <- as.numeric(Pk %*% g)
    val <- val + 0.5 * lambda[k] * sum(g * Pg)
    grads[[k]] <- array(lambda[k] * Pg, dim(cores[[k]]))
  }
  list(value = val, grads = grads)
}
