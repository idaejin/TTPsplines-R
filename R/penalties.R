#' Difference penalty D'D for P-splines.
#' @keywords internal
difference_penalty <- function(p, order = 2) {
  crossprod(diff(diag(p), differences = order))
}

#' Circular (periodic) difference penalty for cyclic P-splines.
#' @keywords internal
#' @noRd
circular_difference_penalty <- function(p, order = 2) {
  p <- as.integer(p)
  order <- as.integer(order)
  if (p < 2L) stop("circular penalty needs p >= 2.", call. = FALSE)
  D1 <- diag(-1, p)
  D1[cbind(seq_len(p), c(seq_len(p)[-1L], 1L))] <- 1
  D <- D1
  if (order > 1L) {
    for (i in seq_len(order - 1L)) D <- D1 %*% D
  }
  crossprod(D)
}

#' Core-wise Kronecker penalty for vec(G_k) with order (a, j, b).
#' @keywords internal
core_penalty <- function(rl, p, rr, penalty_order = 2, cyclic = FALSE) {
  DtD <- if (isTRUE(cyclic)) {
    circular_difference_penalty(p, penalty_order)
  } else {
    difference_penalty(p, penalty_order)
  }
  kronecker(diag(rr), kronecker(DtD, diag(rl)))
}

#' Per-margin TT core penalties (optional cyclic flags).
#' @keywords internal
#' @noRd
tt_core_penalties <- function(ranks, p, penalty_order = 2, cyclic = NULL) {
  d <- length(ranks) - 1L
  cyclic <- normalize_cyclic(cyclic, d)
  lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p, ranks[k + 1L], penalty_order, cyclic = cyclic[k])
  })
}

#' Build core penalties using `attr(basis, "cyclic")` when present.
#' @keywords internal
#' @noRd
tt_core_penalties_from_basis <- function(ranks, basis, penalty_order = 2) {
  p <- ncol(basis[[1]])
  cyclic <- attr(basis, "cyclic")
  tt_core_penalties(ranks, p, penalty_order, cyclic = cyclic)
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
