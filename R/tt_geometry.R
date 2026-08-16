# TT geometry: contractions, interfaces, conditional designs.

initialize_tt_cores <- function(p, ranks, seed = 1, sd = 0.15) {
  set.seed(seed)
  d <- length(ranks) - 1L
  p <- rep(as.integer(p), length.out = d)
  cores <- vector("list", d)
  for (k in seq_len(d)) {
    cores[[k]] <- array(
      stats::rnorm(ranks[k] * p[k] * ranks[k + 1L], 0, sd),
      c(ranks[k], p[k], ranks[k + 1L])
    )
  }
  cores
}

contract_left_step <- function(left, core, Bk) {
  if (exists("contract_left_step_cpp", mode = "function")) {
    return(contract_left_step_cpp(left, core, Bk))
  }
  n <- nrow(Bk)
  rl <- dim(core)[1]
  p <- dim(core)[2]
  rr <- dim(core)[3]
  out <- matrix(0, n, rr)
  for (j in seq_len(p)) {
    Cj <- matrix(core[, j, ], rl, rr)
    out <- out + (left %*% Cj) * Bk[, j]
  }
  out
}

contract_right_step <- function(right, core, Bk) {
  if (exists("contract_right_step_cpp", mode = "function")) {
    return(contract_right_step_cpp(right, core, Bk))
  }
  n <- nrow(Bk)
  rl <- dim(core)[1]
  p <- dim(core)[2]
  rr <- dim(core)[3]
  out <- matrix(0, n, rl)
  for (j in seq_len(p)) {
    Cj <- matrix(core[, j, ], rl, rr)
    out <- out + (right %*% t(Cj)) * Bk[, j]
  }
  out
}

tt_contraction <- function(cores, basis) {
  n <- nrow(basis[[1]])
  cur <- matrix(1, n, 1)
  for (k in seq_along(cores)) {
    cur <- contract_left_step(cur, cores[[k]], basis[[k]])
  }
  drop(cur)
}

left_interfaces <- function(cores, basis) {
  d <- length(cores)
  n <- nrow(basis[[1]])
  L <- vector("list", d)
  L[[1]] <- matrix(1, n, 1)
  cur <- L[[1]]
  if (d >= 2L) {
    for (k in seq_len(d - 1L)) {
      cur <- contract_left_step(cur, cores[[k]], basis[[k]])
      L[[k + 1L]] <- cur
    }
  }
  L
}

right_interfaces <- function(cores, basis) {
  d <- length(cores)
  n <- nrow(basis[[1]])
  R <- vector("list", d)
  R[[d]] <- matrix(1, n, 1)
  cur <- R[[d]]
  if (d >= 2L) {
    for (k in d:2) {
      cur <- contract_right_step(cur, cores[[k]], basis[[k]])
      R[[k - 1L]] <- cur
    }
  }
  R
}

#' Conditional design for vec(core_k); order a (fast), j, b (slow).
#' @keywords internal
tt_design_core <- function(Left, Right, Bk) {
  if (exists("tt_design_core_d_cpp", mode = "function")) {
    return(tt_design_core_d_cpp(Left, Right, Bk))
  }
  n <- nrow(Bk)
  rl <- ncol(Left)
  rr <- ncol(Right)
  p <- ncol(Bk)
  X <- matrix(0, n, rl * p * rr)
  col <- 1L
  for (b in seq_len(rr)) {
    for (j in seq_len(p)) {
      for (a in seq_len(rl)) {
        X[, col] <- Left[, a] * Bk[, j] * Right[, b]
        col <- col + 1L
      }
    }
  }
  X
}

#' TRUE if `order` is left-to-right 1:d (enables design L/R sweep cache).
#' @keywords internal
#' @noRd
.tt_is_ltr_order <- function(order, d) {
  identical(as.integer(order), seq_len(as.integer(d)))
}

#' TRUE if `order` is right-to-left d:1.
#' @keywords internal
#' @noRd
.tt_is_rtl_order <- function(order, d) {
  identical(as.integer(order), rev(seq_len(as.integer(d))))
}

#' Prepare right interfaces once for a left-to-right ALS sweep.
#'
#' During LTR, cores \(G_{k+1},\ldots,G_d\) are not yet updated when core \(k\)
#' is fitted, so the full right chain remains valid for the entire sweep.
#' @keywords internal
#' @noRd
.tt_design_prepare_right <- function(cores, basis) {
  right_interfaces(cores, basis)
}

#' Prepare left interfaces once for a right-to-left ALS sweep.
#' @keywords internal
#' @noRd
.tt_design_prepare_left <- function(cores, basis) {
  left_interfaces(cores, basis)
}

#' Absorb updated core into the running left interface (LTR).
#' @keywords internal
#' @noRd
.tt_design_left_absorb <- function(Left, core, Bk) {
  contract_left_step(Left, core, Bk)
}

#' Absorb updated core into the running right interface (RTL).
#' @keywords internal
#' @noRd
.tt_design_right_absorb <- function(Right, core, Bk) {
  contract_right_step(Right, core, Bk)
}

#' Reconstruct full coefficient tensor (dangerous for large p^d).
#' @keywords internal
tt_full_theta <- function(cores) {
  # Sequential contraction over TT cores into a dense array — only small d,p
  d <- length(cores)
  p <- vapply(cores, function(a) dim(a)[2], integer(1))
  if (prod(p) > 5e5) {
    stop("Full Theta reconstruction refused: prod(p)=", prod(p), call. = FALSE)
  }
  # Convert TT to full tensor via successive multiplications
  # G1: 1 x p1 x r1 -> matrix p1 x r1
  g <- matrix(cores[[1]][1, , ], p[1], dim(cores[[1]])[3])
  if (d == 1L) return(array(g, p))
  for (k in 2:d) {
    rl <- dim(cores[[k]])[1]
    pk <- p[k]
    rr <- dim(cores[[k]])[3]
    # g is (prod p_1..p_{k-1}) x rl
    next_mat <- matrix(0, nrow(g) * pk, rr)
    for (j in seq_len(pk)) {
      Cj <- matrix(cores[[k]][, j, ], rl, rr)
      block <- g %*% Cj
      rows <- ((j - 1L) * nrow(g) + 1L):(j * nrow(g))
      next_mat[rows, ] <- block
    }
    g <- next_mat
  }
  array(g[, 1], p)
}
