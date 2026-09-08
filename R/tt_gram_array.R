# Array-mode Gram and RHS for TT-ALS (Gaussian, no X materialisation).
#
# When data Y live on a complete d-way grid (n_1 x ... x n_d) and the basis
# for each margin is B[[k]] (n_k x p_k), the conditional Gram and RHS for
# core k can be computed WITHOUT forming the n x (r_l * p * r_r) design matrix.
#
# Key identities (unweighted Gaussian):
#
#   S_k = kron(R_uniq' R_uniq, kron(Bk' Bk, L_uniq' L_uniq))
#
#   b_k = triple-mode contraction of Y_centered over (L_uniq, Bk, R_uniq)
#
# where L_uniq (n_left x r_l) and R_uniq (n_right x r_r) are the unique rows
# of the TT left/right interfaces, extracted from the scattered interfaces
# (which contain n_left and n_right repeated blocks respectively, when data
# come from a complete grid with dim1 varying fastest -- as in expand.grid).
#
# The saving is O(n_total) -> O(n_left + p + n_right) for the Gram, avoiding
# the n_total x q_k matrix entirely.

# --------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------

#' Extract unique rows of Left/Right interfaces for array-mode gram.
#'
#' For a complete d-way grid stored with dim1 fastest (expand.grid convention),
#' the left interface L[[k]] has n_left unique rows (repeated n_k * n_right times)
#' and the right interface R[[k]] has n_right unique rows.
#'
#' @param k Margin index (1-based).
#' @param L_all List of left interfaces from [left_interfaces()].
#' @param R_all List of right interfaces from [right_interfaces()].
#' @param n_grid Integer vector of grid sizes (n_1, ..., n_d).
#' @return List with `L` (n_left x r_l), `R` (n_right x r_r).
#' @keywords internal
#' @noRd
.tt_array_extract_interfaces <- function(k, L_all, R_all, n_grid) {
  d       <- length(n_grid)
  n_left  <- if (k == 1L) 1L else prod(n_grid[seq_len(k - 1L)])
  n_right <- if (k == d)  1L else prod(n_grid[(k + 1L):d])
  n_total <- prod(n_grid)
  # L_uniq: first n_left rows of L[[k]]
  L_uniq <- L_all[[k]][seq_len(n_left), , drop = FALSE]
  # R_uniq: rows at positions seq(1, n_total, by = n_left * n_k)
  stride_r <- n_left * n_grid[k]
  R_uniq <- R_all[[k]][seq(1L, n_total, by = stride_r)[seq_len(n_right)], ,
                        drop = FALSE]
  list(L = L_uniq, R = R_uniq)
}

#' Triple-mode contraction: b_k = (L_uniq, Bk, R_uniq)^T * vec(Y_centered).
#'
#' Computes b_k = X_k' y without forming X_k, using the array structure.
#' Equivalent to `crossprod(X_k, y)` but O(n_left * p + n_k * r_l * n_right)
#' instead of O(n_total * q_k).
#'
#' @param L_uniq n_left x r_l matrix.
#' @param Bk n_k x p matrix (marginal B-spline basis).
#' @param R_uniq n_right x r_r matrix.
#' @param Y_arr d-way array, dimensions (n_1,...,n_d), dim1 fastest (expand.grid).
#' @param k Margin index (determines which mode is Bk vs L vs R).
#' @param n_grid Integer grid size vector.
#' @return Numeric vector of length r_l * p * r_r (= ncol(X_k)).
#' @keywords internal
#' @noRd
.tt_array_rhs <- function(L_uniq, Bk, R_uniq, Y_arr, k, n_grid) {
  d       <- length(n_grid)
  n_left  <- nrow(L_uniq)
  n_k     <- n_grid[k]
  n_right <- nrow(R_uniq)
  rl <- ncol(L_uniq); p <- ncol(Bk)
  Y_flat <- as.numeric(array(as.numeric(Y_arr), c(n_left, n_k, n_right)))
  if (exists("tt_array_rhs_cpp", mode = "function")) {
    return(as.numeric(tt_array_rhs_cpp(
      L_uniq, Bk, R_uniq, Y_flat, as.integer(n_left), as.integer(n_k),
      as.integer(n_right)
    )))
  }
  # Fallback R path (used if C++ not loaded)
  Y_flat_arr <- array(Y_flat, c(n_left, n_k, n_right))
  Y1 <- crossprod(L_uniq, matrix(Y_flat_arr, n_left, n_k * n_right))
  Y1_perm <- aperm(array(Y1, c(rl, n_k, n_right)), c(2L, 1L, 3L))
  BtY     <- crossprod(Bk, matrix(Y1_perm, n_k, rl * n_right))
  BtY_perm <- aperm(array(BtY, c(p, rl, n_right)), c(2L, 1L, 3L))
  as.numeric(matrix(BtY_perm, rl * p, n_right) %*% R_uniq)
}

# --------------------------------------------------------------------
# Public: array-mode Gram + RHS (replaces tt_gram_rhs in array mode)
# --------------------------------------------------------------------

#' Array-mode Gram and RHS for one TT core (Gaussian, no X).
#'
#' Computes S = X_k' X_k and b = X_k' y using the Kronecker structure of
#' a complete data grid, without materialising the n x q_k design matrix X_k.
#'
#' Restriction: unweighted Gaussian.  For weighted / GLM use on a complete
#' grid, see [tt_gram_rhs_array_weighted()] (row-tensor Gram).
#'
#' @param k Margin index (1-based).
#' @param Left Left interface for core k.  When `marginal_iface = TRUE`
#'   this is already the unique matrix (n_left x r_l) from
#'   [left_interfaces_marginal()].  When `marginal_iface = FALSE` it is
#'   the full scattered interface (n_total x r_l) and unique rows are
#'   extracted internally.
#' @param Right Analogous right interface.
#' @param Bk n_k x p marginal B-spline basis for margin k.
#' @param Y_centered d-way array of centred responses, dim1 fastest.
#' @param n_grid Integer vector (n_1, ..., n_d).
#' @param marginal_iface Logical. If TRUE, Left/Right are already the unique
#'   (marginal) interfaces; no row extraction needed.
#' @return List with `S` (q_k x q_k), `b` (length q_k), `q` (= q_k), and
#'   `method = "array_kron"`.
#' @keywords internal
#' @noRd
tt_gram_rhs_array <- function(k, Left, Right, Bk, Y_centered, n_grid,
                               marginal_iface = FALSE) {
  d       <- length(n_grid)
  n_total <- prod(n_grid)
  n_left  <- if (k == 1L) 1L else prod(n_grid[seq_len(k - 1L)])
  n_k     <- n_grid[k]
  n_right <- if (k == d) 1L else prod(n_grid[(k + 1L):d])

  if (isTRUE(marginal_iface)) {
    # Left/Right already are the unique marginal interfaces.
    L_uniq <- Left
    R_uniq <- Right
  } else {
    # Scattered layout: extract unique rows.
    # Left: first n_left rows.
    L_uniq <- Left[seq_len(n_left), , drop = FALSE]
    # Right: rows at stride n_left * n_k.
    stride_r <- n_left * n_k
    idx_r    <- seq(1L, n_total, by = stride_r)
    R_uniq   <- Right[idx_r[seq_len(n_right)], , drop = FALSE]
  }

  # Gram via Kronecker product of three small marginal grams.
  S <- kronecker(crossprod(R_uniq),
                 kronecker(crossprod(Bk), crossprod(L_uniq)))
  # RHS via triple-mode contraction (never forms X_k).
  b <- .tt_array_rhs(L_uniq, Bk, R_uniq, Y_centered, k, n_grid)
  q <- ncol(L_uniq) * ncol(Bk) * ncol(R_uniq)
  list(S = S, b = b, q = q, method = "array_kron")
}

# --------------------------------------------------------------------
# Weighted array-mode Gram (GLAM-style row tensors): X_k' W X_k without X_k.
# --------------------------------------------------------------------

#' Row tensor (row-wise self Khatri--Rao): column (a,a') at index a + (a'-1)*c.
#' @keywords internal
#' @noRd
.tt_row_tensor <- function(M) {
  c_ <- ncol(M)
  M[, rep(seq_len(c_), times = c_), drop = FALSE] *
    M[, rep(seq_len(c_), each = c_), drop = FALSE]
}

#' Weighted array-mode Gram and RHS for one TT core (general weights, no X).
#'
#' Computes S = X_k' W X_k and b = X_k' W z on a complete grid for an
#' arbitrary (non-separable) weight array W, without materialising the
#' n x q_k design matrix.  This is the GLAM row-tensor identity applied to
#' the three factors (L, B_k, R) of the conditional design:
#'
#'   S[(a,j,b),(a',j',b')] = sum_i w_i L[iL,a]L[iL,a'] B[ik,j]B[ik,j'] R[iR,b]R[iR,b']
#'
#' i.e. a triple-mode contraction of the weight array with the row tensors
#' L~ (n_left x r_l^2), B~ (n_k x p^2), R~ (n_right x r_r^2).  Cost is
#' O(m r^2) for the dominant contraction (m = prod(n_grid)) versus
#' O(m q_k^2) for the scattered weighted Gram.  Used by Poisson/GLM PIRLS
#' in array mode, where the working weights change every iteration.
#'
#' @param k Margin index (1-based).
#' @param Left,Right Interfaces: marginal (unique-row) when
#'   `marginal_iface = TRUE`, otherwise full scattered (n x r) matrices from
#'   which unique rows are extracted (expand.grid / dim1-fastest layout).
#' @param Bk n_k x p marginal B-spline basis.
#' @param w Weight vector, length prod(n_grid), grid (dim1-fastest) order.
#' @param z Working response vector, same length/order as `w`.
#' @param n_grid Integer vector (n_1, ..., n_d).
#' @param marginal_iface Logical; see `Left`.
#' @return List with `S` (q_k x q_k), `b` (length q_k), `q`, and
#'   `method = "array_kron_weighted"`.
#' @keywords internal
#' @noRd
tt_gram_rhs_array_weighted <- function(k, Left, Right, Bk, w, z, n_grid,
                                       marginal_iface = FALSE) {
  d       <- length(n_grid)
  n_total <- prod(n_grid)
  n_left  <- if (k == 1L) 1L else prod(n_grid[seq_len(k - 1L)])
  n_k     <- n_grid[k]
  n_right <- if (k == d) 1L else prod(n_grid[(k + 1L):d])

  if (isTRUE(marginal_iface)) {
    L_uniq <- Left
    R_uniq <- Right
  } else {
    L_uniq <- Left[seq_len(n_left), , drop = FALSE]
    idx_r  <- seq(1L, n_total, by = n_left * n_k)
    R_uniq <- Right[idx_r[seq_len(n_right)], , drop = FALSE]
  }
  rl <- ncol(L_uniq); p <- ncol(Bk); rr <- ncol(R_uniq)

  Lt <- .tt_row_tensor(L_uniq)   # n_left  x rl^2
  Bt <- .tt_row_tensor(Bk)       # n_k     x p^2
  Rt <- .tt_row_tensor(R_uniq)   # n_right x rr^2

  # Triple-mode contraction of the weight array with the row tensors.
  W_flat <- array(as.numeric(w), c(n_left, n_k, n_right))
  T1  <- crossprod(Lt, matrix(W_flat, n_left, n_k * n_right))   # rl^2 x (n_k n_right)
  T1p <- aperm(array(T1, c(rl * rl, n_k, n_right)), c(2L, 1L, 3L))
  T2  <- crossprod(Bt, matrix(T1p, n_k, rl * rl * n_right))     # p^2 x (rl^2 n_right)
  T2p <- aperm(array(T2, c(p * p, rl * rl, n_right)), c(2L, 1L, 3L))
  T3  <- matrix(T2p, rl * rl * p * p, n_right) %*% Rt           # (rl^2 p^2) x rr^2

  # Reorder (a,a',j,j',b,b') -> ((a,j,b),(a',j',b')), matching the
  # kron(R, kron(B, L)) column ordering of the conditional design.
  S6 <- array(T3, c(rl, rl, p, p, rr, rr))
  q  <- rl * p * rr
  S  <- matrix(aperm(S6, c(1L, 3L, 5L, 2L, 4L, 6L)), q, q)
  S  <- (S + t(S)) / 2

  # RHS: X_k' W z = unweighted triple contraction of the (w * z) array.
  b <- .tt_array_rhs(L_uniq, Bk, R_uniq,
                     array(as.numeric(w) * as.numeric(z), n_grid), k, n_grid)
  list(S = S, b = b, q = q, method = "array_kron_weighted")
}
