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
  n_left  <- nrow(L_uniq);  rl <- ncol(L_uniq)
  n_k     <- n_grid[k];     p  <- ncol(Bk)
  n_right <- nrow(R_uniq);  rr <- ncol(R_uniq)
  # Reshape Y to (n_left x n_k x n_right) -- dim1 fastest
  Y_flat <- array(as.numeric(Y_arr), c(n_left, n_k, n_right))
  # Step 1: contract left mode with L_uniq' -> r_l x n_k x n_right
  Y1 <- crossprod(L_uniq, matrix(Y_flat, n_left, n_k * n_right))
  # Step 2: contract middle mode with Bk' -> p x (r_l * n_right)
  # Permute Y1 to (n_k x r_l x n_right), flatten to (n_k x r_l*n_right)
  Y1_perm <- aperm(array(Y1, c(rl, n_k, n_right)), c(2L, 1L, 3L))
  BtY     <- crossprod(Bk, matrix(Y1_perm, n_k, rl * n_right))   # p x (rl*n_right)
  # Step 3: contract right mode with R_uniq -> (r_l*p) x r_r
  # Permute BtY to (r_l x p x n_right), flatten to (r_l*p x n_right)
  BtY_perm <- aperm(array(BtY, c(p, rl, n_right)), c(2L, 1L, 3L))
  b <- as.numeric(matrix(BtY_perm, rl * p, n_right) %*% R_uniq)   # r_l*p*r_r
  b
}

# --------------------------------------------------------------------
# Public: array-mode Gram + RHS (replaces tt_gram_rhs in array mode)
# --------------------------------------------------------------------

#' Array-mode Gram and RHS for one TT core (Gaussian, no X).
#'
#' Computes S = X_k' X_k and b = X_k' y using the Kronecker structure of
#' a complete data grid, without materialising the n x q_k design matrix X_k.
#'
#' Restriction: unweighted Gaussian only in this version.
#' For weighted / GLM use, fall back to [tt_gram_rhs()].
#'
#' @param k Margin index (1-based).
#' @param Left Full left interface matrix for core k (n_total x r_l) as provided
#'   by ALS (scattered layout, dim1 fastest).
#' @param Right Full right interface matrix for core k (n_total x r_r) as
#'   provided by ALS (scattered layout, dim1 fastest).
#' @param Bk n_k x p marginal B-spline basis for margin k.
#' @param Y_centered d-way array of centred responses, dim1 fastest.
#' @param n_grid Integer vector (n_1, ..., n_d).
#' @return List with `S` (q_k x q_k), `b` (length q_k), `q` (= q_k), and
#'   `method = "array_kron"`.
#' @keywords internal
#' @noRd
tt_gram_rhs_array <- function(k, Left, Right, Bk, Y_centered, n_grid) {
  d <- length(n_grid)
  n_total <- prod(n_grid)
  n_left  <- if (k == 1L) 1L else prod(n_grid[seq_len(k - 1L)])
  n_k     <- n_grid[k]
  n_right <- if (k == d) 1L else prod(n_grid[(k + 1L):d])

  # Unique left rows are exactly the first n_left rows in the scattered layout.
  L_uniq <- Left[seq_len(n_left), , drop = FALSE]

  # Unique right rows are at a fixed stride n_left * n_k in the scattered layout.
  stride_r <- n_left * n_k
  idx_r <- seq(1L, n_total, by = stride_r)
  R_uniq <- Right[idx_r[seq_len(n_right)], , drop = FALSE]

  # Gram via Kronecker product of marginal grams
  S <- kronecker(crossprod(R_uniq),
                 kronecker(crossprod(Bk), crossprod(L_uniq)))
  # RHS via triple-mode contraction
  b <- .tt_array_rhs(L_uniq, Bk, R_uniq, Y_centered, k, n_grid)
  q <- ncol(L_uniq) * ncol(Bk) * ncol(R_uniq)
  list(S = S, b = b, q = q, method = "array_kron")
}
