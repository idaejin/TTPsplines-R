#' Weighted Gram \(X_k' W X_k\) and RHS \(X_k' W z\) from TT interfaces.
#'
#' Does not return the design matrix. See [tt_gram_rhs_cpp()] methods:
#' `blas`, `fused`, `fused_blocked`, `kron`.
#'
#' When `array_data` is supplied (a list with `k`, `L_all`, `R_all`,
#' `Y_centered`, `n_grid`), the Gram and RHS are computed via the array
#' Kronecker trick ([tt_gram_rhs_array()]) without forming the design matrix.
#' This is only valid for unweighted Gaussian data on a complete grid.
#'
#' @param n_threads OpenMP threads for `fused_blocked` observation reduction.
#' @param array_data Optional list for array mode; see above.
#' @keywords internal
#' @noRd
tt_gram_rhs <- function(Left, Right, Bk, z, weight = NULL,
                        method = c("fused_blocked", "blas", "fused", "kron"),
                        block_size = 64L,
                        n_threads = 1L,
                        array_data = NULL) {
  # Array mode: use Kronecker trick (unweighted Gaussian on complete grid)
  if (!is.null(array_data) && is.null(weight)) {
    return(tt_gram_rhs_array(
      k          = array_data$k,
      L_all      = array_data$L_all,
      R_all      = array_data$R_all,
      Bk         = Bk,
      Y_centered = array_data$Y_centered,
      n_grid     = array_data$n_grid
    ))
  }
  method <- match.arg(method)
  if (exists("tt_gram_rhs_cpp", envir = asNamespace("TTPsplines"), inherits = FALSE)) {
    return(tt_gram_rhs_cpp(
      Left, Right, Bk, as.numeric(z),
      weight = if (is.null(weight)) NULL else as.numeric(weight),
      method = method,
      block_size = as.integer(block_size),
      n_threads = as.integer(n_threads)
    ))
  }
  X <- tt_design_core(Left, Right, Bk)
  if (is.null(weight)) {
    list(S = crossprod(X), b = as.numeric(crossprod(X, z)),
         q = ncol(X), method = "blas_R", n_threads = 1L, omp = FALSE)
  } else {
    sw <- sqrt(pmax(as.numeric(weight), 0))
    Xw <- X * sw
    list(S = crossprod(Xw), b = as.numeric(crossprod(Xw, sw * z)),
         q = ncol(X), method = "blas_R", n_threads = 1L, omp = FALSE)
  }
}
