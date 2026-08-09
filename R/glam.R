# Minimal GLAM helpers for Paper 1 fixed-λ grid benchmarks (not the main API).

#' @keywords internal
glam_H <- function(M, A) {
  d <- dim(A)
  Amat <- matrix(A, nrow = d[1])
  array(M %*% Amat, c(nrow(M), d[-1]))
}

#' @keywords internal
glam_Rotate <- function(A) {
  d <- dim(A)
  if (is.null(d) || length(d) == 1L) return(A)
  aperm(A, c(seq_along(d)[-1L], 1L))
}

#' @keywords internal
glam_RH <- function(M, A) glam_Rotate(glam_H(M, A))

#' @keywords internal
glam_linear_predictor <- function(Theta, B_list) {
  out <- Theta
  for (k in seq_along(B_list)) out <- glam_RH(B_list[[k]], out)
  out
}

#' @keywords internal
glam_Bt_y <- function(Y, B_list) {
  out <- Y
  for (k in seq_along(B_list)) out <- glam_RH(t(B_list[[k]]), out)
  out
}

#' @keywords internal
glam_penalty <- function(p_vec, lambda, penalty_order = 2L) {
  d <- length(p_vec)
  lambda <- rep(as.numeric(lambda), length.out = d)
  npar <- prod(p_vec)
  P <- matrix(0, npar, npar)
  for (k in seq_len(d)) {
    Dd <- difference_penalty(p_vec[k], penalty_order)
    blocks <- vector("list", d)
    for (j in seq_len(d)) {
      blocks[[j]] <- if (j == k) Dd else diag(p_vec[j])
    }
    K <- blocks[[1L]]
    for (j in 2:d) K <- kronecker(blocks[[j]], K)
    P <- P + lambda[k] * K
  }
  P
}

#' @keywords internal
glam_gram_unweighted <- function(B_list) {
  G <- crossprod(B_list[[1L]])
  if (length(B_list) == 1L) return(G)
  for (k in 2:length(B_list)) G <- kronecker(crossprod(B_list[[k]]), G)
  G
}

#' Fixed-λ Gaussian GLAM on a d-way grid (Currie array methods).
#'
#' Exposed for Paper 1 compression benchmarks; scattered-data users should
#' prefer [ttpspline()].
#'
#' @param Y d-way array.
#' @param B_list Marginal bases (`n_k × p_k`).
#' @param lambda Fixed smoothing.
#' @export
glam_fit_gaussian <- function(Y, B_list, lambda = 1, penalty_order = 2L) {
  d <- length(B_list)
  stopifnot(length(dim(Y)) == d)
  p_vec <- vapply(B_list, ncol, integer(1))
  lambda <- rep(as.numeric(lambda), length.out = d)
  intercept <- mean(Y)
  Yc <- Y - intercept
  XtX <- glam_gram_unweighted(B_list)
  P <- glam_penalty(p_vec, lambda, penalty_order)
  ridge <- ridge_scale(XtX)
  coef <- drop(solve_spd(XtX + P + ridge * diag(nrow(XtX)),
                         as.numeric(glam_Bt_y(Yc, B_list))))
  Theta <- array(coef, p_vec)
  mu <- intercept + glam_linear_predictor(Theta, B_list)
  list(
    Theta = Theta, intercept = intercept, mu = mu, lambda = lambda,
    npar = prod(p_vec), family = "gaussian", method = "GLAM-fixed"
  )
}
