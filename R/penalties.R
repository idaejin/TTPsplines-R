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
#'
#' This is the *separable* ALS surrogate (own-margin only). For the exact
#' global discrete P-spline penalty use [tt_global_penalty_value()].
#'
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

#' Frobenius inner product of two TT coefficient tensors (equal dimensions).
#' @keywords internal
#' @noRd
tt_frobenius_inner <- function(cores_a, cores_b) {
  d <- length(cores_a)
  if (d != length(cores_b)) {
    stop("tt_frobenius_inner: core length mismatch.", call. = FALSE)
  }
  cur <- matrix(1, 1, 1)
  for (k in seq_len(d)) {
    Ga <- cores_a[[k]]
    Gb <- cores_b[[k]]
    rl <- dim(Ga)[1L]
    p <- dim(Ga)[2L]
    rr <- dim(Ga)[3L]
    if (!identical(as.integer(dim(Gb)), as.integer(dim(Ga)))) {
      stop("tt_frobenius_inner: core ", k, " dim mismatch.", call. = FALSE)
    }
    nxt <- matrix(0, rr, dim(Gb)[3L])
    for (j in seq_len(p)) {
      Ca <- matrix(Ga[, j, ], rl, rr)
      Cb <- matrix(Gb[, j, ], dim(Gb)[1L], dim(Gb)[3L])
      nxt <- nxt + crossprod(Ca, cur %*% Cb)
    }
    cur <- nxt
  }
  as.numeric(cur)
}

#' Apply D'D on the physical mode of TT core `m` (self-adjoint mode map).
#' @keywords internal
#' @noRd
tt_apply_DtD_mode <- function(cores, m, DtD) {
  out <- lapply(cores, function(g) array(as.numeric(g), dim = dim(g)))
  Gm <- out[[m]]
  p <- dim(Gm)[2L]
  if (nrow(DtD) != p || ncol(DtD) != p) {
    stop("tt_apply_DtD_mode: DtD size mismatch.", call. = FALSE)
  }
  rl <- dim(Gm)[1L]
  rr <- dim(Gm)[3L]
  newG <- array(0, dim(Gm))
  for (a in seq_len(rl)) {
    for (b in seq_len(rr)) {
      newG[a, , b] <- as.numeric(DtD %*% Gm[a, , b])
    }
  }
  out[[m]] <- newG
  out
}

#' Left bond Gram after contracting cores 1..upto (coefficient space).
#' @keywords internal
#' @noRd
tt_left_bond_gram <- function(cores, upto) {
  upto <- as.integer(upto)
  G <- matrix(1, 1, 1)
  if (upto < 1L) return(G)
  for (t in seq_len(upto)) {
    Ct <- cores[[t]]
    rl <- dim(Ct)[1L]
    p <- dim(Ct)[2L]
    rr <- dim(Ct)[3L]
    Gn <- matrix(0, rr, rr)
    for (j in seq_len(p)) {
      C <- matrix(Ct[, j, ], rl, rr)
      Gn <- Gn + crossprod(C, G %*% C)
    }
    G <- Gn
  }
  G
}

#' Right bond Gram after contracting cores from..d (coefficient space).
#' @keywords internal
#' @noRd
tt_right_bond_gram <- function(cores, from) {
  d <- length(cores)
  from <- as.integer(from)
  G <- matrix(1, 1, 1)
  if (from > d) return(G)
  for (t in d:from) {
    Ct <- cores[[t]]
    rl <- dim(Ct)[1L]
    p <- dim(Ct)[2L]
    rr <- dim(Ct)[3L]
    Gn <- matrix(0, rl, rl)
    for (j in seq_len(p)) {
      C <- matrix(Ct[, j, ], rl, rr)
      Gn <- Gn + C %*% G %*% t(C)
    }
    G <- Gn
  }
  G
}

#' Exact own-margin restriction A_k^T S_k A_k via left/right bond Grams.
#'
#' Equals kronecker(W_R, kronecker(DtD, W_L)) with the same vec order as
#' [core_penalty()] / [tt_design_core()] (a fast, j, b slow).
#'
#' @keywords internal
#' @noRd
core_penalty_own_exact <- function(cores, k, DtD) {
  d <- length(cores)
  k <- as.integer(k)
  rl <- dim(cores[[k]])[1L]
  p <- dim(cores[[k]])[2L]
  rr <- dim(cores[[k]])[3L]
  if (nrow(DtD) != p || ncol(DtD) != p) {
    stop("core_penalty_own_exact: DtD size mismatch.", call. = FALSE)
  }
  W_L <- tt_left_bond_gram(cores, k - 1L)
  W_R <- tt_right_bond_gram(cores, k + 1L)
  if (!isTRUE(all.equal(dim(W_L), c(rl, rl)))) {
    stop("core_penalty_own_exact: W_L dim mismatch.", call. = FALSE)
  }
  if (!isTRUE(all.equal(dim(W_R), c(rr, rr)))) {
    stop("core_penalty_own_exact: W_R dim mismatch.", call. = FALSE)
  }
  # When interfaces are orthonormal, W_L=I, W_R=I recovers core_penalty().
  kronecker(W_R, kronecker(DtD, W_L))
}

#' Per-margin D'D list matching TT physical sizes.
#' @keywords internal
#' @noRd
tt_DtD_list <- function(cores, penalty_order = 2L, cyclic = NULL) {
  d <- length(cores)
  cyclic <- normalize_cyclic(cyclic, d)
  lapply(seq_len(d), function(m) {
    pm <- dim(cores[[m]])[2L]
    if (isTRUE(cyclic[m])) {
      circular_difference_penalty(pm, penalty_order)
    } else {
      difference_penalty(pm, penalty_order)
    }
  })
}

#' Apply mode-m D'D to a vectorized dense coefficient array.
#' @keywords internal
#' @noRd
.apply_DtD_vec <- function(v, p_vec, m, DtD) {
  arr <- array(v, dim = p_vec)
  d <- length(p_vec)
  # Move mode m to front, apply DtD, restore
  perm <- c(m, setdiff(seq_len(d), m))
  a <- aperm(arr, perm)
  dims <- dim(a)
  mat <- matrix(a, nrow = dims[1L])
  mat2 <- DtD %*% mat
  a2 <- array(mat2, c(nrow(DtD), dims[-1L]))
  # DtD is square p x p for ordinary differences D'D
  invperm <- match(seq_len(d), perm)
  as.numeric(aperm(a2, invperm))
}

#' Dense A_k^T S_λ A_k when prod(p) is modest.
#' @keywords internal
#' @noRd
tt_conditional_penalty_full_dense <- function(cores, k, lambda, DtD_list) {
  d <- length(cores)
  k <- as.integer(k)
  rl <- dim(cores[[k]])[1L]
  p <- dim(cores[[k]])[2L]
  rr <- dim(cores[[k]])[3L]
  mk <- rl * p * rr
  p_vec <- vapply(cores, function(z) dim(z)[2L], integer(1))
  nfull <- as.numeric(prod(p_vec))
  A <- matrix(0, nfull, mk)
  for (i in seq_len(mk)) {
    g <- numeric(mk)
    g[i] <- 1
    cc <- cores
    cc[[k]] <- array(g, c(rl, p, rr))
    A[, i] <- as.numeric(tt_full_theta(cc))
  }
  SA_own <- matrix(0, nfull, mk)
  SA_other <- matrix(0, nfull, mk)
  for (j in seq_len(mk)) {
    v <- A[, j]
    SA_own[, j] <- .apply_DtD_vec(v, p_vec, k, DtD_list[[k]])
    so <- numeric(nfull)
    for (m in seq_len(d)) {
      if (m == k) next
      so <- so + lambda[m] * .apply_DtD_vec(v, p_vec, m, DtD_list[[m]])
    }
    SA_other[, j] <- so
  }
  P_own <- crossprod(A, SA_own)
  P_other <- crossprod(A, SA_other)
  P_own <- 0.5 * (P_own + t(P_own))
  P_other <- 0.5 * (P_other + t(P_other))
  list(
    P_own = P_own,
    P_other = P_other,
    P_full = lambda[k] * P_own + P_other
  )
}

#' TT-inner construction of A_k^T S_λ A_k (no k^d materialization).
#' @keywords internal
#' @noRd
tt_conditional_penalty_full_tt <- function(cores, k, lambda, DtD_list) {
  d <- length(cores)
  k <- as.integer(k)
  rl <- dim(cores[[k]])[1L]
  p <- dim(cores[[k]])[2L]
  rr <- dim(cores[[k]])[3L]
  mk <- rl * p * rr

  # Exact own-margin via Grams (fast).
  P_own <- core_penalty_own_exact(cores, k, DtD_list[[k]])
  P_other <- matrix(0, mk, mk)

  if (d == 1L) {
    return(list(
      P_own = P_own,
      P_other = P_other,
      P_full = lambda[k] * P_own + P_other
    ))
  }

  unit_core <- function(i) {
    g <- numeric(mk)
    g[i] <- 1
    array(g, c(rl, p, rr))
  }

  # Precompute 𝒯_m(Θ(e_j)) cores for m ≠ k
  other_m <- setdiff(seq_len(d), k)
  Tmj <- vector("list", mk)
  for (j in seq_len(mk)) {
    cc <- cores
    cc[[k]] <- unit_core(j)
    Tmj[[j]] <- lapply(other_m, function(m) {
      tt_apply_DtD_mode(cc, m, DtD_list[[m]])
    })
  }

  for (i in seq_len(mk)) {
    cores_i <- cores
    cores_i[[k]] <- unit_core(i)
    for (j in seq_len(i)) {
      s <- 0
      for (ii in seq_along(other_m)) {
        m <- other_m[ii]
        s <- s + lambda[m] * tt_frobenius_inner(cores_i, Tmj[[j]][[ii]])
      }
      P_other[i, j] <- s
      P_other[j, i] <- s
    }
  }
  P_other <- 0.5 * (P_other + t(P_other))
  list(
    P_own = P_own,
    P_other = P_other,
    P_full = lambda[k] * P_own + P_other
  )
}

#' Exact conditional restriction of the global discrete P-spline penalty.
#'
#' Returns matrices such that
#' \eqn{g^\top P_{\mathrm{full}} g = J_{\boldsymbol\lambda}(\Theta(G_k))}
#' with other cores fixed, and
#' \eqn{P_{\mathrm{full}} = \lambda_k P_{\mathrm{own}} + P_{\mathrm{other}}}.
#'
#' @param cores Current TT cores.
#' @param k Margin / core index (1..d).
#' @param lambda Length-`d` (or scalar) smoothing vector.
#' @param penalty_order Difference order.
#' @param cyclic Optional cyclic flags.
#' @param max_dense Maximum \code{prod(p)} for the dense validation path.
#' @return List with `P_own`, `P_other`, `P_full`, `method`.
#' @keywords internal
#' @noRd
tt_conditional_penalty_full <- function(cores, k, lambda,
                                        penalty_order = 2L,
                                        cyclic = NULL,
                                        max_dense = 20000L) {
  d <- length(cores)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (length(lambda) != d) {
    stop("tt_conditional_penalty_full: lambda length mismatch.", call. = FALSE)
  }
  DtD_list <- tt_DtD_list(cores, penalty_order, cyclic)
  if (exists("tt_conditional_penalty_full_cpp", mode = "function")) {
    out <- tt_conditional_penalty_full_cpp(cores, as.integer(k), lambda, DtD_list)
    out$method <- out$method %||% "tt_cpp"
    return(out)
  }
  p_vec <- vapply(cores, function(z) dim(z)[2L], integer(1))
  nfull <- prod(as.numeric(p_vec))
  mk <- length(cores[[k]])
  if (is.finite(nfull) && nfull <= max_dense && nfull * mk <= 5e7) {
    out <- tt_conditional_penalty_full_dense(cores, k, lambda, DtD_list)
    out$method <- "dense"
    return(out)
  }
  out <- tt_conditional_penalty_full_tt(cores, k, lambda, DtD_list)
  out$method <- "tt"
  out
}

#' Global discrete P-spline penalty 0.5 * vec(Θ)' S_λ vec(Θ) via TT.
#' @keywords internal
#' @noRd
tt_global_penalty_value <- function(cores, lambda, penalty_order = 2L,
                                    cyclic = NULL) {
  d <- length(cores)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (length(lambda) != d) {
    stop("tt_global_penalty_value: lambda length mismatch.", call. = FALSE)
  }
  DtD_list <- tt_DtD_list(cores, penalty_order, cyclic)
  if (exists("tt_global_penalty_value_cpp", mode = "function")) {
    return(tt_global_penalty_value_cpp(cores, lambda, DtD_list))
  }
  val <- 0
  for (m in seq_len(d)) {
    val <- val + lambda[m] * tt_frobenius_inner(
      cores, tt_apply_DtD_mode(cores, m, DtD_list[[m]])
    )
  }
  0.5 * val
}

#' Alias matching the manuscript API: P_k^full = A_k^T S_lambda A_k.
#' @keywords internal
#' @noRd
tt_core_penalty_full <- function(...) tt_conditional_penalty_full(...)

#' Gaussian ALS criterion Q = 0.5 * RSS + J_lambda(Theta).
#'
#' Package convention uses the 1/2 factors so that the normal equations
#' (X_k'X_k + P_k^full) g = X_k' y_c are exact critical-point conditions.
#' Minimizing Q is equivalent to minimizing RSS + sum_m lambda_m ||Theta x_m Delta||_F^2.
#'
#' @return List with `rss`, `penalty` (J), `value` (Q), `eta`.
#' @keywords internal
#' @noRd
tt_gaussian_Q <- function(y, cores, intercept, basis, lambda,
                          offset = NULL, weights = NULL,
                          penalty_order = 2L, cyclic = NULL,
                          penalty_mode = "global") {
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  eta <- tt_eta(offset, intercept, cores, basis)
  rss <- sum(w * (y - eta)^2)
  normalize_penalty_mode(penalty_mode)
  pen <- tt_global_penalty_value(
    cores, lambda, penalty_order = penalty_order,
    cyclic = cyclic %||% attr(basis, "cyclic")
  )
  list(rss = rss, penalty = pen, value = 0.5 * rss + pen, eta = eta)
}

#' Value and per-core gradients of the global discrete P-spline penalty.
#'
#' Uses \eqn{\nabla_{g_k} J = P_k^{\mathrm{full}} g_k} at the current cores.
#'
#' @keywords internal
#' @noRd
.tt_penalty_value_grad_full <- function(cores, lambda, penalty_order = 2L,
                                        cyclic = NULL) {
  d <- length(cores)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  val <- tt_global_penalty_value(cores, lambda, penalty_order, cyclic)
  grads <- vector("list", d)
  for (k in seq_len(d)) {
    pen_k <- tt_conditional_penalty_full(
      cores, k, lambda, penalty_order = penalty_order, cyclic = cyclic
    )
    g <- as.numeric(cores[[k]])
    grads[[k]] <- array(as.numeric(pen_k$P_full %*% g), dim(cores[[k]]))
  }
  list(value = val, grads = grads)
}
