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

# ---- Cached left/right TT penalty environments ----------------------------
# P_k^full = P_left + P_own + P_right with
#   P_own   = λ_k kron(R0, kron(DtD_k, L0))
#   P_left  = kron(R0, kron(I_p, LP))     # sum_{m<k} λ_m P_{k,m}
#   P_right = kron(RP, kron(I_p, L0))     # sum_{m>k} λ_m P_{k,m}

#' Ordinary left transfer: L (r_in x r_in) through core -> r_out x r_out.
#' @keywords internal
#' @noRd
.tt_left_transfer_gram <- function(L, Ct) {
  rl <- dim(Ct)[1L]
  p <- dim(Ct)[2L]
  rr <- dim(Ct)[3L]
  Gn <- matrix(0, rr, rr)
  for (j in seq_len(p)) {
    C <- matrix(Ct[, j, ], rl, rr)
    Gn <- Gn + crossprod(C, L %*% C)
  }
  (Gn + t(Gn)) / 2
}

#' Ordinary right transfer: E (r_out x r_out) through core -> r_in x r_in.
#' @keywords internal
#' @noRd
.tt_right_transfer_gram <- function(E, Ct) {
  rl <- dim(Ct)[1L]
  p <- dim(Ct)[2L]
  rr <- dim(Ct)[3L]
  Gn <- matrix(0, rl, rl)
  for (j in seq_len(p)) {
    C <- matrix(Ct[, j, ], rl, rr)
    Gn <- Gn + C %*% E %*% t(C)
  }
  (Gn + t(Gn)) / 2
}

#' Site penalty contribution on a left transfer (adds λ DtD at this core).
#' @keywords internal
#' @noRd
.tt_left_transfer_pen_site <- function(L0, Ct, DtD, lambda_t) {
  rl <- dim(Ct)[1L]
  p <- dim(Ct)[2L]
  rr <- dim(Ct)[3L]
  if (lambda_t == 0) return(matrix(0, rr, rr))
  Ms <- vector("list", p)
  for (j in seq_len(p)) {
    C <- matrix(Ct[, j, ], rl, rr)
    Ms[[j]] <- L0 %*% C
  }
  Gn <- matrix(0, rr, rr)
  for (j in seq_len(p)) {
    Cj <- matrix(Ct[, j, ], rl, rr)
    for (jp in seq_len(p)) {
      dval <- DtD[j, jp]
      if (dval != 0) {
        Gn <- Gn + (lambda_t * dval) * crossprod(Cj, Ms[[jp]])
      }
    }
  }
  (Gn + t(Gn)) / 2
}

#' Site penalty contribution on a right transfer.
#' @keywords internal
#' @noRd
.tt_right_transfer_pen_site <- function(E0, Ct, DtD, lambda_t) {
  rl <- dim(Ct)[1L]
  p <- dim(Ct)[2L]
  rr <- dim(Ct)[3L]
  if (lambda_t == 0) return(matrix(0, rl, rl))
  # Gn = λ sum_{j,j'} DtD[j,j'] C_j E0 C_{j'}^T
  Cs <- vector("list", p)
  for (j in seq_len(p)) Cs[[j]] <- matrix(Ct[, j, ], rl, rr)
  # Precompute E0 C_j^T
  ECt <- vector("list", p)
  for (j in seq_len(p)) ECt[[j]] <- E0 %*% t(Cs[[j]])
  Gn <- matrix(0, rl, rl)
  for (j in seq_len(p)) {
    for (jp in seq_len(p)) {
      dval <- DtD[j, jp]
      if (dval != 0) {
        Gn <- Gn + (lambda_t * dval) * (Cs[[j]] %*% ECt[[jp]])
      }
    }
  }
  (Gn + t(Gn)) / 2
}

#' Absorb core t into left environments (L0, LP) -> (L0', LP').
#' @keywords internal
#' @noRd
.tt_left_env_absorb <- function(L0, LP, Ct, DtD, lambda_t) {
  if (exists("tt_penalty_left_env_absorb_cpp", mode = "function")) {
    return(tt_penalty_left_env_absorb_cpp(L0, LP, Ct, DtD, lambda_t))
  }
  list(
    L0 = .tt_left_transfer_gram(L0, Ct),
    LP = .tt_left_transfer_gram(LP, Ct) +
      .tt_left_transfer_pen_site(L0, Ct, DtD, lambda_t)
  )
}

#' Absorb core t into right environments (R0, RP) on its right bond -> left bond.
#' @keywords internal
#' @noRd
.tt_right_env_absorb <- function(R0, RP, Ct, DtD, lambda_t) {
  list(
    R0 = .tt_right_transfer_gram(R0, Ct),
    RP = .tt_right_transfer_gram(RP, Ct) +
      .tt_right_transfer_pen_site(R0, Ct, DtD, lambda_t)
  )
}

#' Precompute right ordinary / cumulative-penalty environments for a sweep.
#'
#' `R0[[k]]`, `RP[[k]]` are `r_k x r_k` matrices on the bond to the right of
#' core `k` (contraction of cores `k+1..d`). For the last core, both are `1x1`.
#'
#' @keywords internal
#' @noRd
tt_penalty_prepare_right_envs <- function(cores, lambda, DtD_list) {
  d <- length(cores)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (exists("tt_penalty_prepare_right_envs_cpp", mode = "function")) {
    return(tt_penalty_prepare_right_envs_cpp(cores, lambda, DtD_list))
  }
  R0 <- vector("list", d)
  RP <- vector("list", d)
  # Bond after core d
  cur0 <- matrix(1, 1, 1)
  curP <- matrix(0, 1, 1)
  R0[[d]] <- cur0
  RP[[d]] <- curP
  if (d == 1L) {
    return(list(R0 = R0, RP = RP))
  }
  for (t in d:2) {
    absb <- .tt_right_env_absorb(
      cur0, curP, cores[[t]], DtD_list[[t]], lambda[t]
    )
    cur0 <- absb$R0
    curP <- absb$RP
    R0[[t - 1L]] <- cur0
    RP[[t - 1L]] <- curP
  }
  list(R0 = R0, RP = RP)
}

#' Assemble P_own / P_other / P_full from left/right environments at core k.
#'
#' `P_own` is returned **without** λ_k (same contract as
#' [core_penalty_own_exact()]); callers form `λ_k P_own + P_other`.
#' @keywords internal
#' @noRd
tt_penalty_from_envs <- function(L0, LP, R0, RP, DtD_k, lambda_k, p) {
  if (exists("tt_penalty_from_envs_cpp", mode = "function")) {
    return(tt_penalty_from_envs_cpp(L0, LP, R0, RP, DtD_k, lambda_k, as.integer(p)))
  }
  Ip <- diag(p)
  P_own <- kronecker(R0, kronecker(DtD_k, L0))
  P_left <- kronecker(R0, kronecker(Ip, LP))
  P_right <- kronecker(RP, kronecker(Ip, L0))
  P_other <- P_left + P_right
  P_full <- as.numeric(lambda_k) * P_own + P_other
  list(
    P_own = 0.5 * (P_own + t(P_own)),
    P_other = 0.5 * (P_other + t(P_other)),
    P_full = 0.5 * (P_full + t(P_full)),
    method = "tt_env"
  )
}

#' Exact conditional P_k^full via left/right cumulative penalty environments.
#'
#' Algebraically equivalent to [tt_conditional_penalty_full_tt()] / the dense
#' reference, but avoids the O(m_k^2) unit-core loop.
#'
#' @keywords internal
#' @noRd
tt_conditional_penalty_full_env <- function(cores, k, lambda, DtD_list) {
  d <- length(cores)
  k <- as.integer(k)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (exists("tt_conditional_penalty_full_env_cpp", mode = "function")) {
    return(tt_conditional_penalty_full_env_cpp(cores, k, lambda, DtD_list))
  }
  right <- tt_penalty_prepare_right_envs(cores, lambda, DtD_list)
  L0 <- matrix(1, 1, 1)
  LP <- matrix(0, 1, 1)
  if (k > 1L) {
    for (t in seq_len(k - 1L)) {
      absb <- .tt_left_env_absorb(
        L0, LP, cores[[t]], DtD_list[[t]], lambda[t]
      )
      L0 <- absb$L0
      LP <- absb$LP
    }
  }
  p <- dim(cores[[k]])[2L]
  tt_penalty_from_envs(
    L0, LP, right$R0[[k]], right$RP[[k]],
    DtD_list[[k]], lambda[k], p
  )
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
                                        max_dense = 20000L,
                                        method = c("auto", "env", "tt_cpp", "tt", "dense")) {
  d <- length(cores)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (length(lambda) != d) {
    stop("tt_conditional_penalty_full: lambda length mismatch.", call. = FALSE)
  }
  method <- match.arg(method)
  DtD_list <- tt_DtD_list(cores, penalty_order, cyclic)

  # Default / env: cumulative left–own–right environments (no unit-core loop).
  if (method %in% c("auto", "env")) {
    if (exists("tt_conditional_penalty_full_env_cpp", mode = "function")) {
      out <- tt_conditional_penalty_full_env_cpp(
        cores, as.integer(k), lambda, DtD_list
      )
      out$method <- out$method %||% "tt_env_cpp"
      return(out)
    }
    out <- tt_conditional_penalty_full_env(cores, k, lambda, DtD_list)
    out$method <- "tt_env"
    return(out)
  }

  # Legacy / diagnostic paths
  if (identical(method, "tt_cpp")) {
    if (!exists("tt_conditional_penalty_full_cpp", mode = "function")) {
      stop("tt_conditional_penalty_full_cpp not available.", call. = FALSE)
    }
    out <- tt_conditional_penalty_full_cpp(cores, as.integer(k), lambda, DtD_list)
    out$method <- out$method %||% "tt_cpp"
    return(out)
  }
  if (identical(method, "dense")) {
    p_vec <- vapply(cores, function(z) dim(z)[2L], integer(1))
    nfull <- prod(as.numeric(p_vec))
    mk <- length(cores[[k]])
    if (!is.finite(nfull) || nfull > max_dense || nfull * mk > 5e7) {
      stop("dense path infeasible for these dimensions.", call. = FALSE)
    }
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
                          penalty_mode = "global",
                          linear = NULL, beta = NULL) {
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  eta <- tt_eta(offset, intercept, cores, basis, linear = linear, beta = beta)
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
