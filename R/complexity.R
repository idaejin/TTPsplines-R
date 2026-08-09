#' TT complexity layers (storage != intrinsic != EDF).
#'
#' Separates concepts that must not be collapsed into a single “complexity”:
#'
#' \describe{
#'   \item{Full tensor coefficients}{\eqn{N_{\mathrm{full}}=\prod_k p_k}.}
#'   \item{TT stored parameters}{\eqn{N_{\mathrm{TT}}=\sum_k r_{k-1}p_k r_k}
#'     (entries stored in the cores; \eqn{r_0=r_d=1}).}
#'   \item{TT intrinsic dimension}{\eqn{\dim(\mathcal{M}_r)=N_{\mathrm{TT}}-\sum_{k=1}^{d-1}r_k^2}
#'     after removing TT gauge redundancy (regular case).}
#'   \item{Compression ratio}{\eqn{CR=N_{\mathrm{full}}/N_{\mathrm{TT}}}.}
#'   \item{EDF}{Joint linearized effective degrees of freedom after penalization;
#'     \strong{not} equal to \eqn{N_{\mathrm{TT}}}.}
#' }
#'
#' Methodological reading: rank \eqn{r} is \emph{structural}
#' capacity; \eqn{\boldsymbol{\lambda}} / EDF is \emph{smoothness} complexity.
#'
#' @param fit A [ttpspline()] object, or a list with `d`, `k` (or `basis_dim` /
#'   `p`), and `rank` (full chain or as accepted by [tt_rank()]).
#' @param p Optional basis sizes (scalar or length-`d`); overrides `fit$k`.
#' @param rank Optional rank (scalar / internal / full chain).
#' @param d Optional dimension (with `p` and `rank` allows calling without a fit).
#' @param edf Optional EDF to attach (defaults to `fit$edf` when present).
#'
#' @return An object of class `"tt_complexity"` with the fields above
#'   (aliases `npar_dense` / `npar_tt` kept for compatibility).
#'
#' @examples
#' # Spectacular storage compression (no fit needed)
#' tt_complexity(d = 8, p = 10, rank = 3)
#'
#' @export
tt_complexity <- function(fit = NULL, p = NULL, rank = NULL, d = NULL,
                          edf = NULL) {
  if (!is.null(fit)) {
    if (is.null(d)) d <- fit$d
    if (is.null(p)) {
      p <- fit$k %||% fit$basis_dim %||% fit$p
    }
    if (is.null(rank)) rank <- fit$rank
    if (is.null(edf) && !is.null(fit$edf)) edf <- fit$edf
  }
  if (is.null(d) || is.null(p) || is.null(rank)) {
    stop("Provide a fit, or d + p + rank.", call. = FALSE)
  }
  d <- as.integer(d)
  ranks <- if (length(as.integer(rank)) == d + 1L) {
    as.integer(rank)
  } else {
    tt_rank(rank, d = d)
  }
  p_vec <- rep(as.integer(p), length.out = d)

  n_full <- as.numeric(prod(p_vec))
  n_tt <- tt_npar(p_vec, ranks)
  n_gauge <- tt_gauge_dim(ranks)
  n_intr <- n_tt - n_gauge
  cr <- n_full / max(n_tt, 1)

  internal <- ranks[-c(1L, length(ranks))]
  rank_uniform <- if (length(unique(internal)) == 1L) as.integer(internal[1]) else NA_integer_

  out <- list(
    # geometry
    d = d,
    basis_dim = if (length(unique(p_vec)) == 1L) p_vec[1] else p_vec,
    p = p_vec,
    rank = rank_uniform,
    rank_chain = ranks,
    rank_internal = internal,
    # 1–2 storage / compression
    n_full = n_full,
    n_tt_stored = as.numeric(n_tt),
    compression_ratio = cr,
    # 4 gauge / intrinsic
    n_gauge = as.numeric(n_gauge),
    n_tt_intrinsic = as.numeric(n_intr),
    # 3 statistical (optional)
    edf = if (is.null(edf)) NA_real_ else as.numeric(edf),
    # aliases (backward compatible)
    npar_dense = n_full,
    npar_tt = as.numeric(n_tt),
    npar_tt_intrinsic = as.numeric(n_intr)
  )
  class(out) <- "tt_complexity"
  out
}

#' Gauge redundancy dimension \eqn{\sum_{k=1}^{d-1} r_k^2}.
#' @keywords internal
tt_gauge_dim <- function(ranks) {
  ranks <- as.integer(ranks)
  d <- length(ranks) - 1L
  if (d < 2L) return(0)
  # interfaces 1..d-1 (exclude boundary ranks r_0=r_d=1)
  sum(as.numeric(ranks[2:d])^2)
}

#' @export
print.tt_complexity <- function(x, ...) {
  cat("TTPsplines complexity layers\n")
  cat("(storage != intrinsic TT dimension != EDF)\n\n")

  cat("--- Geometry ---\n")
  cat(sprintf("  d:                         %d\n", x$d))
  bd <- x$basis_dim
  if (length(bd) == 1L) {
    cat(sprintf("  Basis size p (per margin): %d\n", bd))
  } else {
    cat(sprintf("  Basis sizes p_k:           %s\n", paste(bd, collapse = ",")))
  }
  if (is.finite(x$rank)) {
    cat(sprintf("  TT rank (uniform r):       %d\n", x$rank))
  } else {
    cat(sprintf("  TT ranks (internal):       %s\n",
                paste(x$rank_internal, collapse = "-")))
  }
  cat(sprintf("  Rank chain:                %s\n", paste(x$rank_chain, collapse = "-")))

  cat("\n--- 1. Parametric / storage ---\n")
  cat(sprintf("  Full tensor coefficients:  %s\n",
              format(x$n_full, big.mark = ",", scientific = FALSE)))
  cat(sprintf("  TT stored parameters:      %s\n",
              format(x$n_tt_stored, big.mark = ",", scientific = FALSE)))
  cat(sprintf("  Compression ratio:         %.2fx\n", x$compression_ratio))

  cat("\n--- 2. TT manifold (gauge) ---\n")
  cat(sprintf("  Gauge redundancy:          %s\n",
              format(x$n_gauge, big.mark = ",", scientific = FALSE)))
  cat(sprintf("  TT intrinsic dimension:    %s\n",
              format(x$n_tt_intrinsic, big.mark = ",", scientific = FALSE)))
  cat("  (intrinsic = stored parameters - gauge redundancy)\n")

  cat("\n--- 3. Statistical effective complexity ---\n")
  if (is.finite(x$edf)) {
    cat(sprintf("  Joint EDF:                 %.2f\n", x$edf))
    cat(sprintf("  EDF / stored parameters:   %.2f\n",
                x$edf / max(x$n_tt_stored, 1)))
    cat(sprintf("  EDF / intrinsic TT dim.:   %.2f\n",
                x$edf / max(x$n_tt_intrinsic, 1)))
  } else {
    cat("  Joint EDF:                 NA (fit with compute_edf=TRUE, or pass edf=)\n")
  }

  cat("\nInterpretation:\n")
  cat("  rank   = structural capacity\n")
  cat("  lambda = directional smoothness\n")
  cat("  EDF    = effective fitted flexibility\n")
  invisible(x)
}
