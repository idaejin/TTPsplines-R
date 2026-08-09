#' Coefficient compression diagnostics for a TT P-spline fit.
#'
#' @param fit A [ttpspline()] object (or list with `d`, `k`, `rank`, …).
#' @return A list / printable object with dense vs TT parameter counts.
#' @export
tt_complexity <- function(fit) {
  d <- fit$d
  p <- fit$k
  ranks <- fit$rank
  npar_tt <- fit$npar_tt
  npar_dense <- fit$npar_dense
  out <- list(
    d = d,
    basis_dim = p,
    rank_chain = ranks,
    npar_dense = npar_dense,
    npar_tt = npar_tt,
    compression_ratio = npar_dense / max(npar_tt, 1)
  )
  class(out) <- "tt_complexity"
  out
}

#' @export
print.tt_complexity <- function(x, ...) {
  cat("TTPsplines complexity\n")
  cat(sprintf("  Dimensions d:              %d\n", x$d))
  cat(sprintf("  Basis size (per margin):   %d\n", x$basis_dim))
  cat(sprintf("  Rank chain:                %s\n", paste(x$rank_chain, collapse = "-")))
  cat(sprintf("  Full tensor coefficients:  %s\n", format(x$npar_dense, big.mark = ",")))
  cat(sprintf("  TT parameters:             %s\n", format(x$npar_tt, big.mark = ",")))
  cat(sprintf("  Compression ratio:         %.2f\n", x$compression_ratio))
  invisible(x)
}
