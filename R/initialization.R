#' Initialize TT cores (reproducible across optimizers).
#'
#' Use the same `init` object for ALS / L-BFGS / Adam so optimizer benchmarks
#' differ only by algorithm, not by random start (gauge aside).
#'
#' @param d Number of modes, or a covariate matrix `X` (uses `ncol(X)`).
#' @param rank Scalar or rank specification; see [tt_rank()].
#' @param k Basis size per margin (equal margins in v0).
#' @param seed RNG seed.
#' @param sd Initialization standard deviation.
#' @param ... Ignored (reserved).
#' @return List of TT core arrays with attribute `ranks`.
#' @export
tt_initialize <- function(d, rank = 3, k = 10, seed = 1, sd = 0.15, ...) {
  if (is.matrix(d) || is.data.frame(d)) {
    d <- ncol(as.matrix(d))
  }
  d <- as.integer(d)
  ranks <- tt_rank(rank, d = d)
  cores <- initialize_tt_cores(k, ranks, seed = seed, sd = sd)
  attr(cores, "ranks") <- ranks
  attr(cores, "k") <- as.integer(k)
  attr(cores, "seed") <- as.integer(seed)
  cores
}
