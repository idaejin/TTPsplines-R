#' Build / validate a TT rank chain.
#'
#' @param rank Scalar internal rank, or length-`d-1` vector of internal ranks,
#'   or a full length-`d+1` chain starting and ending with 1.
#' @param d Number of modes (covariates).
#' @return Integer vector of length `d + 1` with boundaries 1.
#' @examples
#' tt_rank(3, d = 5)
#' tt_rank(c(2, 4, 3), d = 4)
#' @export
tt_rank <- function(rank, d) {
  d <- as.integer(d)
  stopifnot(d >= 2L)
  r <- as.integer(rank)
  if (length(r) == 1L) {
    return(c(1L, rep(r, d - 1L), 1L))
  }
  if (length(r) == d - 1L) {
    return(c(1L, r, 1L))
  }
  if (length(r) == d + 1L) {
    if (r[1] != 1L || r[length(r)] != 1L) {
      stop("Full rank chain must start and end with 1.")
    }
    return(r)
  }
  stop("rank must be scalar, length d-1, or length d+1.")
}

#' Number of free TT parameters for equal basis size p.
#' @keywords internal
tt_npar <- function(p, ranks) {
  d <- length(ranks) - 1L
  p <- rep(as.integer(p), length.out = d)
  sum(vapply(seq_len(d), function(k) {
    ranks[k] * p[k] * ranks[k + 1L]
  }, numeric(1)))
}

#' Dense tensor coefficient count.
#' @keywords internal
dense_npar <- function(p, d = length(p)) {
  p <- rep(as.integer(p), length.out = d)
  prod(p)
}
