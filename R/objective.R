#' Training penalized objective (same definition for ALS and L-BFGS).
#'
#' @param fit A [ttpspline()] object.
#' @param X Training covariates used to build bases (required).
#' @param y Optional response (defaults to `fit$y`).
#' @return List with `value`, `nll_or_sse`, `penalty`, `rss`, `eta`.
#' @export
tt_objective <- function(fit, X, y = NULL) {
  if (is.null(y)) y <- fit$y
  if (is.null(y)) stop("Need y (or fit$y).", call. = FALSE)
  y <- as.numeric(y)
  X <- as.matrix(X)
  basis <- eval_marginal_bases(X, fit$knots, fit$degree)
  cores <- fit$cores
  intercept <- fit$intercept
  lambda <- as.numeric(fit$lambda)
  d <- length(cores)
  p <- ncol(basis[[1]])
  ranks <- integer(d + 1L)
  ranks[1] <- dim(cores[[1]])[1]
  for (k in seq_len(d)) ranks[k + 1L] <- dim(cores[[k]])[3]
  po <- fit$penalty_order %||% 2L
  penalties <- lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p, ranks[k + 1L], po)
  })
  fam <- fit$family %||% normalize_family(fit$family_key %||% "gaussian")
  out <- tt_glm_penalized_objective(y, cores, intercept, basis, penalties, lambda, fam)
  list(
    value = out$value,
    nll_or_sse = out$nll,
    penalty = out$penalty,
    rss = if (identical(family_key(fam), "gaussian")) 2 * out$nll else NA_real_,
    eta = out$eta
  )
}
