#' Training penalized objective (same definition for ALS and L-BFGS).
#'
#' Gaussian:
#' \deqn{\tfrac12\|y-\alpha-\hat f\|^2 + \tfrac12\sum_k \lambda_k\, g_k^\top P_k g_k.}
#'
#' @param fit A [ttpspline()] object.
#' @param X Training covariates used to build bases (required).
#' @param y Optional response (defaults to `fit$y`).
#' @return List with `value`, `nll_or_sse`, `penalty`, `rss`.
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
  f <- tt_contraction(cores, basis)
  eta <- intercept + f
  key <- fit$family_key %||% "gaussian"
  if (identical(key, "gaussian")) {
    rss <- sum((y - eta)^2)
    nll <- 0.5 * rss
  } else {
    fam <- fit$family
    mu <- invlink_eta(fam, eta)
    nll <- switch(
      key,
      poisson = {
        mu <- pmax(mu, 1e-12)
        sum(mu - y * log(mu))
      },
      bernoulli = {
        mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
        -sum(y * log(mu) + (1 - y) * log(1 - mu))
      },
      stop("Unsupported family for tt_objective.", call. = FALSE)
    )
    rss <- NA_real_
  }
  pen <- .tt_penalty_value_grad(cores, penalties, lambda)$value
  list(
    value = nll + pen,
    nll_or_sse = nll,
    penalty = pen,
    rss = rss,
    eta = eta
  )
}
