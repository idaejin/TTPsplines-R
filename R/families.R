#' Normalize a GLM family argument to a stats family object.
#' @keywords internal
normalize_family <- function(family) {
  if (inherits(family, "family")) return(family)
  if (is.character(family)) {
    family <- tolower(family)
    return(switch(
      family,
      gaussian = stats::gaussian(),
      poisson = stats::poisson(),
      binomial = stats::binomial(),
      bernoulli = stats::binomial(),
      stop("Unknown family: ", family, call. = FALSE)
    ))
  }
  stop("family must be a stats::family or character string.", call. = FALSE)
}

family_key <- function(family) {
  fam <- normalize_family(family)
  nm <- fam$family
  if (identical(nm, "binomial")) return("bernoulli")
  nm
}

#' Working weights / response for PIRLS (canonical links).
#'
#' For Bernoulli, `W` and `z` share one floored working variance
#' (`binomial_weight_floor`); μ clipping uses `binomial_mu_eps`.
#'
#' @keywords internal
glm_working <- function(family, y, eta, control = NULL) {
  key <- family_key(family)
  if (identical(key, "gaussian")) {
    return(list(mu = eta, weight = rep(1, length(y)), z = y,
                var_work = rep(1, length(y))))
  }
  if (identical(key, "bernoulli")) {
    mu_eps <- if (!is.null(control$binomial_mu_eps)) control$binomial_mu_eps else 1e-5
    w_floor <- if (!is.null(control$binomial_weight_floor)) {
      control$binomial_weight_floor
    } else {
      1e-4
    }
    mu <- stats::plogis(eta)
    mu <- pmin(pmax(mu, mu_eps), 1 - mu_eps)
    var_raw <- mu * (1 - mu)
    var_work <- pmax(var_raw, w_floor)
    w <- var_work
    z <- eta + (y - mu) / var_work
    return(list(mu = mu, weight = w, z = z, var_raw = var_raw, var_work = var_work))
  }
  if (identical(key, "poisson")) {
    et <- pmin(pmax(eta, -20), 20)
    mu <- pmax(exp(et), 1e-8)
    z <- et + (y - mu) / mu
    return(list(mu = mu, weight = mu, z = z, var_work = mu))
  }
  stop("Unsupported family key: ", key)
}

#' Model deviance (Gaussian uses RSS). Optional observation `weights`.
#' @keywords internal
glm_deviance <- function(family, y, mu, weights = NULL) {
  key <- family_key(family)
  w <- normalize_weights(weights, length(y))
  if (identical(key, "gaussian")) {
    return(sum(w * (y - mu)^2))
  }
  if (identical(key, "bernoulli")) {
    mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
    return(-2 * sum(w * (y * log(mu) + (1 - y) * log(1 - mu))))
  }
  if (identical(key, "poisson")) {
    mu <- pmax(mu, 1e-12)
    term <- ifelse(y > 0, y * log(y / mu), 0) - (y - mu)
    return(2 * sum(w * term))
  }
  stop("Unsupported family")
}

init_intercept <- function(family, y, offset = NULL, weights = NULL) {
  key <- family_key(family)
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  sw <- sum(w)
  if (sw <= 0) stop("sum(weights) must be positive.", call. = FALSE)
  if (identical(key, "gaussian")) return(sum(w * (y - offset)) / sw)
  if (identical(key, "bernoulli")) {
    p <- sum(w * y) / sw
    p <- pmin(pmax(p, 0.05), 0.95)
    return(qlogis(p) - sum(w * offset) / sw)
  }
  if (identical(key, "poisson")) {
    rate <- sum(w * y / pmax(exp(offset), 1e-12)) / sw
    return(log(max(rate, 1e-8)))
  }
  sum(w * (y - offset)) / sw
}

invlink_eta <- function(family, eta) {
  key <- family_key(family)
  if (identical(key, "gaussian")) return(eta)
  if (identical(key, "bernoulli")) return(plogis(eta))
  if (identical(key, "poisson")) return(exp(pmin(pmax(eta, -20), 20)))
  eta
}

#' Full linear predictor η = offset + intercept + TT contraction.
#' @keywords internal
#' @noRd
tt_eta <- function(offset, intercept, cores, basis) {
  as.numeric(offset) + as.numeric(intercept) + tt_contraction(cores, basis)
}

#' Normalize a length-n offset (NULL → zeros).
#' @keywords internal
#' @noRd
normalize_offset <- function(offset, n) {
  n <- as.integer(n)
  if (is.null(offset)) return(rep(0, n))
  offset <- as.numeric(offset)
  if (length(offset) == 1L) offset <- rep(offset, n)
  if (length(offset) != n) {
    stop("`offset` must have length 1 or length(y) = ", n, ".", call. = FALSE)
  }
  if (anyNA(offset)) stop("`offset` contains NA.", call. = FALSE)
  offset
}

#' Normalize observation weights (NULL → ones). Must be finite and ≥ 0.
#' @keywords internal
#' @noRd
normalize_weights <- function(weights, n) {
  n <- as.integer(n)
  if (is.null(weights)) return(rep(1, n))
  weights <- as.numeric(weights)
  if (length(weights) == 1L) weights <- rep(weights, n)
  if (length(weights) != n) {
    stop("`weights` must have length 1 or length(y) = ", n, ".", call. = FALSE)
  }
  if (anyNA(weights)) stop("`weights` contains NA.", call. = FALSE)
  if (any(weights < 0)) stop("`weights` must be non-negative.", call. = FALSE)
  if (sum(weights) <= 0) stop("sum(weights) must be positive.", call. = FALSE)
  weights
}

#' True GLM + TT penalty objective at (cores, intercept).
#' @keywords internal
tt_glm_penalized_objective <- function(y, cores, intercept, basis, penalties,
                                       lambda, family, offset = NULL,
                                       weights = NULL) {
  key <- family_key(family)
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  eta <- tt_eta(offset, intercept, cores, basis)
  if (identical(key, "gaussian")) {
    rss <- sum(w * (y - eta)^2)
    nll <- 0.5 * rss
  } else if (identical(key, "bernoulli")) {
    mu <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12)
    nll <- -sum(w * (y * log(mu) + (1 - y) * log(1 - mu)))
  } else if (identical(key, "poisson")) {
    mu <- pmax(exp(pmin(pmax(eta, -20), 20)), 1e-12)
    nll <- sum(w * (mu - y * log(mu)))
  } else {
    stop("Unsupported family in tt_glm_penalized_objective", call. = FALSE)
  }
  pen <- .tt_penalty_value_grad(cores, penalties, lambda)$value
  list(value = nll + pen, nll = nll, penalty = pen, eta = eta)
}

#' Linear blend of TT cores + intercept (parameter space).
#'
#' No gauge alignment is applied: current R ALS does not re-orthogonalize
#' cores between outer iterates, so old/candidate share a continuous ALS path.
#' If future ALS adds canonicalization, align before blending.
#'
#' @keywords internal
tt_blend_params <- function(cores_old, intercept_old, cores_new, intercept_new, alpha) {
  cores <- lapply(seq_along(cores_old), function(k) {
    cores_old[[k]] + alpha * (cores_new[[k]] - cores_old[[k]])
  })
  list(
    cores = cores,
    intercept = intercept_old + alpha * (intercept_new - intercept_old)
  )
}
