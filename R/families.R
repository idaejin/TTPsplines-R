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
#' @keywords internal
glm_working <- function(family, y, eta) {
  key <- family_key(family)
  if (identical(key, "gaussian")) {
    return(list(mu = eta, weight = rep(1, length(y)), z = y))
  }
  if (identical(key, "bernoulli")) {
    mu <- plogis(eta)
    mu <- pmin(pmax(mu, 1e-5), 1 - 1e-5)
    var <- mu * (1 - mu)
    w <- pmax(var, 1e-4)
    z <- eta + (y - mu) / var
    return(list(mu = mu, weight = w, z = z))
  }
  if (identical(key, "poisson")) {
    et <- pmin(pmax(eta, -20), 20)
    mu <- pmax(exp(et), 1e-8)
    z <- et + (y - mu) / mu
    return(list(mu = mu, weight = mu, z = z))
  }
  stop("Unsupported family key: ", key)
}

#' Model deviance (Gaussian uses RSS).
#' @keywords internal
glm_deviance <- function(family, y, mu) {
  key <- family_key(family)
  if (identical(key, "gaussian")) {
    return(sum((y - mu)^2))
  }
  if (identical(key, "bernoulli")) {
    mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
    return(-2 * sum(y * log(mu) + (1 - y) * log(1 - mu)))
  }
  if (identical(key, "poisson")) {
    mu <- pmax(mu, 1e-12)
    term <- ifelse(y > 0, y * log(y / mu), 0) - (y - mu)
    return(2 * sum(term))
  }
  stop("Unsupported family")
}

init_intercept <- function(family, y) {
  key <- family_key(family)
  if (identical(key, "gaussian")) return(mean(y))
  if (identical(key, "bernoulli")) {
    p <- pmin(pmax(mean(y), 0.05), 0.95)
    return(qlogis(p))
  }
  if (identical(key, "poisson")) return(log(max(mean(y), 0.1)))
  mean(y)
}

invlink_eta <- function(family, eta) {
  key <- family_key(family)
  if (identical(key, "gaussian")) return(eta)
  if (identical(key, "bernoulli")) return(plogis(eta))
  if (identical(key, "poisson")) return(exp(pmin(pmax(eta, -20), 20)))
  eta
}
