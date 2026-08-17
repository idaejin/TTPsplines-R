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

#' Parametric contribution linear %*% beta (0 if absent).
#' @keywords internal
#' @noRd
tt_linear_contrib <- function(linear, beta) {
  if (is.null(linear) || is.null(beta) || length(beta) == 0L ||
      ncol(as.matrix(linear)) == 0L) {
    return(0)
  }
  as.numeric(as.matrix(linear) %*% as.numeric(beta))
}

#' Full linear predictor η = offset + intercept + linear β + smooth + TT.
#' @keywords internal
#' @noRd
tt_eta <- function(offset, intercept, cores, basis,
                   linear = NULL, beta = NULL, smooth = NULL) {
  as.numeric(offset) + as.numeric(intercept) +
    tt_linear_contrib(linear, beta) +
    tt_smooth_contrib(smooth) +
    tt_contraction(cores, basis)
}

#' Normalize unpenalized parametric design (NULL → NULL).
#'
#' Do **not** include an intercept column: [ttps()] already estimates one.
#'
#' @keywords internal
#' @noRd
normalize_linear <- function(linear, n) {
  n <- as.integer(n)
  if (is.null(linear)) return(NULL)
  # Idempotent: ttps() + ALS/PIRLS both call this; warn at most once
  if (isTRUE(attr(linear, "ttps_linear_ok")) && is.matrix(linear)) {
    if (nrow(linear) != n) {
      stop("`linear` must have ", n, " rows (same as length(y)).", call. = FALSE)
    }
    return(linear)
  }
  if (is.vector(linear) && !is.list(linear)) {
    linear <- matrix(as.numeric(linear), ncol = 1L)
  }
  linear <- as.matrix(linear)
  storage.mode(linear) <- "double"
  if (nrow(linear) != n) {
    stop("`linear` must have ", n, " rows (same as length(y)).", call. = FALSE)
  }
  if (anyNA(linear)) stop("`linear` contains NA.", call. = FALSE)
  if (ncol(linear) == 0L) return(NULL)
  # Detect an intercept column (by name or literal ones)
  ones <- which(vapply(seq_len(ncol(linear)), function(j) {
    nm <- colnames(linear)[j]
    if (!is.null(nm) && identical(nm, "(Intercept)")) return(TRUE)
    max(abs(linear[, j] - 1)) < 1e-10
  }, logical(1)))
  if (length(ones)) {
    warning(
      "`linear` includes an intercept column (",
      paste(colnames(linear)[ones], collapse = ", "),
      "); omit it - ttps() already has an intercept. ",
      "Prefer model.matrix(~ 0 + ..., data).",
      call. = FALSE
    )
  }
  if (is.null(colnames(linear))) {
    colnames(linear) <- paste0("linear", seq_len(ncol(linear)))
  }
  attr(linear, "ttps_linear_ok") <- TRUE
  linear
}

#' Weighted OLS update of (intercept, beta) given TT surface f.
#'
#' Solves min_α,β ||√w (target − offset − f − α − linear β)||².
#' When `linear` is NULL, reduces to a weighted mean for the intercept.
#'
#' @keywords internal
#' @noRd
tt_update_intercept_beta <- function(target, offset, f, linear = NULL,
                                     weights = NULL) {
  target <- as.numeric(target)
  n <- length(target)
  offset <- normalize_offset(offset, n)
  f <- if (length(f) == 1L) rep(as.numeric(f), n) else as.numeric(f)
  if (length(f) != n) stop("f must have length n.", call. = FALSE)
  w <- normalize_weights(weights, n)
  r <- target - offset - f
  if (is.null(linear) || ncol(linear) == 0L) {
    return(list(
      intercept = sum(w * r) / max(sum(w), 1e-12),
      beta = numeric(0)
    ))
  }
  L <- as.matrix(linear)
  Xb <- cbind(1, L)
  sw <- sqrt(w)
  Xw <- Xb * sw
  rw <- r * sw
  fit <- tryCatch(stats::lm.fit(Xw, rw), error = function(e) NULL)
  if (is.null(fit)) {
    xtx <- crossprod(Xw)
    diag(xtx) <- diag(xtx) + 1e-8
    coef <- as.numeric(solve(xtx, crossprod(Xw, rw)))
  } else {
    coef <- as.numeric(fit$coefficients)
    coef[is.na(coef)] <- 0
  }
  beta <- coef[-1L]
  names(beta) <- colnames(L)
  list(intercept = coef[1L], beta = beta)
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

#' True GLM + TT penalty objective at cores, intercept, and optional beta.
#' @keywords internal
tt_glm_penalized_objective <- function(y, cores, intercept, basis, penalties,
                                       lambda, family, offset = NULL,
                                       weights = NULL,
                                       penalty_mode = "global",
                                       penalty_order = 2L,
                                       linear = NULL, beta = NULL,
                                       smooth = NULL) {
  key <- family_key(family)
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  eta <- tt_eta(offset, intercept, cores, basis,
                linear = linear, beta = beta, smooth = smooth)
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
  # Always classical J_λ(Θ); penalty_mode kept for call compatibility only.
  normalize_penalty_mode(penalty_mode)
  pen <- tt_global_penalty_value(
    cores, lambda, penalty_order = penalty_order,
    cyclic = attr(basis, "cyclic")
  ) + tt_smooth_penalty_value(smooth)
  list(value = nll + pen, nll = nll, penalty = pen, eta = eta)
}

#' Linear blend of TT cores + intercept (+ optional beta).
#'
#' No gauge alignment is applied: current R ALS does not re-orthogonalize
#' cores between outer iterates, so old/candidate share a continuous ALS path.
#' If future ALS adds canonicalization, align before blending.
#'
#' @keywords internal
tt_blend_params <- function(cores_old, intercept_old, cores_new, intercept_new,
                            alpha, beta_old = NULL, beta_new = NULL,
                            smooth_old = NULL, smooth_new = NULL) {
  cores <- lapply(seq_along(cores_old), function(k) {
    cores_old[[k]] + alpha * (cores_new[[k]] - cores_old[[k]])
  })
  out <- list(
    cores = cores,
    intercept = intercept_old + alpha * (intercept_new - intercept_old)
  )
  if (!is.null(beta_old) || !is.null(beta_new)) {
    b0 <- if (is.null(beta_old)) numeric(0) else as.numeric(beta_old)
    b1 <- if (is.null(beta_new)) numeric(0) else as.numeric(beta_new)
    if (length(b0) == 0L && length(b1) == 0L) {
      out$beta <- numeric(0)
    } else if (length(b0) == length(b1)) {
      out$beta <- b0 + alpha * (b1 - b0)
      if (!is.null(names(b1))) names(out$beta) <- names(b1)
      else if (!is.null(names(b0))) names(out$beta) <- names(b0)
    } else {
      out$beta <- b1
    }
  }
  if (!is.null(smooth_old) || !is.null(smooth_new)) {
    s0 <- smooth_old
    s1 <- smooth_new
    if (is.null(s0)) {
      out$smooth <- s1
    } else if (is.null(s1) || length(s0) != length(s1)) {
      out$smooth <- s1
    } else {
      out$smooth <- lapply(seq_along(s0), function(j) {
        sm <- s1[[j]]
        sm$gamma <- as.numeric(s0[[j]]$gamma) +
          alpha * (as.numeric(s1[[j]]$gamma) - as.numeric(s0[[j]]$gamma))
        sm
      })
      names(out$smooth) <- names(s1)
      attr(out$smooth, "ttps_smooth_ok") <- TRUE
    }
  }
  out
}
