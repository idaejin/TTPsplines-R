# Global L-BFGS over packed TT cores (same penalized objective as ALS).

.tt_pack_cores <- function(cores) {
  as.numeric(unlist(lapply(cores, as.vector), use.names = FALSE))
}

.tt_unpack_cores <- function(theta, template) {
  out <- template
  pos <- 1L
  for (k in seq_along(template)) {
    n <- length(template[[k]])
    out[[k]][] <- theta[pos:(pos + n - 1L)]
    pos <- pos + n
  }
  out
}

# .tt_penalty_value_grad lives in penalties.R (shared with PIRLS objective)

.tt_sse_grad_cores <- function(resid, cores, basis) {
  # resid = y - intercept - f  (or working residual for GLM score)
  d <- length(cores)
  L <- left_interfaces(cores, basis)
  R <- right_interfaces(cores, basis)
  grads <- vector("list", d)
  for (k in seq_len(d)) {
    Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
    grads[[k]] <- array(-as.numeric(crossprod(Xk, resid)), dim(cores[[k]]))
  }
  grads
}

.tt_gaussian_objective <- function(theta, y, intercept, basis, template,
                                   penalties, lambda, offset = NULL,
                                   weights = NULL,
                                   penalty_mode = "global",
                                   penalty_order = 2L) {
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  cores <- .tt_unpack_cores(theta, template)
  f <- tt_contraction(cores, basis)
  resid <- y - offset - intercept - f
  sse <- 0.5 * sum(w * resid^2)
  normalize_penalty_mode(penalty_mode)
  pen <- .tt_penalty_value_grad_full(
    cores, lambda, penalty_order = penalty_order,
    cyclic = attr(basis, "cyclic")
  )
  list(
    value = sse + pen$value,
    grad = .tt_pack_cores(
      mapply(function(gs, gp) gs + gp,
             .tt_sse_grad_cores(w * resid, cores, basis),
             pen$grads, SIMPLIFY = FALSE)
    ),
    cores = cores,
    eta = offset + intercept + f,
    resid = resid
  )
}

.tt_glm_objective <- function(theta, y, intercept, basis, template,
                              penalties, lambda, fam, offset = NULL,
                              weights = NULL,
                              penalty_mode = "global",
                              penalty_order = 2L) {
  offset <- normalize_offset(offset, length(y))
  w <- normalize_weights(weights, length(y))
  cores <- .tt_unpack_cores(theta, template)
  eta <- offset + intercept + tt_contraction(cores, basis)
  mu <- invlink_eta(fam, eta)
  key <- family_key(fam)
  nll <- switch(
    key,
    poisson = {
      mu <- pmax(mu, 1e-12)
      sum(w * (mu - y * log(mu)))
    },
    bernoulli = {
      mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
      -sum(w * (y * log(mu) + (1 - y) * log(1 - mu)))
    },
    stop("LBFGS GLM only poisson/binomial in v0.", call. = FALSE)
  )
  normalize_penalty_mode(penalty_mode)
  pen <- .tt_penalty_value_grad_full(
    cores, lambda, penalty_order = penalty_order,
    cyclic = attr(basis, "cyclic")
  )
  grads_nll <- .tt_sse_grad_cores(w * (y - mu), cores, basis)
  list(
    value = nll + pen$value,
    grad = .tt_pack_cores(
      mapply(function(gs, gp) gs + gp, grads_nll, pen$grads, SIMPLIFY = FALSE)
    ),
    cores = cores,
    eta = eta,
    mu = mu
  )
}

.tt_lbfgs_optimize_cores <- function(y, basis, ranks, lambda, control,
                                     penalty_order, init_cores,
                                     family = NULL, intercept0 = NULL,
                                     offset = NULL, weights = NULL) {
  offset <- normalize_offset(offset, length(y))
  weights <- normalize_weights(weights, length(y))
  d <- length(basis)
  p <- ncol(basis[[1]])
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- tt_core_penalties_from_basis(ranks, basis, penalty_order)
  template <- cores
  theta0 <- .tt_pack_cores(cores)
  penalty_mode <- normalize_penalty_mode(control$penalty_mode %||% "global")
  cyclic <- attr(basis, "cyclic")

  is_gauss <- is.null(family) || identical(family_key(family), "gaussian")
  if (is_gauss) {
    intercept <- if (is.null(intercept0)) {
      sum(weights * (y - offset)) / sum(weights)
    } else {
      intercept0
    }
    fn <- function(th) {
      .tt_gaussian_objective(
        th, y, intercept, basis, template, penalties, lambda,
        offset = offset, weights = weights,
        penalty_mode = penalty_mode, penalty_order = penalty_order
      )$value
    }
    gr <- function(th) {
      .tt_gaussian_objective(
        th, y, intercept, basis, template, penalties, lambda,
        offset = offset, weights = weights,
        penalty_mode = penalty_mode, penalty_order = penalty_order
      )$grad
    }
  } else {
    fam <- normalize_family(family)
    intercept <- if (is.null(intercept0)) {
      init_intercept(fam, y, offset = offset, weights = weights)
    } else {
      intercept0
    }
    fn <- function(th) {
      .tt_glm_objective(
        th, y, intercept, basis, template, penalties, lambda, fam,
        offset = offset, weights = weights,
        penalty_mode = penalty_mode, penalty_order = penalty_order
      )$value
    }
    gr <- function(th) {
      .tt_glm_objective(
        th, y, intercept, basis, template, penalties, lambda, fam,
        offset = offset, weights = weights,
        penalty_mode = penalty_mode, penalty_order = penalty_order
      )$grad
    }
  }

  opt <- stats::optim(
    theta0, fn, gr,
    method = "L-BFGS-B",
    control = list(
      maxit = as.integer(control$lbfgs_maxit),
      factr = 1e7,
      pgtol = 1e-8,
      trace = if (isTRUE(control$trace)) 1L else 0L
    )
  )
  cores <- .tt_unpack_cores(opt$par, template)
  list(
    cores = cores,
    intercept = intercept,
    penalties = penalties,
    penalty_mode = penalty_mode,
    opt = opt,
    is_gauss = is_gauss,
    family = if (is_gauss) NULL else normalize_family(family)
  )
}

#' One conditional cGCV pass over all cores (shared with outer LBFGS/Adam).
#' @keywords internal
tt_cgcv_update_lambdas <- function(y, cores, intercept, basis, penalties, lambda,
                                   control, weight = NULL, z = NULL, offset = NULL,
                                   penalty_order = 2L) {
  d <- length(cores)
  ranks <- integer(d + 1L)
  ranks[1] <- dim(cores[[1]])[1]
  for (k in seq_len(d)) ranks[k + 1L] <- dim(cores[[k]])[3]
  p <- ncol(basis[[1]])
  offset <- normalize_offset(offset, length(y))
  target <- if (is.null(z)) y - offset - intercept else z - offset - intercept
  n_eval <- 0L
  use_spec <- isTRUE(control$use_spectral_gcv)
  use_full <- TRUE
  cyclic <- attr(basis, "cyclic")
  for (k in seq_len(d)) {
    L <- left_interfaces(cores, basis)
    R <- right_interfaces(cores, basis)
    Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
    pen_k <- tt_conditional_penalty_full(
      cores, k, lambda, penalty_order = penalty_order, cyclic = cyclic
    )
    Pk <- pen_k$P_own
    P0 <- pen_k$P_other
    penalties[[k]] <- Pk
    ws <- make_core_workspace(
      target, Xk, Pk, lambda[k],
      control$lambda_bounds, control$tol_lambda,
      weight = weight, use_spectral = use_spec, P0 = P0
    )
    upd <- update_lambda("cGCV", ws)
    cores[[k]] <- array(upd$g, c(ranks[k], p, ranks[k + 1L]))
    lambda[k] <- upd$lambda
    n_eval <- n_eval + upd$n_eval
  }
  list(cores = cores, lambda = lambda, n_eval = n_eval, penalties = penalties)
}

#' L-BFGS fit (Gaussian or GLM) with fixed λ or outer cGCV.
#' @keywords internal
tt_lbfgs_fit <- function(y, basis, ranks, lambda_spec, control,
                         penalty_order = 2, init_cores = NULL,
                         family = NULL, intercept0 = NULL, offset = NULL,
                         weights = NULL) {
  method <- lambda_spec$method
  lambda <- lambda_spec$values %||% lambda_spec$lambda0
  offset <- normalize_offset(offset, length(y))
  weights <- normalize_weights(weights, length(y))
  t0 <- proc.time()[["elapsed"]]
  n_eval <- 0L
  n_outer <- 0L
  history <- list()
  opt_convergence <- NA_integer_

  is_gauss <- is.null(family) || identical(family_key(family), "gaussian")

  run_once <- function(lam, cores0, intercept0 = NULL) {
    .tt_lbfgs_optimize_cores(
      y, basis, ranks, lam, control, penalty_order, cores0,
      family = family, intercept0 = intercept0, offset = offset,
      weights = weights
    )
  }

  if (identical(method, "fixed")) {
    fit0 <- run_once(lambda, init_cores, intercept0 = intercept0)
    cores <- fit0$cores
    intercept <- fit0$intercept
    penalties <- fit0$penalties
    n_iter <- fit0$opt$counts[["function"]]
    opt_convergence <- fit0$opt$convergence
    # 0 = converged; 1 = maxit — accept if finite objective
    converged <- opt_convergence %in% c(0L, 1L)
    n_outer <- 1L
  } else {
    # Outer alternation: LBFGS cores <-> conditional cGCV λ
    cores <- init_cores
    intercept <- if (!is.null(intercept0)) {
      intercept0
    } else if (is_gauss) {
      sum(weights * (y - offset)) / sum(weights)
    } else {
      init_intercept(normalize_family(family), y, offset = offset, weights = weights)
    }
    penalties <- NULL
    prev_lam <- lambda
    converged <- FALSE
    n_iter <- 0L
    for (outer in seq_len(control$outer_maxit)) {
      fit0 <- run_once(lambda, cores, intercept0 = intercept)
      cores <- fit0$cores
      intercept <- fit0$intercept
      penalties <- fit0$penalties
      n_iter <- n_iter + fit0$opt$counts[["function"]]
      opt_convergence <- fit0$opt$convergence

      if (is_gauss) {
        upd <- tt_cgcv_update_lambdas(
          y, cores, intercept, basis, penalties, lambda, control,
          weight = weights, offset = offset, penalty_order = penalty_order
        )
      } else {
        fam <- normalize_family(family)
        eta_cur <- offset + intercept + tt_contraction(cores, basis)
        work <- glm_working(fam, y, eta_cur)
        upd <- tt_cgcv_update_lambdas(
          y, cores, intercept, basis, penalties, lambda, control,
          weight = work$weight * weights, z = work$z, offset = offset,
          penalty_order = penalty_order
        )
      }
      cores <- upd$cores
      lambda <- upd$lambda
      if (!is.null(upd$penalties)) penalties <- upd$penalties
      n_eval <- n_eval + upd$n_eval
      n_outer <- outer
      history[[outer]] <- list(outer = outer, lambda = lambda, value = fit0$opt$value)
      if (control$trace) {
        cat(sprintf("  LBFGS-cGCV outer %2d | λ=%s\n",
                    outer, paste(sprintf("%.3g", lambda), collapse = ",")))
      }
      dlog <- max(abs(log(lambda) - log(pmax(prev_lam, 1e-12))))
      if (dlog < control$outer_tol && outer > 1L) {
        converged <- TRUE
        break
      }
      prev_lam <- lambda
    }
    if (n_outer >= control$outer_maxit) converged <- TRUE
  }

  if (is_gauss) {
    eta <- offset + intercept + tt_contraction(cores, basis)
    mu <- eta
    deviance <- sum(weights * (y - eta)^2)
  } else {
    fam <- normalize_family(family)
    eta <- offset + intercept + tt_contraction(cores, basis)
    mu <- invlink_eta(fam, eta)
    deviance <- glm_deviance(fam, y, mu, weights = weights)
  }

  list(
    cores = cores,
    intercept = intercept,
    lambda = lambda,
    ranks = ranks,
    eta = eta,
    mu = mu,
    deviance = deviance,
    n_sweeps = NA_integer_,
    n_pirls = NA_integer_,
    n_criterion_evals = n_eval,
    n_opt_iter = as.integer(n_iter),
    n_outer = as.integer(n_outer),
    history = history,
    penalties = penalties,
    elapsed = proc.time()[["elapsed"]] - t0,
    converged = isTRUE(converged) && is.finite(deviance),
    convergence = list(
      overall = isTRUE(converged) && is.finite(deviance),
      pirls = NA,
      als = NA,
      reason = if (isTRUE(converged)) {
        if (isTRUE(opt_convergence == 1L)) {
          "L-BFGS maxit (finite objective)"
        } else {
          "L-BFGS"
        }
      } else {
        "L-BFGS did not report convergence"
      },
      n_pirls = NA_integer_,
      n_als_sweeps = NA_integer_,
      n_step_halvings = 0L,
      n_opt_iter = as.integer(n_iter)
    ),
    method_lambda = method,
    optimizer = "LBFGS",
    backend = "R"
  )
}
