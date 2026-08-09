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

.tt_penalty_value_grad <- function(cores, penalties, lambda) {
  val <- 0
  grads <- vector("list", length(cores))
  for (k in seq_along(cores)) {
    g <- as.numeric(cores[[k]])
    Pk <- penalties[[k]]
    Pg <- as.numeric(Pk %*% g)
    val <- val + 0.5 * lambda[k] * sum(g * Pg)
    grads[[k]] <- array(lambda[k] * Pg, dim(cores[[k]]))
  }
  list(value = val, grads = grads)
}

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
                                   penalties, lambda) {
  cores <- .tt_unpack_cores(theta, template)
  f <- tt_contraction(cores, basis)
  resid <- y - intercept - f
  sse <- 0.5 * sum(resid^2)
  pen <- .tt_penalty_value_grad(cores, penalties, lambda)
  list(
    value = sse + pen$value,
    grad = .tt_pack_cores(
      mapply(function(gs, gp) gs + gp, .tt_sse_grad_cores(resid, cores, basis),
             pen$grads, SIMPLIFY = FALSE)
    ),
    cores = cores,
    eta = intercept + f,
    resid = resid
  )
}

.tt_glm_objective <- function(theta, y, intercept, basis, template,
                              penalties, lambda, fam) {
  cores <- .tt_unpack_cores(theta, template)
  eta <- intercept + tt_contraction(cores, basis)
  mu <- invlink_eta(fam, eta)
  # canonical exponential-family score: -(y - mu) for poisson/bernoulli
  key <- family_key(fam)
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
    stop("LBFGS GLM only poisson/binomial in v0.", call. = FALSE)
  )
  pen <- .tt_penalty_value_grad(cores, penalties, lambda)
  score_resid <- y - mu # grad_eta NLL = -(y-mu); SSE-style helper uses -X'resid
  # .tt_sse_grad_cores(resid) returns -X'resid; we need -X'(y-mu) = X'(mu-y)
  # so pass resid = y - mu to get -X'(y-mu) = desired NLL grad.
  grads_nll <- .tt_sse_grad_cores(y - mu, cores, basis)
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
                                     family = NULL, intercept0 = NULL) {
  d <- length(basis)
  p <- ncol(basis[[1]])
  if (is.null(init_cores)) {
    cores <- initialize_tt_cores(p, ranks, seed = control$seed, sd = control$init_sd)
  } else {
    cores <- init_cores
  }
  penalties <- lapply(seq_len(d), function(k) {
    core_penalty(ranks[k], p, ranks[k + 1L], penalty_order)
  })
  template <- cores
  theta0 <- .tt_pack_cores(cores)

  is_gauss <- is.null(family) || identical(family_key(family), "gaussian")
  if (is_gauss) {
    intercept <- if (is.null(intercept0)) mean(y) else intercept0
    fn <- function(th) {
      .tt_gaussian_objective(th, y, intercept, basis, template, penalties, lambda)$value
    }
    gr <- function(th) {
      .tt_gaussian_objective(th, y, intercept, basis, template, penalties, lambda)$grad
    }
  } else {
    fam <- normalize_family(family)
    intercept <- if (is.null(intercept0)) init_intercept(fam, y) else intercept0
    fn <- function(th) {
      .tt_glm_objective(th, y, intercept, basis, template, penalties, lambda, fam)$value
    }
    gr <- function(th) {
      .tt_glm_objective(th, y, intercept, basis, template, penalties, lambda, fam)$grad
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
    opt = opt,
    is_gauss = is_gauss,
    family = if (is_gauss) NULL else normalize_family(family)
  )
}

#' One conditional cGCV pass over all cores (shared with outer LBFGS/Adam).
#' @keywords internal
tt_cgcv_update_lambdas <- function(y, cores, intercept, basis, penalties, lambda,
                                   control, weight = NULL, z = NULL) {
  d <- length(cores)
  ranks <- integer(d + 1L)
  ranks[1] <- dim(cores[[1]])[1]
  for (k in seq_len(d)) ranks[k + 1L] <- dim(cores[[k]])[3]
  p <- ncol(basis[[1]])
  target <- if (is.null(z)) y - intercept else z - intercept
  n_eval <- 0L
  use_spec <- isTRUE(control$use_spectral_gcv)
  for (k in seq_len(d)) {
    L <- left_interfaces(cores, basis)
    R <- right_interfaces(cores, basis)
    Xk <- tt_design_core(L[[k]], R[[k]], basis[[k]])
    ws <- make_core_workspace(
      target, Xk, penalties[[k]], lambda[k],
      control$lambda_bounds, control$tol_lambda,
      weight = weight, use_spectral = use_spec
    )
    upd <- update_lambda("cGCV", ws)
    cores[[k]] <- array(upd$g, c(ranks[k], p, ranks[k + 1L]))
    lambda[k] <- upd$lambda
    n_eval <- n_eval + upd$n_eval
  }
  list(cores = cores, lambda = lambda, n_eval = n_eval)
}

#' L-BFGS fit (Gaussian or GLM) with fixed λ or outer cGCV.
#' @keywords internal
tt_lbfgs_fit <- function(y, basis, ranks, lambda_spec, control,
                         penalty_order = 2, init_cores = NULL,
                         family = NULL) {
  method <- lambda_spec$method
  lambda <- lambda_spec$values %||% lambda_spec$lambda0
  t0 <- proc.time()[["elapsed"]]
  n_eval <- 0L
  n_outer <- 0L
  history <- list()

  is_gauss <- is.null(family) || identical(family_key(family), "gaussian")

  run_once <- function(lam, cores0, intercept0 = NULL) {
    .tt_lbfgs_optimize_cores(
      y, basis, ranks, lam, control, penalty_order, cores0,
      family = family, intercept0 = intercept0
    )
  }

  if (identical(method, "fixed")) {
    fit0 <- run_once(lambda, init_cores)
    cores <- fit0$cores
    intercept <- fit0$intercept
    penalties <- fit0$penalties
    n_iter <- fit0$opt$counts[["function"]]
    converged <- fit0$opt$convergence == 0L
    n_outer <- 1L
  } else {
    # Outer alternation: LBFGS cores <-> conditional cGCV λ
    cores <- init_cores
    intercept <- if (is_gauss) mean(y) else init_intercept(normalize_family(family), y)
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

      if (is_gauss) {
        upd <- tt_cgcv_update_lambdas(
          y, cores, intercept, basis, penalties, lambda, control
        )
      } else {
        fam <- normalize_family(family)
        eta_cur <- intercept + tt_contraction(cores, basis)
        work <- glm_working(fam, y, eta_cur)
        upd <- tt_cgcv_update_lambdas(
          y, cores, intercept, basis, penalties, lambda, control,
          weight = work$weight, z = work$z
        )
      }
      cores <- upd$cores
      lambda <- upd$lambda
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
    eta <- intercept + tt_contraction(cores, basis)
    mu <- eta
    deviance <- sum((y - eta)^2)
  } else {
    fam <- normalize_family(family)
    eta <- intercept + tt_contraction(cores, basis)
    mu <- invlink_eta(fam, eta)
    deviance <- glm_deviance(fam, y, mu)
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
    converged = isTRUE(converged),
    method_lambda = method,
    optimizer = "LBFGS",
    backend = "R"
  )
}
