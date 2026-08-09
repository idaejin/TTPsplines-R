# Gradient descent on the same penalized likelihood as L-BFGS (not PIRLS).

#' Gradient descent over packed TT cores (analytical grad shared with L-BFGS).
#'
#' Direct minimization of the penalized Gaussian / GLM objective. Optional
#' Armijo backtracking line search (`control$gd_linesearch`).
#'
#' @keywords internal
#' @noRd
.tt_gd_optimize_cores <- function(y, basis, ranks, lambda, control,
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
  theta <- .tt_pack_cores(cores)

  is_gauss <- is.null(family) || identical(family_key(family), "gaussian")
  if (is_gauss) {
    intercept <- if (is.null(intercept0)) mean(y) else intercept0
    eval_obj <- function(th) {
      .tt_gaussian_objective(th, y, intercept, basis, template, penalties, lambda)
    }
  } else {
    fam <- normalize_family(family)
    intercept <- if (is.null(intercept0)) init_intercept(fam, y) else intercept0
    eval_obj <- function(th) {
      .tt_glm_objective(th, y, intercept, basis, template, penalties, lambda, fam)
    }
  }

  lr0 <- as.numeric(control$gd_lr %||% 1e-2)
  maxit <- as.integer(control$gd_maxit %||% 5000L)
  gtol <- as.numeric(control$gd_tol %||% 1e-7)
  do_ls <- isTRUE(control$gd_linesearch %||% TRUE)
  rho <- as.numeric(control$gd_step_factor %||% control$step_factor %||% 0.5)
  step_min <- as.numeric(control$gd_step_min %||% control$step_min %||% 1e-12)
  c_armijo <- as.numeric(control$gd_armijo_c %||% 1e-4)
  ftol <- as.numeric(control$tol %||% 1e-8)

  cur <- eval_obj(theta)
  f <- cur$value
  g <- cur$grad
  n_eval <- 1L
  n_grad <- 1L
  n_ls <- 0L
  hist <- list()
  reason <- "maxit"
  converged <- FALSE

  for (it in seq_len(maxit)) {
    gnorm <- max(abs(g))
    if (!is.finite(gnorm) || !is.finite(f)) {
      reason <- "non-finite objective/gradient"
      converged <- FALSE
      break
    }
    if (gnorm < gtol) {
      reason <- "gradient tolerance"
      converged <- TRUE
      break
    }

    # descent direction
    dir <- -g
    dir_deriv <- sum(g * dir) # = -||g||^2
    if (!(dir_deriv < 0)) {
      reason <- "not a descent direction"
      converged <- FALSE
      break
    }

    accepted <- FALSE
    gamma <- if (do_ls) lr0 else lr0
    f_new <- NA_real_
    g_new <- g
    th_new <- theta

    if (do_ls) {
      while (gamma + 1e-18 >= step_min) {
        n_ls <- n_ls + 1L
        th_try <- theta + gamma * dir
        try_obj <- eval_obj(th_try)
        n_eval <- n_eval + 1L
        # Armijo: f(x + γd) ≤ f + c γ <∇f, d>
        if (is.finite(try_obj$value) &&
            try_obj$value <= f + c_armijo * gamma * dir_deriv) {
          th_new <- th_try
          f_new <- try_obj$value
          g_new <- try_obj$grad
          n_grad <- n_grad + 1L
          accepted <- TRUE
          break
        }
        gamma <- gamma * rho
      }
    } else {
      th_try <- theta + gamma * dir
      try_obj <- eval_obj(th_try)
      n_eval <- n_eval + 1L
      if (is.finite(try_obj$value)) {
        th_new <- th_try
        f_new <- try_obj$value
        g_new <- try_obj$grad
        n_grad <- n_grad + 1L
        accepted <- TRUE
      }
    }

    if (!accepted) {
      reason <- "line search failed"
      converged <- FALSE
      break
    }

    rel <- abs(f - f_new) / max(1, abs(f))
    theta <- th_new
    f <- f_new
    g <- g_new

    if (isTRUE(control$trace)) {
      cat(sprintf(
        "  GD %4d | L=%.6g | ||g||_inf=%.3g | step=%.3g\n",
        it, f, max(abs(g)), gamma
      ))
    }
    if (isTRUE(control$trace) || it %% 50L == 0L) {
      hist[[length(hist) + 1L]] <- list(
        iter = it, value = f, grad_inf = max(abs(g)), step = gamma
      )
    }

    if (it > 1L && rel < ftol && max(abs(g)) < max(gtol * 10, 1e-5)) {
      reason <- "relative objective change"
      converged <- TRUE
      break
    }
    if (it == maxit) reason <- "maxit"
  }

  cores <- .tt_unpack_cores(theta, template)
  list(
    cores = cores,
    intercept = intercept,
    penalties = penalties,
    value = f,
    grad_inf = max(abs(g)),
    n_iter = it,
    n_eval = n_eval,
    n_grad = n_grad,
    n_ls = n_ls,
    converged = converged,
    reason = reason,
    history = hist,
    is_gauss = is_gauss,
    family = if (is_gauss) NULL else normalize_family(family)
  )
}

#' Gradient-descent fit (Gaussian or GLM) with fixed λ or outer cGCV.
#'
#' Same penalized objective / analytical gradient as [tt_lbfgs_fit()].
#'
#' @keywords internal
#' @noRd
tt_gd_fit <- function(y, basis, ranks, lambda_spec, control,
                      penalty_order = 2, init_cores = NULL,
                      family = NULL, intercept0 = NULL) {
  method <- lambda_spec$method
  lambda <- lambda_spec$values %||% lambda_spec$lambda0
  t0 <- proc.time()[["elapsed"]]
  n_eval <- 0L
  n_outer <- 0L
  history <- list()

  is_gauss <- is.null(family) || identical(family_key(family), "gaussian")

  run_once <- function(lam, cores0, intercept0 = NULL) {
    .tt_gd_optimize_cores(
      y, basis, ranks, lam, control, penalty_order, cores0,
      family = family, intercept0 = intercept0
    )
  }

  if (identical(method, "fixed")) {
    fit0 <- run_once(lambda, init_cores, intercept0 = intercept0)
    cores <- fit0$cores
    intercept <- fit0$intercept
    penalties <- fit0$penalties
    n_iter <- fit0$n_iter
    converged <- isTRUE(fit0$converged)
    reason <- fit0$reason
    n_outer <- 1L
    n_eval <- fit0$n_eval
    history <- fit0$history
  } else {
    cores <- init_cores
    intercept <- if (!is.null(intercept0)) {
      intercept0
    } else if (is_gauss) {
      mean(y)
    } else {
      init_intercept(normalize_family(family), y)
    }
    penalties <- NULL
    prev_lam <- lambda
    converged <- FALSE
    n_iter <- 0L
    reason <- "maxit"
    for (outer in seq_len(control$outer_maxit)) {
      fit0 <- run_once(lambda, cores, intercept0 = intercept)
      cores <- fit0$cores
      intercept <- fit0$intercept
      penalties <- fit0$penalties
      n_iter <- n_iter + fit0$n_iter
      n_eval <- n_eval + fit0$n_eval
      reason <- fit0$reason

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
      history[[outer]] <- list(outer = outer, lambda = lambda, value = fit0$value)
      dlog <- max(abs(log(lambda) - log(pmax(prev_lam, 1e-12))))
      if (dlog < control$outer_tol && outer > 1L) {
        converged <- TRUE
        reason <- "outer lambda convergence"
        break
      }
      prev_lam <- lambda
    }
    if (n_outer >= control$outer_maxit) converged <- isTRUE(fit0$converged)
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
    converged = isTRUE(converged) && is.finite(deviance),
    convergence = list(
      overall = isTRUE(converged) && is.finite(deviance),
      pirls = NA,
      als = NA,
      reason = reason,
      n_pirls = NA_integer_,
      n_als_sweeps = NA_integer_,
      n_step_halvings = 0L,
      n_opt_iter = as.integer(n_iter)
    ),
    method_lambda = method,
    optimizer = "GD",
    backend = "R"
  )
}
