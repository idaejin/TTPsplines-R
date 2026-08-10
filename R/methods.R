#' @export
#' @export
print.ttpspline <- function(x, ...) {
  cat("Tensor-Train P-spline fit\n")
  cat(sprintf("  Family: %s (%s)\n", x$family$family, x$family$link))
  cat(sprintf("  n = %d, d = %d, k = %d\n", x$n, x$d, x$k))
  cat(sprintf("  TT rank: %s | chain: %s\n",
              .tt_rank_label(x), paste(x$rank, collapse = "-")))
  opt_used <- x$optimizer_used %||% x$optimizer %||% "ALS"
  opt_req <- x$optimizer_requested %||% opt_used
  if (!identical(opt_req, opt_used)) {
    cat(sprintf("  Optimizer: %s (requested: %s) | backend: %s\n",
                opt_used, opt_req, x$backend))
  } else {
    cat(sprintf("  Optimizer: %s | backend: %s\n", opt_used, x$backend))
  }
  cat(sprintf("  Lambda (%s): %s\n",
              x$lambda_method,
              paste(sprintf("%.4g", x$lambda), collapse = ", ")))
  cat(sprintf("  TT stored: %s | full: %s | CR: %.1fx\n",
              format(x$npar_tt, big.mark = ","),
              format(x$npar_dense, big.mark = ","),
              x$compression_ratio))
  if (is.finite(x$edf)) {
    cat(sprintf("  EDF: %.2f (%.0f%% of stored parameters)\n",
                x$edf, 100 * x$edf / max(x$npar_tt, 1)))
  }
  cat(sprintf("  Deviance/RSS: %.6g | time: %.3fs\n",
              x$deviance, x$timing))
  invisible(x)
}

#' Parametric coef table for intercept + `linear` (glm-style).
#'
#' SEs treat the fitted TT surface as an offset: Fisher / WLS information
#' for `Z = [1 | linear]` at the fitted eta, conditional on TT rank and lambda.
#' Rank-deficient columns (e.g. full DOW dummies) get `NA` SE like `lm`/`glm`.
#'
#' @keywords internal
#' @noRd
.tt_parametric_coef_table <- function(object) {
  n <- as.integer(object$n %||% length(object$y))
  beta <- object$beta %||% numeric(0)
  L <- object$linear
  has_L <- !is.null(L) && length(beta) > 0L

  Z <- matrix(1, n, 1L)
  colnames(Z) <- "(Intercept)"
  est <- c(`(Intercept)` = as.numeric(object$intercept))

  if (has_L) {
    L <- as.matrix(L)
    if (nrow(L) != n) {
      stop("`linear` rows do not match n.", call. = FALSE)
    }
    if (is.null(colnames(L))) {
      colnames(L) <- names(beta) %||% paste0("linear", seq_len(ncol(L)))
    }
    Z <- cbind(Z, L)
    b <- as.numeric(beta)
    names(b) <- colnames(L)
    est <- c(est, b)
  }
  # Align names with design columns
  names(est) <- colnames(Z)

  eta <- object$linear.predictors
  if (is.null(eta) || length(eta) != n) {
    off <- object$offset %||% rep(0, n)
    basis <- .tt_fit_basis(object)
    eta <- as.numeric(
      tt_eta(off, object$intercept, object$cores, basis,
             linear = if (has_L) L else NULL,
             beta = if (has_L) beta else NULL)
    )
  }

  ww <- .tt_inference_weights(object, eta)
  w <- pmax(as.numeric(ww$weight), 0)
  key <- object$family_key %||% family_key(object$family)

  p_lin <- ncol(Z)
  n_obs <- sum(w > 0)
  edf_tt <- if (is.finite(object$edf)) as.numeric(object$edf) else 0
  df_res <- max(n_obs - edf_tt - p_lin, 1)

  if (identical(key, "gaussian")) {
    rss <- sum(w * (as.numeric(object$y) - as.numeric(eta))^2)
    dispersion <- rss / df_res
    if (!is.finite(dispersion) || dispersion <= 0) {
      dispersion <- max(rss / max(n_obs - 1, 1), .Machine$double.eps)
    }
    use_t <- TRUE
  } else {
    dispersion <- 1
    use_t <- FALSE
  }

  sw <- sqrt(w)
  Zw <- Z * sw
  qZ <- qr(Zw, tol = 1e-10)
  rank <- as.integer(qZ$rank)
  pivot <- as.integer(qZ$pivot)
  se <- rep(NA_real_, p_lin)

  if (rank > 0L) {
    R1 <- qr.R(qZ)[seq_len(rank), seq_len(rank), drop = FALSE]
    Ri <- tryCatch(backsolve(R1, diag(rank)), error = function(e) NULL)
    if (!is.null(Ri)) {
      cov_piv <- dispersion * tcrossprod(Ri)
      se_piv <- sqrt(pmax(diag(cov_piv), 0))
      se[pivot[seq_len(rank)]] <- se_piv
    }
  }

  zval <- est / se
  if (use_t) {
    pval <- 2 * stats::pt(-abs(zval), df = df_res)
    cn <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  } else {
    pval <- 2 * stats::pnorm(-abs(zval))
    cn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }

  coef <- cbind(Estimate = est, `Std. Error` = se,
                stat = zval, p = pval)
  colnames(coef) <- cn
  rownames(coef) <- names(est)

  list(
    coefficients = coef,
    dispersion = dispersion,
    df.residual = if (use_t) df_res else NA_real_,
    rank = rank,
    aliased = is.na(se),
    use_t = use_t,
    note = paste0(
      "Parametric SEs treat the fitted TT surface as an offset ",
      "(conditional on TT rank and lambda; not joint with TT cores)."
    )
  )
}

#' @export
summary.ttpspline <- function(object, ...) {
  tab <- .tt_parametric_coef_table(object)
  out <- object
  out$coefficients <- tab$coefficients
  out$dispersion <- tab$dispersion
  out$df.residual <- tab$df.residual
  out$aliased <- tab$aliased
  out$parametric_rank <- tab$rank
  out$parametric_note <- tab$note
  out$use_t <- tab$use_t
  class(out) <- c("summary.ttpspline", "ttpspline")
  out
}

#' @export
print.summary.ttpspline <- function(x, digits = max(3L, getOption("digits") - 3L),
                                    signif.stars = getOption("show.signif.stars"),
                                    ...) {
  cat("Tensor-Train P-spline fit\n\n")
  cat(sprintf("Family:                 %s\n", x$family$family))
  cat(sprintf("Link:                   %s\n", x$family$link))
  cat(sprintf("Observations:           %d\n", x$n))
  cat(sprintf("Dimensions d:           %d\n", x$d))
  cat(sprintf("Basis size k:           %d\n", x$k))
  cat(sprintf("Degree:                 %d\n", x$degree))
  cat(sprintf("Penalty order:          %d\n", x$penalty_order))
  cat(sprintf("Penalty mode:           %s\n",
              x$penalty_mode %||% x$control$penalty_mode %||% "global"))
  cat(sprintf("TT rank:                %s\n", .tt_rank_label(x)))
  cat(sprintf("Rank chain:             %s\n", paste(x$rank, collapse = "-")))
  cat(sprintf("Full tensor coeffs:     %s\n",
              format(x$npar_dense, big.mark = ",", scientific = FALSE)))
  cat(sprintf("TT stored parameters:   %s\n",
              format(x$npar_tt, big.mark = ",", scientific = FALSE)))
  if (!is.null(x$npar_tt_intrinsic) && is.finite(x$npar_tt_intrinsic)) {
    cat(sprintf("TT intrinsic dimension: %s\n",
                format(x$npar_tt_intrinsic, big.mark = ",", scientific = FALSE)))
  }
  cat(sprintf("Compression ratio:      %.2fx\n", x$compression_ratio))
  cat("\n")
  opt_used <- x$optimizer_used %||% x$optimizer %||% "ALS"
  opt_req <- x$optimizer_requested %||% opt_used
  opt_reason <- x$optimizer_reason %||% NA_character_
  cat(sprintf("Requested optimizer:    %s\n", opt_req))
  cat(sprintf("Selected optimizer:     %s\n", opt_used))
  if (!is.na(opt_reason) && nzchar(opt_reason)) {
    cat(sprintf("Reason:                 %s\n", opt_reason))
  }
  cat(sprintf("Backend:                %s\n", x$backend))
  cat(sprintf("Lambda method:          %s\n", x$lambda_method))
  cat(sprintf("Lambda:                 %s\n",
              paste(sprintf("%.6g", x$lambda), collapse = ", ")))
  if (!is.null(x$lambda_boundary) && length(x$lambda_boundary) &&
      (identical(x$lambda_method, "cGCV") || isTRUE(x$lambda_at_boundary))) {
    cat(sprintf("Lambda boundary:        %s\n",
                paste(x$lambda_boundary, collapse = ", ")))
    if (!is.null(x$lambda_bounds) && length(x$lambda_bounds) == 2L) {
      cat(sprintf("Lambda search bounds:   [%s, %s]\n",
                  format(x$lambda_bounds[1], digits = 4),
                  format(x$lambda_bounds[2], digits = 4)))
    }
    if (isTRUE(x$lambda_at_boundary)) {
      cat("Note: lambda at search bound - for init diagnostics see ttps_multistart();\n")
      cat("      stable hits across starts suggest margin/null/r, not a bad seed.\n")
    }
  }

  # ---- Parametric coefficients (glm-style table) ----
  if (!is.null(x$coefficients)) {
    cat("\nParametric coefficients:\n")
    stats::printCoefmat(
      x$coefficients,
      digits = digits,
      signif.stars = signif.stars,
      na.print = "NA",
      ...
    )
    if (isTRUE(x$use_t) && is.finite(x$dispersion) && is.finite(x$df.residual)) {
      cat(sprintf(
        "\n(Dispersion parameter phi = %s on %s degrees of freedom)\n",
        format(x$dispersion, digits = digits),
        format(trunc(x$df.residual))
      ))
    } else if (!isTRUE(x$use_t)) {
      cat("\n(Dispersion parameter for Poisson/Bernoulli family taken to be 1)\n")
    }
    if (!is.null(x$parametric_note)) {
      cat(x$parametric_note, "\n", sep = "")
    }
    if (any(x$aliased %||% FALSE)) {
      cat("Aliased coefficients have Std. Error = NA (rank-deficient design).\n")
    }
  }

  cat("\n")
  cat(sprintf("Deviance / RSS:         %.6g\n", x$deviance))
  if (is.finite(x$edf)) {
    cat(sprintf("EDF (joint linearized): %.2f\n", x$edf))
    cat(sprintf("EDF / stored parameters: %.2f\n",
                x$edf / max(x$npar_tt, 1)))
    if (!is.null(x$npar_tt_intrinsic) && is.finite(x$npar_tt_intrinsic) &&
        x$npar_tt_intrinsic > 0) {
      cat(sprintf("EDF / intrinsic TT dim.: %.2f\n",
                  x$edf / x$npar_tt_intrinsic))
    }
  } else {
    note <- x$edf_note %||% "not computed"
    cat(sprintf("EDF:                    NA (%s)\n", note))
  }
  cat("\nInference (TT surface):\n")
  inf <- x$inference
  if (is.null(inf) && is.environment(x$._inf) && !is.null(x$._inf$data)) {
    inf <- x$._inf$data
  }
  if (!is.null(inf) && isTRUE(inf$ready)) {
    cat(sprintf("  Covariance:              Bayesian penalized-spline\n"))
    cat(sprintf("  Conditional on rank:     yes\n"))
    cat(sprintf("  Conditional on lambda:   yes\n"))
    cat(sprintf("  Smoothing uncertainty:   not included\n"))
    cat(sprintf("  Rank uncertainty:        not included\n"))
    cat(sprintf("  Scale estimator:         %s\n",
                inf$scale_estimator %||% "RSS/(n-edf)"))
    cat(sprintf("  Gauge:                   %s\n",
                inf$gauge %||% "left-orthogonal"))
  } else {
    cat("  Not prepared (call vcov(fit) or predict(..., se.fit=TRUE)\n")
    cat("  for conditional Bayesian/frequentist SE on the TT surface;\n")
    cat("  Gaussian/Poisson/Bernoulli). Level-1 only: conditional on TT rank\n")
    cat("  and fitted lambda; pointwise intervals (not simultaneous bands).\n")
  }
  cat(sprintf("\nConverged:              %s\n", x$converged))
  cat(sprintf("ALS sweeps:             %s\n",
              if (is.na(x$n_sweeps)) "NA" else as.character(x$n_sweeps)))
  cat(sprintf("PIRLS iterations:       %s\n", .fmt_pirls_iters(x)))
  cat(sprintf("Optimizer iterations:   %s\n",
              if (is.null(x$n_opt_iter) || is.na(x$n_opt_iter)) {
                "NA"
              } else {
                sprintf("%s (cumulative)", as.character(x$n_opt_iter))
              }))
  cat(sprintf("Outer iterations:       %s\n",
              if (is.null(x$n_outer) || is.na(x$n_outer)) "NA"
              else as.character(x$n_outer)))
  cat(sprintf("Criterion evaluations:  %s\n", .fmt_criterion_evals(x)))
  cat(sprintf("Wall time (s):          %.4f\n", x$timing))
  invisible(x)
}

#' @export
fitted.ttpspline <- function(object, ...) object$fitted.values

# predict.ttpspline is defined in inference.R (supports se.fit / intervals)

#' @export
residuals.ttpspline <- function(object,
                                type = c("response", "deviance", "pearson"),
                                ...) {
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  if (is.null(y)) y <- mu + object$residuals
  if (identical(type, "response")) return(object$residuals)
  key <- object$family_key
  if (identical(type, "pearson")) {
    if (identical(key, "gaussian")) return(object$residuals)
    if (identical(key, "bernoulli")) {
      return((y - mu) / sqrt(pmax(mu * (1 - mu), 1e-8)))
    }
    if (identical(key, "poisson")) {
      return((y - mu) / sqrt(pmax(mu, 1e-8)))
    }
  }
  if (identical(type, "deviance")) {
    # signed deviance residuals (GLM style, approximate)
    if (identical(key, "gaussian")) return(object$residuals)
    if (identical(key, "bernoulli")) {
      mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
      d <- 2 * (y * log(pmax(y, 1e-12) / mu) + (1 - y) * log(pmax(1 - y, 1e-12) / (1 - mu)))
      d[y == 0] <- 2 * log(1 / (1 - mu[y == 0]))
      d[y == 1] <- 2 * log(1 / mu[y == 1])
      return(sign(y - mu) * sqrt(pmax(d, 0)))
    }
    if (identical(key, "poisson")) {
      mu <- pmax(mu, 1e-12)
      d <- 2 * (ifelse(y > 0, y * log(y / mu), 0) - (y - mu))
      return(sign(y - mu) * sqrt(pmax(d, 0)))
    }
  }
  object$residuals
}

#' @export
coef.ttpspline <- function(object, full = FALSE, ...) {
  if (isTRUE(full)) {
    Theta <- tt_full_theta(object$cores)
    return(list(
      intercept = object$intercept,
      beta = object$beta %||% numeric(0),
      Theta = Theta,
      cores = object$cores
    ))
  }
  list(
    intercept = object$intercept,
    beta = object$beta %||% numeric(0),
    cores = object$cores,
    rank = object$rank,
    lambda = object$lambda
  )
}

#' @export
deviance.ttpspline <- function(object, ...) object$deviance

#' @export
plot.ttpspline <- function(x,
                           type = c("slice", "convergence"),
                           dims = c(1, 2),
                           at = "median",
                           n_grid = 40,
                           ...) {
  type <- match.arg(type)
  if (identical(type, "convergence")) {
    h <- x$history
    if (is.data.frame(h) && "deviance" %in% names(h)) {
      plot(h$pirls, h$deviance, type = "b", xlab = "PIRLS", ylab = "deviance",
           main = "PIRLS path", ...)
    } else if (is.list(h) && length(h) && !is.null(h[[1]]$rss)) {
      rss <- vapply(h, function(z) z$rss, numeric(1))
      plot(seq_along(rss), rss, type = "b", xlab = "ALS sweep", ylab = "RSS",
           main = "ALS path", ...)
    } else {
      message("No convergence history to plot.")
    }
    return(invisible(x))
  }

  dims <- as.integer(dims)
  stopifnot(length(dims) == 2L, all(dims >= 1L), all(dims <= x$d))
  # Build evaluation grid on dims; fix others
  g1 <- seq(x$x_range[1, dims[1]], x$x_range[2, dims[1]], length.out = n_grid)
  g2 <- seq(x$x_range[1, dims[2]], x$x_range[2, dims[2]], length.out = n_grid)
  fixed <- numeric(x$d)
  for (j in seq_len(x$d)) {
    if (j %in% dims) next
    if (identical(at, "median")) {
      fixed[j] <- mean(x$x_range[, j]) # midpoint of training range
    } else if (is.numeric(at) && length(at) == x$d) {
      fixed[j] <- at[j]
    } else {
      fixed[j] <- mean(x$x_range[, j])
    }
  }
  Xg <- as.matrix(expand.grid(g1, g2))
  Xfull <- matrix(0, nrow(Xg), x$d)
  Xfull[, dims[1]] <- Xg[, 1]
  Xfull[, dims[2]] <- Xg[, 2]
  for (j in seq_len(x$d)) {
    if (!j %in% dims) Xfull[, j] <- fixed[j]
  }
  eta <- {
    if (!is.null(x$linear) && length(x$beta %||% numeric(0)) > 0L) {
      lin_g <- matrix(colMeans(x$linear), nrow(Xfull), ncol(x$linear),
                      byrow = TRUE)
      colnames(lin_g) <- colnames(x$linear)
      predict(x, Xfull, type = "link", linear = lin_g)
    } else {
      predict(x, Xfull, type = "link")
    }
  }
  z <- matrix(eta, n_grid, n_grid)
  image(g1, g2, z,
        xlab = paste0("x", dims[1]), ylab = paste0("x", dims[2]),
        main = sprintf("TT slice (dims %d,%d)", dims[1], dims[2]), ...)
  invisible(x)
}

#' Human-readable TT rank label (uniform r, or anisotropic).
#' @keywords internal
#' @noRd
.tt_rank_label <- function(x) {
  ranks <- x$rank
  if (is.null(ranks) || length(ranks) < 3L) return("?")
  internal <- ranks[-c(1L, length(ranks))]
  if (length(unique(internal)) == 1L) {
    as.character(internal[1])
  } else {
    paste0("anisotropic (", paste(internal, collapse = "-"), ")")
  }
}

#' Format PIRLS iteration count for summary (NA = not used by this solver).
#' @keywords internal
#' @noRd
.fmt_pirls_iters <- function(x) {
  n <- x$n_pirls
  if (!is.null(n) && !is.na(n)) return(as.character(n))
  opt <- x$optimizer_used %||% x$optimizer %||% ""
  key <- x$family_key %||% ""
  if (identical(key, "gaussian") || identical(opt, "ALS")) {
    return("NA (not used; Gaussian ALS)")
  }
  if (opt %in% c("LBFGS", "GD", "Adam", "Damped-Newton-ALS", "LBFGS-ALS")) {
    return(sprintf("NA (not used; %s)", opt))
  }
  "NA (not used)"
}

#' Format smoothing-criterion evaluation count for summary.
#'
#' With fixed λ there is no GCV/FS/REML search; ALS may still have counted
#' one linear solve per core update in older builds. Never label those as
#' "lambda evaluations".
#'
#' @keywords internal
#' @noRd
.fmt_criterion_evals <- function(x) {
  meth <- x$lambda_method %||% "fixed"
  if (identical(meth, "fixed")) {
    return("NA (fixed lambda; not searched)")
  }
  n <- x$n_criterion_evals
  if (is.null(n) || (length(n) == 1L && is.na(n))) return("NA")
  as.character(as.integer(n))
}

