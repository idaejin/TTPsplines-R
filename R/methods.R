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
    cat(sprintf("  EDF: %.2f (!= N_TT; %.0f%% of stored params)\n",
                x$edf, 100 * x$edf / max(x$npar_tt, 1)))
  }
  cat(sprintf("  Deviance/RSS: %.6g | time: %.3fs\n",
              x$deviance, x$timing))
  invisible(x)
}

#' @export
summary.ttpspline <- function(object, ...) {
  structure(object, class = c("summary.ttpspline", "ttpspline"))
}

#' @export
print.summary.ttpspline <- function(x, ...) {
  cat("Tensor-Train P-spline fit\n\n")
  cat(sprintf("Family:                 %s\n", x$family$family))
  cat(sprintf("Link:                   %s\n", x$family$link))
  cat(sprintf("Observations:           %d\n", x$n))
  cat(sprintf("Dimensions d:           %d\n", x$d))
  cat(sprintf("Basis size k:           %d\n", x$k))
  cat(sprintf("Degree:                 %d\n", x$degree))
  cat(sprintf("Penalty order:          %d\n", x$penalty_order))
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
  cat("\n")
  cat(sprintf("Deviance / RSS:         %.6g\n", x$deviance))
  if (is.finite(x$edf)) {
    cat(sprintf("EDF (joint linearized): %.2f\n", x$edf))
    cat(sprintf("EDF / TT stored:        %.2f\n", x$edf / max(x$npar_tt, 1)))
  } else {
    note <- x$edf_note %||% "not computed"
    cat(sprintf("EDF:                    NA (%s)\n", note))
  }
  cat(sprintf("Converged:              %s\n", x$converged))
  cat(sprintf("ALS sweeps:             %s\n",
              if (is.na(x$n_sweeps)) "NA" else as.character(x$n_sweeps)))
  cat(sprintf("PIRLS iterations:       %s\n", .fmt_pirls_iters(x)))
  cat(sprintf("Optimizer iterations:   %s\n",
              if (is.null(x$n_opt_iter) || is.na(x$n_opt_iter)) "NA"
              else as.character(x$n_opt_iter)))
  cat(sprintf("Outer iterations:       %s\n",
              if (is.null(x$n_outer) || is.na(x$n_outer)) "NA"
              else as.character(x$n_outer)))
  cat(sprintf("Lambda evaluations:     %s\n",
              if (is.null(x$n_criterion_evals)) "NA" else as.character(x$n_criterion_evals)))
  cat(sprintf("Wall time (s):          %.4f\n", x$timing))
  invisible(x)
}

#' @export
fitted.ttpspline <- function(object, ...) object$fitted.values

#' @export
predict.ttpspline <- function(object,
                              newdata = NULL,
                              type = c("link", "response"),
                              ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    Xnew <- as.matrix(newdata)
    if (ncol(Xnew) != object$d) {
      stop("newdata must have ", object$d, " columns.", call. = FALSE)
    }
    if (anyNA(Xnew)) stop("NA in newdata not supported.", call. = FALSE)
    basis <- eval_marginal_bases(Xnew, object$knots, object$degree)
    eta <- object$intercept + tt_contraction(object$cores, basis)
  }
  if (identical(type, "link") || identical(object$family_key, "gaussian")) {
    return(eta)
  }
  invlink_eta(object$family, eta)
}

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
    return(list(intercept = object$intercept, Theta = Theta, cores = object$cores))
  }
  list(
    intercept = object$intercept,
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
  eta <- predict(x, Xfull, type = "link")
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

