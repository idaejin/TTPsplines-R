#' Data-driven TT rank selection via K-fold CV (+ optional 1-SE rule).
#'
#' Chooses a uniform TT rank \(r\) (chain `(1,r,...,r,1)`) using out-of-sample
#' predictive loss only — never simulation truth. Fitting stays separate:
#' [ttpspline()] still requires an explicit `rank`; use [tt_rank_refit()] after
#' selection to fit on all data at the chosen rank.
#'
#' Inner smoothing (`lambda` fixed or `"cGCV"`) is estimated **only** on each
#' fold's training portion (no validation leakage into cGCV).
#'
#' @param y,X Response and covariate matrix (as in [ttpspline()]).
#' @param ranks Integer vector of candidate **uniform** ranks (default `1:5`).
#' @param family,k,degree,penalty_order,lambda,optimizer,backend,control,knots,offset,weights
#'   Passed through to [ttpspline()] on each fold (and on final refit).
#' @param folds Number of CV folds (default 5).
#' @param rule `"1se"` (default, parsimonious) or `"min"` (minimum mean CV).
#' @param metric `"auto"` (family-aware) or one of `"rmse"`, `"poisson_deviance"`,
#'   `"logloss"`.
#' @param seed Optional RNG seed for fold assignment (reproducible splits).
#' @param keep_fits If `TRUE`, store fold fits (large); default `FALSE`.
#' @param rank_chain Reserved for future non-uniform rank search; must be `NULL`
#'   in this version.
#' @param ... Currently unused (reserved).
#'
#' @return An object of class `"tt_rank_selection"`.
#'
#' @section Complexity layers:
#' Rank \(r\) = structural capacity; \(\lambda\) = directional smoothness;
#' EDF = effective fitted flexibility. These are not interchangeable:
#' \(r \neq \lambda \neq \mathrm{EDF}\).
#'
#' @seealso [tt_rank_refit()], [ttpspline()], [tt_rank_profile()] (in-sample
#'   rank diagnostic, not CV).
#'
#' @examples
#' data(friedman)
#' X <- as.matrix(friedman[1:200, paste0("x", 1:5)])
#' y <- friedman$y[1:200]
#' sel <- tt_rank_select(
#'   y, X, ranks = 1:3, k = 5, lambda = 1, folds = 3, rule = "1se",
#'   seed = 1, control = tt_control(max_sweeps = 4, compute_edf = FALSE)
#' )
#' sel
#' fit <- tt_rank_refit(sel)
#'
#' @export
tt_rank_select <- function(y,
                           X,
                           ranks = 1:5,
                           family = stats::gaussian(),
                           k = 10,
                           degree = 3,
                           penalty_order = 2,
                           lambda = "cGCV",
                           optimizer = "auto",
                           backend = "auto",
                           control = tt_control(),
                           folds = 5L,
                           rule = c("1se", "min"),
                           metric = c("auto", "rmse", "poisson_deviance", "logloss"),
                           seed = NULL,
                           keep_fits = FALSE,
                           knots = NULL,
                           offset = NULL,
                           weights = NULL,
                           rank_chain = NULL,
                           ...) {
  cl <- match.call()
  if (!is.null(rank_chain)) {
    stop(
      "`rank_chain` is reserved for a future non-uniform search; ",
      "pass uniform candidates via `ranks` for now.",
      call. = FALSE
    )
  }
  dots <- list(...)
  if (length(dots)) {
    stop("Unused arguments in tt_rank_select(): ",
         paste(names(dots), collapse = ", "), call. = FALSE)
  }

  rule <- match.arg(rule)
  metric_arg <- match.arg(metric)
  fam <- normalize_family(family)
  key <- family_key(fam)
  metric <- .tt_resolve_cv_metric(metric_arg, key)

  y <- as.numeric(y)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- length(y)
  if (nrow(X) != n) stop("`y` and `X` must have the same number of rows.", call. = FALSE)
  if (anyNA(y) || anyNA(X)) stop("NA in y/X not supported.", call. = FALSE)

  ranks <- sort(unique(as.integer(ranks)))
  if (!length(ranks) || any(ranks < 1L)) {
    stop("`ranks` must be positive integers.", call. = FALSE)
  }
  folds <- as.integer(folds)[1L]
  if (folds < 2L || folds > n) {
    stop("`folds` must satisfy 2 <= folds <= n.", call. = FALSE)
  }

  offset_full <- normalize_offset(offset, n)
  weights_full <- normalize_weights(weights, n)
  fold_id <- .tt_make_fold_id(n, folds, seed)

  # CV fits: skip EDF (costly / irrelevant for selection)
  ctrl <- control
  if (is.null(ctrl$compute_edf) || isTRUE(ctrl$compute_edf)) {
    ctrl$compute_edf <- FALSE
  }

  fit_args <- list(
    family = fam,
    k = k,
    degree = degree,
    penalty_order = penalty_order,
    lambda = lambda,
    optimizer = optimizer,
    backend = backend,
    control = ctrl,
    knots = knots
  )
  lambda_method <- if (is.character(lambda) && identical(lambda, "cGCV")) {
    "cGCV"
  } else {
    "fixed"
  }

  d <- ncol(X)
  p <- as.integer(k)[1L]
  loss_mat <- matrix(NA_real_, nrow = length(ranks), ncol = folds,
                     dimnames = list(as.character(ranks), paste0("fold", seq_len(folds))))
  time_mat <- matrix(NA_real_, nrow = length(ranks), ncol = folds,
                     dimnames = dimnames(loss_mat))
  conv_mat <- matrix(FALSE, nrow = length(ranks), ncol = folds,
                     dimnames = dimnames(loss_mat))
  lambda_list <- vector("list", length(ranks))
  names(lambda_list) <- as.character(ranks)
  fold_fits <- if (isTRUE(keep_fits)) {
    vector("list", length(ranks))
  } else {
    NULL
  }
  if (isTRUE(keep_fits)) names(fold_fits) <- as.character(ranks)

  for (i in seq_along(ranks)) {
    r <- ranks[i]
    lambda_list[[i]] <- vector("list", folds)
    if (isTRUE(keep_fits)) fold_fits[[i]] <- vector("list", folds)
    for (f in seq_len(folds)) {
      test <- fold_id == f
      train <- !test
      off_tr <- offset_full[train]
      off_te <- offset_full[test]
      w_tr <- weights_full[train]
      w_te <- weights_full[test]
      t0 <- proc.time()[["elapsed"]]
      fit <- tryCatch(
        ttpspline(
          y[train], X[train, , drop = FALSE],
          family = fam,
          rank = r,
          k = k,
          degree = degree,
          penalty_order = penalty_order,
          lambda = lambda,
          optimizer = optimizer,
          backend = backend,
          control = ctrl,
          knots = knots,
          offset = off_tr,
          weights = w_tr
        ),
        error = function(e) e
      )
      elapsed <- proc.time()[["elapsed"]] - t0
      time_mat[i, f] <- elapsed
      if (inherits(fit, "error")) {
        loss_mat[i, f] <- Inf
        conv_mat[i, f] <- FALSE
        lambda_list[[i]][[f]] <- NA_real_
        next
      }
      mu <- tryCatch(
        predict(fit, newdata = X[test, , drop = FALSE],
                type = "response", offset = off_te),
        error = function(e) NULL
      )
      if (is.null(mu) || anyNA(mu) || !all(is.finite(mu))) {
        loss_mat[i, f] <- Inf
        conv_mat[i, f] <- FALSE
        lambda_list[[i]][[f]] <- as.numeric(fit$lambda)
        next
      }
      loss_mat[i, f] <- .tt_cv_loss(y[test], mu, metric = metric, family = fam,
                                   weights = w_te)
      conv_mat[i, f] <- isTRUE(fit$converged)
      lambda_list[[i]][[f]] <- as.numeric(fit$lambda)
      if (isTRUE(keep_fits)) fold_fits[[i]][[f]] <- fit
    }
  }

  cv_results <- .tt_summarize_rank_cv(
    ranks = ranks,
    loss_mat = loss_mat,
    time_mat = time_mat,
    conv_mat = conv_mat,
    lambda_list = lambda_list,
    lambda_method = lambda_method,
    d = d,
    p = p
  )

  rank_min <- .tt_rank_min_cv(cv_results)
  rank_1se <- .tt_rank_1se(cv_results, rank_min)
  selected <- if (identical(rule, "min")) rank_min else rank_1se

  structure(
    list(
      ranks = ranks,
      metric = metric,
      metric_requested = metric_arg,
      folds = folds,
      fold_id = fold_id,
      cv_results = cv_results,
      loss_by_fold = loss_mat,
      time_by_fold = time_mat,
      converged_by_fold = conv_mat,
      lambda_by_fold = lambda_list,
      rank_min = rank_min,
      rank_1se = rank_1se,
      selected_rank = selected,
      rule = rule,
      family = fam,
      family_key = key,
      lambda = lambda,
      lambda_method = lambda_method,
      fit_args = fit_args,
      y = y,
      X = X,
      offset = offset_full,
      weights = weights_full,
      call = cl,
      timings = list(
        total_s = sum(time_mat, na.rm = TRUE),
        by_rank = rowSums(time_mat, na.rm = TRUE)
      ),
      keep_fits = isTRUE(keep_fits),
      fits = fold_fits
    ),
    class = "tt_rank_selection"
  )
}

#' Refit [ttpspline()] on all data at the selected TT rank.
#'
#' Uses the fitting arguments stored on a `"tt_rank_selection"` object.
#' Does **not** change [ttpspline()] defaults: this is the explicit final fit
#' after model selection.
#'
#' @param object A `"tt_rank_selection"` from [tt_rank_select()].
#' @param y,X Optional override data (default: data stored on `object`).
#' @param rank Rank to refit (default: `object$selected_rank`).
#' @param control Optional [tt_control()] override (default: selection control,
#'   but `compute_edf` restored to the stored value / `TRUE` if you pass a new
#'   control).
#' @param ... Passed to [ttpspline()] (override stored args).
#' @return A `"ttpspline"` fit.
#' @export
tt_rank_refit <- function(object,
                          y = NULL,
                          X = NULL,
                          rank = NULL,
                          control = NULL,
                          ...) {
  stopifnot(inherits(object, "tt_rank_selection"))
  if (is.null(y)) y <- object$y
  if (is.null(X)) X <- object$X
  if (is.null(y) || is.null(X)) {
    stop("y/X not stored; pass them explicitly to tt_rank_refit().", call. = FALSE)
  }
  if (is.null(rank)) rank <- object$selected_rank
  args <- object$fit_args
  if (!is.null(control)) args$control <- control
  # Refit may want EDF; if still FALSE from CV, leave as stored unless overridden
  extra <- list(...)
  args$y <- y
  args$X <- X
  args$rank <- rank
  args$offset <- if ("offset" %in% names(extra)) {
    extra$offset
  } else {
    object$offset
  }
  extra$offset <- NULL
  args$weights <- if ("weights" %in% names(extra)) {
    extra$weights
  } else {
    object$weights
  }
  extra$weights <- NULL
  if (length(extra)) args <- utils::modifyList(args, extra)
  do.call(ttpspline, args)
}

# ---- internals -------------------------------------------------------------

#' @keywords internal
#' @noRd
.tt_make_fold_id <- function(n, folds, seed = NULL) {
  if (!is.null(seed)) {
    # isolate RNG state
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(as.integer(seed))
  }
  sample(rep(seq_len(as.integer(folds)), length.out = as.integer(n)))
}

#' @keywords internal
#' @noRd
.tt_resolve_cv_metric <- function(metric, family_key) {
  if (!identical(metric, "auto")) return(metric)
  if (identical(family_key, "gaussian")) return("rmse")
  if (identical(family_key, "poisson")) return("poisson_deviance")
  if (identical(family_key, "bernoulli")) return("logloss")
  stop("No default CV metric for family_key=", family_key, call. = FALSE)
}

#' Mean validation loss for one fold (lower is better).
#' @keywords internal
#' @noRd
.tt_cv_loss <- function(y, mu, metric, family, weights = NULL) {
  y <- as.numeric(y)
  mu <- as.numeric(mu)
  n <- length(y)
  w <- normalize_weights(weights, n)
  sw <- sum(w)
  if (identical(metric, "rmse")) {
    return(sqrt(sum(w * (y - mu)^2) / sw))
  }
  if (identical(metric, "poisson_deviance")) {
    return(as.numeric(glm_deviance(family, y, mu, weights = w) / sw))
  }
  if (identical(metric, "logloss")) {
    mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
    return(-sum(w * (y * log(mu) + (1 - y) * log(1 - mu))) / sw)
  }
  stop("Unknown metric: ", metric, call. = FALSE)
}

#' @keywords internal
#' @noRd
.tt_summarize_rank_cv <- function(ranks, loss_mat, time_mat, conv_mat,
                                  lambda_list, lambda_method, d, p) {
  rows <- lapply(seq_along(ranks), function(i) {
    r <- ranks[i]
    losses <- as.numeric(loss_mat[i, ])
    times <- as.numeric(time_mat[i, ])
    conv <- as.logical(conv_mat[i, ])
    finite <- is.finite(losses)
    K <- sum(finite)
    mean_cv <- if (K > 0) mean(losses[finite]) else Inf
    sd_cv <- if (K > 1) stats::sd(losses[finite]) else NA_real_
    se_cv <- if (K > 1) sd_cv / sqrt(K) else NA_real_
    chain <- tt_rank(r, d = d)
    npar_tt <- as.integer(tt_npar(p, chain))
    npar_dense <- as.integer(dense_npar(p, d))
    npar_intr <- as.integer(npar_tt - tt_gauge_dim(chain))
    lam_fold <- lambda_list[[i]]
    lam_mean <- NA_real_
    if (identical(lambda_method, "cGCV")) {
      lam_mat <- do.call(rbind, lapply(lam_fold, function(v) {
        if (is.null(v) || all(is.na(v))) return(rep(NA_real_, d))
        rep(as.numeric(v), length.out = d)
      }))
      lam_mean <- colMeans(lam_mat, na.rm = TRUE)
    }
    data.frame(
      rank = r,
      mean_cv = mean_cv,
      sd_cv = sd_cv,
      se_cv = se_cv,
      npar_tt = npar_tt,
      npar_dense = npar_dense,
      intrinsic_dim = npar_intr,
      compression = npar_dense / max(npar_tt, 1),
      convergence_rate = mean(conv),
      mean_time = mean(times, na.rm = TRUE),
      total_time = sum(times, na.rm = TRUE),
      n_finite_folds = K,
      lambda_mean = if (length(lam_mean) == 1L) lam_mean else mean(lam_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.tt_rank_min_cv <- function(cv_results, tol = 1e-12) {
  m <- cv_results$mean_cv
  best <- min(m)
  # tie → smallest rank (table is sorted by ranks)
  cv_results$rank[which(m <= best + tol)[1L]]
}

#' Smallest rank with mean_cv <= mean_cv(r_min) + SE(r_min).
#' @keywords internal
#' @noRd
.tt_rank_1se <- function(cv_results, rank_min, tol = 1e-12) {
  row_min <- cv_results[cv_results$rank == rank_min, , drop = FALSE]
  thresh <- row_min$mean_cv + (row_min$se_cv %||% 0)
  if (!is.finite(thresh)) return(rank_min)
  ok <- cv_results$mean_cv <= thresh + tol
  if (!any(ok)) return(rank_min)
  min(cv_results$rank[ok])
}

#' @export
print.tt_rank_selection <- function(x, digits = 3, ...) {
  cat("Tensor-Train rank selection\n\n")
  cat(sprintf("Family:             %s\n", x$family$family))
  cat(sprintf("Metric:             %s\n", x$metric))
  cat(sprintf("Folds:              %d\n", x$folds))
  lam_lab <- if (identical(x$lambda_method, "cGCV")) {
    "cGCV"
  } else if (is.numeric(x$lambda)) {
    paste(sprintf("%.4g", x$lambda), collapse = ",")
  } else {
    as.character(x$lambda)
  }
  cat(sprintf("Lambda:             %s\n\n", lam_lab))

  tab <- x$cv_results
  show <- data.frame(
    rank = tab$rank,
    `CV error` = round(tab$mean_cv, digits),
    SE = round(tab$se_cv, digits),
    `TT params` = tab$npar_tt,
    Compression = sprintf("%.1fx", tab$compression),
    check.names = FALSE
  )
  print(show, row.names = FALSE)
  cat("\n")
  cat(sprintf("Minimum-CV rank:    %d\n", x$rank_min))
  cat(sprintf("1-SE rank:          %d\n", x$rank_1se))
  cat(sprintf("Selected rank:      %d\n", x$selected_rank))
  cat(sprintf("Rule:               %s\n",
              if (identical(x$rule, "1se")) "1-SE" else "minimum-CV"))
  invisible(x)
}

#' @export
summary.tt_rank_selection <- function(object, ...) {
  structure(object, class = c("summary.tt_rank_selection", "tt_rank_selection"))
}

#' @export
print.summary.tt_rank_selection <- function(x, digits = 4, ...) {
  print.tt_rank_selection(x, digits = digits)
  tab <- x$cv_results
  cat("\nComplexity / cost\n")
  show <- data.frame(
    rank = tab$rank,
    npar_tt = tab$npar_tt,
    intrinsic = tab$intrinsic_dim,
    CR = round(tab$compression, 1),
    conv = round(tab$convergence_rate, 2),
    mean_t = round(tab$mean_time, 3),
    total_t = round(tab$total_time, 3),
    stringsAsFactors = FALSE
  )
  print(show, row.names = FALSE)
  cat(sprintf("\nTotal CV wall time: %.3fs\n", x$timings$total_s))
  cat("SE is fold-to-fold variability of the CV loss, not a formal rank test.\n")
  if (identical(x$lambda_method, "cGCV")) {
    cat("Lambda was selected by cGCV on training folds only (no validation leakage).\n")
    cat(sprintf("Mean training lambda (avg over dims/folds): %s\n",
                paste(sprintf("%.4g", tab$lambda_mean), collapse = ", ")))
  }
  invisible(x)
}

#' @export
plot.tt_rank_selection <- function(x,
                                   type = c("error", "compression", "tradeoff"),
                                   ...) {
  type <- match.arg(type)
  tab <- x$cv_results
  r <- tab$rank
  if (identical(type, "error")) {
    y <- tab$mean_cv
    se <- tab$se_cv
    se[!is.finite(se)] <- 0
    ylim <- range(c(y - se, y + se), finite = TRUE)
    plot(r, y, type = "b", pch = 19, xlab = "TT rank r",
         ylab = sprintf("CV %s", x$metric),
         main = "CV error vs rank", ylim = ylim, ...)
    if (any(is.finite(se) & se > 0)) {
      graphics::arrows(r, y - se, r, y + se, length = 0.05, angle = 90, code = 3)
    }
    graphics::abline(v = x$rank_min, col = "gray40", lty = 2)
    graphics::abline(v = x$rank_1se, col = "darkred", lty = 3)
    graphics::legend("topright",
                     legend = c("mean CV", "min-CV rank", "1-SE rank"),
                     lty = c(1, 2, 3), pch = c(19, NA, NA),
                     col = c("black", "gray40", "darkred"), bty = "n", cex = 0.8)
  } else if (identical(type, "compression")) {
    plot(r, tab$compression, type = "b", pch = 19,
         xlab = "TT rank r", ylab = "Compression (dense / TT)",
         main = "Compression vs rank", ...)
    graphics::abline(v = x$selected_rank, col = "darkred", lty = 3)
  } else {
    plot(tab$npar_tt, tab$mean_cv, type = "b", pch = 19,
         xlab = "TT stored parameters", ylab = sprintf("CV %s", x$metric),
         main = "Prediction vs complexity", ...)
    idx <- match(x$selected_rank, r)
    if (!is.na(idx)) {
      graphics::points(tab$npar_tt[idx], tab$mean_cv[idx], pch = 19, col = "darkred", cex = 1.4)
    }
  }
  invisible(x)
}
