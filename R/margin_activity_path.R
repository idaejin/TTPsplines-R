#' Margin Activity Path: score and optionally select TT margins
#'
#' Screens continuous covariates for a Tensor-Train P-spline by measuring
#' **partial-range activity** along an isotropic \(\lambda\) path, ranking
#' margins, and (optionally) choosing the nested top-\(m\) subset by \(K\)-fold
#' CV with a one-standard-error or min-CV rule. Final smoothing still uses
#' [ttps()] with `lambda = "cGCV"` (or a user choice) on the selected columns.
#'
#' The method name is **Margin Activity Path** (avoid the acronym "MAP", which
#' collides with *maximum a posteriori*). Protocol:
#' `inst/benchmarks/margin_path/PROTOCOL_MAP.md`.
#'
#' @param y,X Response and \(n\times d\) covariate matrix (as in [ttps()]).
#' @param lambda_path Numeric vector of isotropic \(\lambda\) values for the
#'   screening path (default `10^seq(2, -1, length.out = 7)`).
#' @param rank,k,degree,penalty_order,family,optimizer,backend,control,offset,weights
#'   Passed to [ttps()] on screening fits (and on nested CV / refit when
#'   applicable). For subsets with \(d=1\), a univariate `smooth.spline` GCV
#'   fallback is used because [ttps()] requires \(d\ge 2\).
#' @param reference How to fix \(x_{-j}\) when forming partials: `"median"`
#'   (default) or `"mean"`.
#' @param n_grid Grid size for partial-range evaluation (default `40`).
#' @param grid Optional numeric vector used as the evaluation grid for each
#'   partial when supplied; otherwise each margin uses an interior quantile
#'   grid of length `n_grid` (or a fixed interior grid on approximately unit
#'   interval covariates).
#' @param select Selection rule after ranking: `"1se"` (default), `"min"`, or
#'   `"none"` (return scores / ranking only).
#' @param folds Number of CV folds (default `5`). Ignored if `foldid` is set.
#' @param foldid Optional fold ids (length \(n\)), as in [tt_rank_select()].
#' @param select_lambda Smoothing for nested CV and final refit (default
#'   `"cGCV"`).
#' @param refit If `TRUE` (default) and `select != "none"`, refit [ttps()] on
#'   the selected margins with `select_lambda`.
#' @param seed Optional RNG seed for fold assignment when `foldid` is `NULL`.
#' @param verbose Print progress messages.
#'
#' @return An object of class `"tt_margin_activity_path"` with components
#'   `scores`, `activity`, `lambda_path`, `order`, `selected`, `m_hat`,
#'   `cv` (data frame or `NULL`), `fit` (or `NULL`), `call`, and metadata.
#'
#' @seealso [ttps()], [tt_rank_select()], [plot.tt_margin_activity_path()]
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' d <- 6
#' X <- matrix(runif(n * d), n, d)
#' colnames(X) <- paste0("x", seq_len(d))
#' y <- sin(2 * pi * X[, 1]) + sin(2 * pi * X[, 2]) + rnorm(n, sd = 0.3)
#' path <- tt_margin_activity_path(
#'   y, X, rank = 2, k = 5,
#'   lambda_path = 10^c(1, 0, -1),
#'   select = "1se", folds = 3, seed = 1,
#'   control = tt_control(max_sweeps = 4, compute_edf = FALSE, outer_maxit = 4),
#'   verbose = FALSE
#' )
#' path
#' path$selected
#'
#' @export
tt_margin_activity_path <- function(y,
                                    X,
                                    lambda_path = 10^seq(2, -1, length.out = 7),
                                    rank = 2L,
                                    k = 5L,
                                    degree = 3L,
                                    penalty_order = 2L,
                                    family = stats::gaussian(),
                                    optimizer = "auto",
                                    backend = "auto",
                                    control = tt_control(),
                                    offset = NULL,
                                    weights = NULL,
                                    reference = c("median", "mean"),
                                    n_grid = 40L,
                                    grid = NULL,
                                    select = c("1se", "min", "none"),
                                    folds = 5L,
                                    foldid = NULL,
                                    select_lambda = "cGCV",
                                    refit = TRUE,
                                    seed = NULL,
                                    verbose = TRUE) {
  reference <- match.arg(reference)
  select <- match.arg(select)
  X <- as.matrix(X)
  y <- as.numeric(y)
  stopifnot(length(y) == nrow(X), ncol(X) >= 1L)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  }
  d <- ncol(X)
  n <- nrow(X)
  lambda_path <- as.numeric(lambda_path)
  stopifnot(length(lambda_path) >= 1L, all(is.finite(lambda_path)),
            all(lambda_path > 0))

  ctrl_screen <- control
  if (is.null(ctrl_screen$compute_edf)) {
    ctrl_screen$compute_edf <- FALSE
  }

  # ---- activity path -------------------------------------------------------
  activity <- matrix(NA_real_, length(lambda_path), d)
  colnames(activity) <- colnames(X)
  rownames(activity) <- paste0("lambda=", format(lambda_path, digits = 4))

  if (isTRUE(verbose)) {
    message("Margin Activity Path: isotropic screening (", length(lambda_path),
            " lambda values, d = ", d, ")")
  }
  for (i in seq_along(lambda_path)) {
    fit_i <- .tt_map_fit(
      y, X,
      rank = rank, k = k, degree = degree, penalty_order = penalty_order,
      lambda = lambda_path[[i]], family = family, optimizer = optimizer,
      backend = backend, control = ctrl_screen, offset = offset, weights = weights
    )
    for (j in seq_len(d)) {
      activity[i, j] <- .tt_partial_range(
        fit_i, X, j,
        reference = reference, n_grid = n_grid, grid = grid
      )
    }
    if (isTRUE(verbose)) {
      message(sprintf("  lambda = %.4g done", lambda_path[[i]]))
    }
  }

  scores <- colMeans(activity)
  ord <- order(scores, decreasing = TRUE)
  names(scores) <- colnames(X)

  out <- list(
    scores = scores,
    activity = activity,
    lambda_path = lambda_path,
    order = ord,
    selected = integer(0),
    selected_names = character(0),
    m_hat = NA_integer_,
    m_gap = .tt_map_gap_m(scores[ord]),
    cv = NULL,
    fit = NULL,
    rule = select,
    reference = reference,
    rank = as.integer(rank),
    k = as.integer(k),
    select_lambda = select_lambda,
    call = match.call(),
    variable_names = colnames(X)
  )
  class(out) <- "tt_margin_activity_path"

  if (identical(select, "none")) {
    return(out)
  }

  # ---- nested CV on top-m --------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  if (is.null(foldid)) {
    folds <- as.integer(folds)
    stopifnot(folds >= 2L, folds <= n)
    foldid <- sample(rep(seq_len(folds), length.out = n))
  } else {
    foldid <- as.integer(foldid)
    stopifnot(length(foldid) == n)
  }
  K <- length(unique(foldid))

  ctrl_cv <- control
  if (is.null(ctrl_cv$compute_edf)) ctrl_cv$compute_edf <- FALSE

  ms <- 0:d
  cv_mse <- cv_se <- numeric(length(ms))
  margins_lab <- character(length(ms))

  if (isTRUE(verbose)) {
    message("Nested CV over top-m models (K = ", K, ", rule = ", select, ")")
  }
  for (ii in seq_along(ms)) {
    m <- ms[[ii]]
    cols <- if (m == 0L) integer(0) else ord[seq_len(m)]
    res <- .tt_map_cv_mse(
      y, X, cols, foldid,
      rank = rank, k = k, degree = degree, penalty_order = penalty_order,
      lambda = select_lambda, family = family, optimizer = optimizer,
      backend = backend, control = ctrl_cv, offset = offset, weights = weights
    )
    cv_mse[[ii]] <- res$mse
    cv_se[[ii]] <- res$se
    margins_lab[[ii]] <- if (m == 0L) "(none)" else
      paste(colnames(X)[cols], collapse = ",")
    if (isTRUE(verbose)) {
      message(sprintf("  m=%2d  CV-MSE=%.4f  (SE=%.4f)", m, res$mse, res$se))
    }
  }

  cv_tab <- data.frame(
    m = ms,
    cv_mse = cv_mse,
    cv_se = cv_se,
    margins = margins_lab,
    stringsAsFactors = FALSE
  )
  m_min <- ms[[which.min(cv_mse)]]
  cv_star <- cv_mse[[which.min(cv_mse)]]
  se_star <- cv_se[[which.min(cv_mse)]]

  if (identical(select, "1se")) {
    ok <- which(cv_mse <= cv_star + se_star)
    m_hat <- min(ms[ok])
  } else {
    m_hat <- m_min
  }

  selected <- if (m_hat == 0L) integer(0) else ord[seq_len(m_hat)]
  out$cv <- cv_tab
  out$m_hat <- as.integer(m_hat)
  out$m_min <- as.integer(m_min)
  out$selected <- as.integer(selected)
  out$selected_names <- colnames(X)[selected]
  out$foldid <- foldid

  if (isTRUE(refit) && length(selected) >= 1L) {
    if (isTRUE(verbose)) {
      message("Refit on selected margins: ",
              paste(out$selected_names, collapse = ", "))
    }
    out$fit <- .tt_map_fit(
      y, X[, selected, drop = FALSE],
      rank = rank, k = k, degree = degree, penalty_order = penalty_order,
      lambda = select_lambda, family = family, optimizer = optimizer,
      backend = backend, control = ctrl_cv, offset = offset, weights = weights
    )
  } else if (isTRUE(refit) && length(selected) == 0L) {
    out$fit <- list(
      type = "intercept",
      mu = mean(y),
      predict = function(Xnew) rep(mean(y), NROW(Xnew))
    )
  }

  out
}

#' @export
print.tt_margin_activity_path <- function(x, ...) {
  cat("Margin Activity Path\n")
  cat("  d = ", length(x$scores),
      ", path length = ", length(x$lambda_path),
      ", rule = ", x$rule, "\n", sep = "")
  sc <- sort(x$scores, decreasing = TRUE)
  cat("  Scores (sorted):\n")
  print(round(sc, 4))
  if (!identical(x$rule, "none") && !is.null(x$m_hat)) {
    cat("  Selected m = ", x$m_hat, ": ",
        if (length(x$selected_names)) {
          paste(x$selected_names, collapse = ", ")
        } else {
          "(none)"
        },
        "\n", sep = "")
  }
  if (!is.null(x$fit) && inherits(x$fit, "ttps")) {
    cat("  Final fit: ttps object on selected margins\n")
  }
  invisible(x)
}

#' Plot Margin Activity Path diagnostics
#'
#' @param x A `"tt_margin_activity_path"` object.
#' @param which Which panels: `1` activity path, `2` sorted scores, `3` CV path
#'   (if available). Default all available.
#' @param ... Unused.
#' @export
plot.tt_margin_activity_path <- function(x, which = NULL, ...) {
  has_cv <- !is.null(x$cv)
  if (is.null(which)) {
    which <- if (has_cv) 1:3 else 1:2
  }
  which <- as.integer(which)
  n_panels <- length(which)
  op <- graphics::par(mfrow = c(1L, n_panels), mar = c(4.5, 4.5, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  d <- length(x$scores)
  active <- seq_len(d) %in% x$selected
  if (!length(x$selected) && identical(x$rule, "none")) {
    cols <- rep("#1F4E79", d)
  } else {
    cols <- ifelse(active, "#1F4E79", "#BBBBBB")
  }
  ord <- x$order

  if (1L %in% which) {
    graphics::matplot(
      log10(x$lambda_path), x$activity,
      type = "l", lwd = 2, lty = 1, col = cols,
      xlab = expression(log[10](lambda)),
      ylab = "partial range",
      main = "Activity path"
    )
  }
  if (2L %in% which) {
    graphics::barplot(
      x$scores[ord],
      names.arg = names(x$scores)[ord],
      col = cols[ord], border = NA, las = 2, cex.names = 0.75,
      ylab = expression(S[j]),
      main = "Sorted scores"
    )
  }
  if (3L %in% which) {
    if (!has_cv) {
      graphics::plot.new()
      graphics::title("CV path (not computed)")
    } else {
      graphics::plot(
        x$cv$m, x$cv$cv_mse, type = "b", lwd = 2, pch = 16, col = "#1F4E79",
        xlab = "m (top-m margins)", ylab = "CV-MSE",
        main = paste0("CV (", x$rule, "), m = ", x$m_hat)
      )
      graphics::arrows(
        x$cv$m, x$cv$cv_mse - x$cv$cv_se,
        x$cv$m, x$cv$cv_mse + x$cv$cv_se,
        length = 0.05, angle = 90, code = 3,
        col = grDevices::adjustcolor("#1F4E79", 0.6)
      )
      if (is.finite(x$m_hat)) {
        graphics::abline(v = x$m_hat, lty = 2, col = "#C0392B")
      }
    }
  }
  invisible(x)
}

# ---- internals -------------------------------------------------------------

.tt_map_ref <- function(X, reference) {
  if (identical(reference, "mean")) {
    apply(X, 2L, mean)
  } else {
    apply(X, 2L, stats::median)
  }
}

.tt_partial_range <- function(fit, X, j, reference = "median",
                              n_grid = 40L, grid = NULL) {
  ref <- .tt_map_ref(X, reference)
  n_grid <- as.integer(n_grid)
  xj <- X[, j]
  if (!is.null(grid)) {
    g <- as.numeric(grid)
  } else {
    # Prefer unit-interval interior grid when covariates look like [0,1]
    rng <- range(xj, finite = TRUE)
    if (is.finite(rng[1]) && is.finite(rng[2]) &&
        rng[1] >= -1e-8 && rng[2] <= 1 + 1e-8) {
      g <- seq(0.05, 0.95, length.out = n_grid)
    } else {
      g <- seq(
        stats::quantile(xj, 0.05, names = FALSE, type = 7),
        stats::quantile(xj, 0.95, names = FALSE, type = 7),
        length.out = n_grid
      )
    }
  }
  Xm <- matrix(ref, length(g), ncol(X), byrow = TRUE)
  colnames(Xm) <- colnames(X)
  Xm[, j] <- g
  pred <- .tt_map_predict(fit, Xm)
  diff(range(pred))
}

.tt_map_fit <- function(y, X, rank, k, degree, penalty_order, lambda,
                        family, optimizer, backend, control,
                        offset = NULL, weights = NULL) {
  p <- ncol(X)
  if (p < 1L) {
    mu <- mean(y)
    return(list(
      type = "intercept",
      mu = mu,
      predict = function(Xnew) rep(mu, NROW(Xnew))
    ))
  }
  if (p == 1L) {
    ss <- stats::smooth.spline(X[, 1L], y, cv = NA)
    return(list(
      type = "univariate",
      ss = ss,
      predict = function(Xnew) {
        as.numeric(stats::predict(ss, x = as.matrix(Xnew)[, 1L])$y)
      }
    ))
  }
  r_use <- max(1L, as.integer(rank))
  ttps(
    y, X,
    rank = r_use,
    k = k,
    degree = degree,
    penalty_order = penalty_order,
    lambda = lambda,
    family = family,
    optimizer = optimizer,
    backend = backend,
    control = control,
    offset = offset,
    weights = weights
  )
}

.tt_map_predict <- function(fit, Xnew) {
  if (is.list(fit) && !is.null(fit$type) &&
      fit$type %in% c("intercept", "univariate")) {
    return(as.numeric(fit$predict(Xnew)))
  }
  as.numeric(stats::predict(fit, Xnew))
}

.tt_map_cv_mse <- function(y, X, cols, foldid, ...) {
  K <- length(unique(foldid))
  pe <- numeric(K)
  for (k in seq_len(K)) {
    tr <- foldid != k
    te <- !tr
    if (length(cols) == 0L) {
      mu <- mean(y[tr])
      pe[[k]] <- mean((mu - y[te])^2)
    } else {
      fi <- .tt_map_fit(y[tr], X[tr, cols, drop = FALSE], ...)
      pe[[k]] <- mean((.tt_map_predict(fi, X[te, cols, drop = FALSE]) - y[te])^2)
    }
  }
  list(
    mse = mean(pe),
    se = if (K >= 2L) stats::sd(pe) / sqrt(K) else NA_real_,
    fold_mse = pe
  )
}

.tt_map_gap_m <- function(S_ord) {
  S_ord <- as.numeric(S_ord)
  d <- length(S_ord)
  if (d < 1L) return(0L)
  gaps <- c(S_ord[-d] - S_ord[-1L], S_ord[[d]])
  as.integer(which.max(gaps))
}
