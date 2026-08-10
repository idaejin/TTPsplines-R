#' Multi-start TT-P-spline fits (nonconvex diagnostics).
#'
#' Runs [ttps()] from several random TT initializations and returns the best
#' fit together with a start-wise summary. Useful when ALS/cGCV may land in
#' poor basins or when \(\hat\lambda\) hits search bounds: stable boundaries
#' across starts suggest a statistical margin/null/\(r\) effect; unstable ones
#' suggest initialization sensitivity.
#'
#' This does **not** change the default single-start [ttps()] path. For
#' predictive rank choice with multi-start inside CV folds, use
#' [tt_rank_select()] with `n_starts > 1`.
#'
#' @param y Numeric response.
#' @param X Covariate matrix (`n × d`).
#' @param family,rank,k,degree,penalty_order,lambda,optimizer,backend,control,monitor,knots,offset,linear,weights,cyclic,period
#'   Passed to [ttps()].
#' @param n_starts Number of random TT initializations (default `10`).
#' @param seed Base RNG seed; start `s` uses `seed + s - 1`.
#' @param select Criterion for the retained fit: `"objective"` (penalized
#'   training objective via [tt_objective()], default) or `"deviance"`.
#' @param prefer_converged If `TRUE` (default), a converged start beats a
#'   non-converged one even with a slightly worse criterion.
#' @param verbose Print a one-line progress message per start.
#'
#' @return An object of class `"ttps_multistart"` with components:
#'   \describe{
#'     \item{best}{The selected `"ttpspline"` fit.}
#'     \item{starts}{`data.frame` with one row per start.}
#'     \item{n_starts, select, seed}{Settings.}
#'     \item{frac_boundary}{Named vector: fraction of finite starts hitting
#'       each margin's search bound (`NA` if not cGCV).}
#'   }
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' X <- matrix(runif(n * 2), n, 2)
#' y <- sin(2 * pi * X[, 1]) + rnorm(n, 0, 0.2)
#' ms <- ttps_multistart(
#'   y, X, rank = 2, k = 5, lambda = 1, n_starts = 3, seed = 1,
#'   control = tt_control(max_sweeps = 4, compute_edf = FALSE)
#' )
#' ms
#' ms$best$deviance
#'
#' @export
#' @seealso [ttps()], [tt_initialize()], [tt_rank_select()]
ttps_multistart <- function(y,
                            X,
                            family = stats::gaussian(),
                            rank = 3,
                            k = 10,
                            degree = 3,
                            penalty_order = 2,
                            lambda = "cGCV",
                            optimizer = "auto",
                            backend = "auto",
                            control = tt_control(),
                            monitor = FALSE,
                            knots = NULL,
                            offset = NULL,
                            linear = NULL,
                            weights = NULL,
                            cyclic = NULL,
                            period = NULL,
                            n_starts = 10L,
                            seed = 1L,
                            select = c("objective", "deviance"),
                            prefer_converged = TRUE,
                            verbose = FALSE) {
  select <- match.arg(select)
  n_starts <- as.integer(n_starts)[1L]
  seed <- as.integer(seed)[1L]
  if (!is.finite(n_starts) || n_starts < 1L) {
    stop("`n_starts` must be a positive integer.", call. = FALSE)
  }
  if (!inherits(control, "tt_control")) {
    control <- do.call(tt_control, as.list(control))
  }
  X <- as.matrix(X)
  d <- ncol(X)

  rows <- vector("list", n_starts)
  fits <- vector("list", n_starts)
  best_idx <- NA_integer_
  best_crit <- Inf
  best_conv <- FALSE

  for (s in seq_len(n_starts)) {
    s_seed <- as.integer(seed + s - 1L)
    init <- tt_initialize(
      X, rank = rank, k = k, seed = s_seed,
      sd = control$init_sd %||% 0.15
    )
    ctrl_s <- control
    ctrl_s$seed <- s_seed
    if (isTRUE(verbose)) {
      message(sprintf("ttps_multistart: start %d/%d (seed=%d)", s, n_starts, s_seed))
    }
    fit <- tryCatch(
      ttps(
        y, X,
        family = family,
        rank = rank,
        k = k,
        degree = degree,
        penalty_order = penalty_order,
        lambda = lambda,
        optimizer = optimizer,
        backend = backend,
        init = init,
        control = ctrl_s,
        monitor = monitor,
        knots = knots,
        offset = offset,
        linear = linear,
        weights = weights,
        cyclic = cyclic,
        period = period
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      rows[[s]] <- data.frame(
        start = s,
        seed = s_seed,
        converged = FALSE,
        objective = Inf,
        deviance = Inf,
        timing = NA_real_,
        lambda = NA_character_,
        lambda_at_boundary = NA,
        lambda_boundary = NA_character_,
        error = conditionMessage(fit),
        stringsAsFactors = FALSE
      )
      next
    }
    fits[[s]] <- fit
    obj <- tryCatch(
      as.numeric(tt_objective(fit, X, y)$value),
      error = function(e) NA_real_
    )
    if (!is.finite(obj)) obj <- Inf
    dev <- as.numeric(fit$deviance)
    if (!is.finite(dev)) dev <- Inf
    crit <- if (identical(select, "objective")) obj else dev
    conv <- isTRUE(fit$converged)
    lam <- as.numeric(fit$lambda)
    bound_lab <- if (!is.null(fit$lambda_boundary) && length(fit$lambda_boundary)) {
      paste(fit$lambda_boundary, collapse = ",")
    } else {
      NA_character_
    }
    rows[[s]] <- data.frame(
      start = s,
      seed = s_seed,
      converged = conv,
      objective = obj,
      deviance = dev,
      timing = as.numeric(fit$timing %||% NA_real_),
      lambda = paste(sprintf("%.6g", lam), collapse = ","),
      lambda_at_boundary = isTRUE(fit$lambda_at_boundary),
      lambda_boundary = bound_lab,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
    # selection: prefer converged if requested
    better <- FALSE
    if (isTRUE(prefer_converged)) {
      if (conv && !best_conv) {
        better <- TRUE
      } else if ((conv == best_conv) && is.finite(crit) && crit < best_crit) {
        better <- TRUE
      }
    } else if (is.finite(crit) && crit < best_crit) {
      better <- TRUE
    }
    if (better) {
      best_idx <- s
      best_crit <- crit
      best_conv <- conv
    }
  }

  starts <- do.call(rbind, rows)
  rownames(starts) <- NULL

  # Per-margin boundary fractions among successful finite starts
  frac_boundary <- .ttps_multistart_frac_boundary(fits, d)

  if (is.na(best_idx) || is.null(fits[[best_idx]])) {
    stop(
      "ttps_multistart: every start failed. Inspect `$starts` after catching ",
      "this error is not possible; try a simpler model or check data.",
      call. = FALSE
    )
  }

  best <- fits[[best_idx]]
  best$multistart <- list(
    n_starts = n_starts,
    selected_start = best_idx,
    select = select,
    seed = seed
  )

  structure(
    list(
      best = best,
      starts = starts,
      n_starts = n_starts,
      select = select,
      seed = seed,
      selected_start = best_idx,
      frac_boundary = frac_boundary,
      call = match.call()
    ),
    class = "ttps_multistart"
  )
}

#' @keywords internal
#' @noRd
.ttps_multistart_frac_boundary <- function(fits, d) {
  ok <- vapply(fits, function(f) {
    !is.null(f) && !inherits(f, "error") &&
      identical(f$lambda_method %||% "", "cGCV") &&
      !is.null(f$lambda_boundary) && length(f$lambda_boundary) == d
  }, logical(1))
  if (!any(ok)) {
    out <- rep(NA_real_, d)
    names(out) <- paste0("margin", seq_len(d))
    return(out)
  }
  mat <- do.call(rbind, lapply(fits[ok], function(f) {
    as.character(f$lambda_boundary)
  }))
  frac <- vapply(seq_len(d), function(j) {
    mean(mat[, j] %in% c("lower", "upper", "both"))
  }, numeric(1))
  names(frac) <- paste0("margin", seq_len(d))
  frac
}

#' @export
print.ttps_multistart <- function(x, ...) {
  cat("TTPsplines multi-start summary\n")
  cat(sprintf("  Starts: %d | select: %s | base seed: %d\n",
              x$n_starts, x$select, x$seed))
  cat(sprintf("  Selected start: %d\n", x$selected_start))
  st <- x$starts
  n_ok <- sum(is.finite(st$objective) & is.finite(st$deviance) &
                (is.na(st$error) | !nzchar(st$error)))
  cat(sprintf("  Finite fits: %d / %d | converged: %d\n",
              n_ok, nrow(st), sum(st$converged, na.rm = TRUE)))
  objs <- st$objective[is.finite(st$objective)]
  if (length(objs)) {
    cat(sprintf("  Objective range: [%.6g, %.6g]\n", min(objs), max(objs)))
  }
  if (any(st$lambda_at_boundary, na.rm = TRUE)) {
    cat(sprintf("  Starts at lambda boundary: %d / %d\n",
                sum(st$lambda_at_boundary, na.rm = TRUE),
                sum(!is.na(st$lambda_at_boundary))))
  }
  fb <- x$frac_boundary
  if (!is.null(fb) && any(is.finite(fb))) {
    cat(sprintf("  Frac. boundary by margin: %s\n",
                paste(sprintf("%s=%.2f", names(fb), fb), collapse = ", ")))
    if (any(fb >= 0.8, na.rm = TRUE)) {
      cat("  Note: stable boundary hits across starts → interpret as margin/null/r,\n")
      cat("        not a single bad ALS initialization.\n")
    } else if (any(fb > 0 & fb < 0.8, na.rm = TRUE)) {
      cat("  Note: unstable boundary hits → check initialization / raise n_starts.\n")
    }
  }
  cat("\nBest fit:\n")
  print(x$best)
  invisible(x)
}

#' @export
summary.ttps_multistart <- function(object, ...) {
  structure(object, class = c("summary.ttps_multistart", "ttps_multistart"))
}

#' @export
print.summary.ttps_multistart <- function(x, ...) {
  print.ttps_multistart(x, ...)
  cat("\nStart table:\n")
  print(x$starts[, c("start", "seed", "converged", "objective", "deviance",
                     "lambda_at_boundary", "lambda_boundary")],
        row.names = FALSE)
  invisible(x)
}
