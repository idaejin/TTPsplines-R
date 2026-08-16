#' Leave-one-margin-out tests for dropping TT covariates
#'
#' For each candidate margin \(j\), compare the full TT fit on all columns of
#' `X` to the reduced fit that **omits** column \(j\). Large improvement from
#' keeping \(j\) argues against dropping it.
#'
#' **Null hypothesis (per margin):** margin \(j\) does not improve fit given
#' the others (candidate to eliminate). Soft rule: `p_value > alpha` marks
#' `drop_candidate = TRUE`.
#'
#' Two approximate procedures:
#'
#' - `"nested"` (default): \(\Delta_j = D_{-j}-D_{\mathrm{full}}\) with an F-like
#'   reference using joint linearized EDF when available. *Approximate* for TT.
#' - `"permute"`: permute column \(j\) under the full design `B` times; compare
#'   the leave-one-out \(\Delta_j\) to the null importance
#'   \(D_{\mathrm{full}}(x_j^{\pi})-D_{\mathrm{full}}\).
#'
#' Prefer [tt_margin_activity_path()] for formal subset choice; use this as a
#' complementary drop diagnostic.
#'
#' @param y,X Response and covariate matrix (as in [ttps()]).
#' @param margins Integer indices or character names to test (default: all).
#' @param rank,k,degree,penalty_order,lambda,family,optimizer,backend,control,offset,weights
#'   Passed to [ttps()] (univariate `smooth.spline` when a reduced model has
#'   one column).
#' @param method `"nested"` or `"permute"`.
#' @param B Permutations when `method = "permute"` (default `99`).
#' @param alpha Soft drop threshold (default `0.05`).
#' @param seed RNG seed for permutations.
#' @param verbose Progress messages.
#'
#' @return Class `"tt_margin_drop_test"` with `results`, `keep`,
#'   `drop_candidates`, and `fit_full`.
#'
#' @seealso [tt_margin_activity_path()], [tt_edf()], [ttps()]
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' X <- matrix(runif(n * 4), n, 4)
#' colnames(X) <- paste0("x", 1:4)
#' y <- sin(2 * pi * X[, 1]) + 0.5 * sin(2 * pi * X[, 2]) + rnorm(n, sd = 0.3)
#' tst <- tt_margin_drop_test(
#'   y, X, rank = 2, k = 5, lambda = 1, method = "nested",
#'   control = tt_control(max_sweeps = 4, compute_edf = TRUE, backend = "R"),
#'   verbose = FALSE
#' )
#' tst
#'
#' @export
tt_margin_drop_test <- function(y,
                                X,
                                margins = NULL,
                                rank = 2L,
                                k = 5L,
                                degree = 3L,
                                penalty_order = 2L,
                                lambda = 1,
                                family = stats::gaussian(),
                                optimizer = "auto",
                                backend = "auto",
                                control = tt_control(compute_edf = TRUE),
                                offset = NULL,
                                weights = NULL,
                                method = c("nested", "permute"),
                                B = 99L,
                                alpha = 0.05,
                                seed = NULL,
                                verbose = TRUE) {
  method <- match.arg(method)
  X <- as.matrix(X)
  y <- as.numeric(y)
  stopifnot(length(y) == nrow(X), ncol(X) >= 2L)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  d <- ncol(X)
  n <- nrow(X)

  if (is.null(margins)) {
    margins <- seq_len(d)
  } else if (is.character(margins)) {
    margins <- match(margins, colnames(X))
    if (anyNA(margins)) stop("Unknown margin name(s).", call. = FALSE)
  } else {
    margins <- as.integer(margins)
  }
  stopifnot(all(margins >= 1L), all(margins <= d))

  ctrl <- control
  if (is.null(ctrl$compute_edf)) ctrl$compute_edf <- TRUE

  fit_args <- list(
    rank = rank, k = k, degree = degree, penalty_order = penalty_order,
    lambda = lambda, family = family, optimizer = optimizer,
    backend = backend, control = ctrl, offset = offset, weights = weights
  )

  if (isTRUE(verbose)) message("Fitting full TT model (d = ", d, ") ...")
  fit_full <- do.call(.tt_map_fit, c(list(y = y, X = X), fit_args))
  D_full <- .tt_drop_deviance(fit_full, y)
  edf_full <- .tt_drop_edf_joint(fit_full)
  edf_m <- .tt_drop_edf_margin(fit_full)

  if (!is.null(seed)) set.seed(seed)
  B <- as.integer(B)

  rows <- vector("list", length(margins))
  for (ii in seq_along(margins)) {
    j <- margins[[ii]]
    nm <- colnames(X)[[j]]
    if (isTRUE(verbose)) {
      message(sprintf("  Testing margin %s (%d/%d)", nm, ii, length(margins)))
    }
    cols_red <- setdiff(seq_len(d), j)
    fit_red <- do.call(
      .tt_map_fit,
      c(list(y = y, X = X[, cols_red, drop = FALSE]), fit_args)
    )
    D_red <- .tt_drop_deviance(fit_red, y)
    edf_red <- .tt_drop_edf_joint(fit_red)
    delta <- as.numeric(D_red - D_full)
    edf_j <- if (!is.null(edf_m) && length(edf_m) >= j) {
      as.numeric(edf_m[[j]])
    } else {
      NA_real_
    }

    if (identical(method, "nested")) {
      if (!is.finite(delta) || delta <= 0) {
        # Reduced model no worse (or better) => no evidence to keep j
        Fstat <- 0
        pval <- 1
      } else {
        nu <- edf_full - edf_red
        if (!is.finite(nu) || nu <= 0) {
          nu <- if (is.finite(edf_j) && edf_j > 0) edf_j else 1
        }
        nu <- max(as.numeric(nu), 1e-8)
        den <- D_full / max(n - edf_full, 1e-8)
        if (!is.finite(den) || den <= 0) den <- D_full / max(n - 1, 1)
        Fstat <- (delta / nu) / den
        pval <- if (is.finite(Fstat) && Fstat >= 0) {
          stats::pf(Fstat, df1 = nu, df2 = max(n - edf_full, 1e-8),
                    lower.tail = FALSE)
        } else {
          1
        }
      }
      stat <- Fstat
    } else {
      # Compare leave-one-out Delta to null importance of permuting x_j
      if (!is.finite(delta) || delta <= 0) {
        stat <- 0
        pval <- 1
      } else {
        imp_obs <- delta
        imp_perm <- numeric(B)
        for (b in seq_len(B)) {
          Xp <- X
          Xp[, j] <- sample(X[, j])
          fit_p <- do.call(.tt_map_fit, c(list(y = y, X = Xp), fit_args))
          imp_perm[[b]] <- .tt_drop_deviance(fit_p, y) - D_full
        }
        pval <- (1 + sum(imp_perm >= imp_obs - 1e-12)) / (B + 1)
        stat <- imp_obs
      }
    }

    rows[[ii]] <- data.frame(
      margin = nm,
      margin_index = j,
      delta_deviance = delta,
      edf_full = edf_full,
      edf_reduced = edf_red,
      edf_margin = edf_j,
      stat = stat,
      p_value = pval,
      drop_candidate = is.finite(pval) && pval > alpha,
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  drop_idx <- results$margin_index[which(results$drop_candidate)]
  keep_idx <- setdiff(seq_len(d), drop_idx)

  out <- list(
    results = results,
    method = method,
    alpha = alpha,
    B = if (identical(method, "permute")) B else NA_integer_,
    keep = keep_idx,
    keep_names = colnames(X)[keep_idx],
    drop_candidates = drop_idx,
    drop_candidate_names = colnames(X)[drop_idx],
    fit_full = fit_full,
    call = match.call(),
    variable_names = colnames(X)
  )
  class(out) <- "tt_margin_drop_test"
  out
}

#' @export
print.tt_margin_drop_test <- function(x, digits = 4, ...) {
  cat("TT margin drop test (", x$method, ")\n", sep = "")
  cat("  H0: margin j adds no improvement ",
      "(drop candidate if p > alpha=", x$alpha, ")\n", sep = "")
  print(x$results[, c("margin", "delta_deviance", "edf_margin", "stat",
                      "p_value", "drop_candidate")], digits = digits)
  cat("  Keep: ", paste(x$keep_names, collapse = ", "), "\n", sep = "")
  cat("  Drop candidates: ",
      if (length(x$drop_candidate_names)) {
        paste(x$drop_candidate_names, collapse = ", ")
      } else {
        "(none)"
      },
      "\n", sep = "")
  invisible(x)
}

.tt_drop_deviance <- function(fit, y) {
  if (inherits(fit, "ttpspline") || inherits(fit, "ttps")) {
    if (is.finite(fit$deviance)) return(as.numeric(fit$deviance))
    return(sum((as.numeric(y) - as.numeric(fit$fitted.values))^2))
  }
  if (is.list(fit) && identical(fit$type, "intercept")) {
    return(sum((as.numeric(y) - fit$mu)^2))
  }
  if (is.list(fit) && identical(fit$type, "univariate") && !is.null(fit$ss)) {
    mu <- as.numeric(stats::predict(fit$ss, x = fit$ss$x)$y)
    return(sum((as.numeric(y) - mu)^2))
  }
  stop("Unknown fit type in .tt_drop_deviance().", call. = FALSE)
}

.tt_drop_edf_joint <- function(fit) {
  if (inherits(fit, "ttpspline") || inherits(fit, "ttps")) {
    if (is.finite(fit$edf)) return(as.numeric(fit$edf))
    ed <- tryCatch(tt_edf(fit)$joint, error = function(e) NA_real_)
    return(as.numeric(ed))
  }
  if (is.list(fit) && identical(fit$type, "univariate") && !is.null(fit$ss)) {
    return(as.numeric(fit$ss$df))
  }
  if (is.list(fit) && identical(fit$type, "intercept")) return(1)
  NA_real_
}

.tt_drop_edf_margin <- function(fit) {
  if (inherits(fit, "ttpspline") || inherits(fit, "ttps")) {
    if (!is.null(fit$edf_margin)) return(fit$edf_margin)
    ed <- tryCatch(tt_edf(fit)$margin, error = function(e) NULL)
    return(ed)
  }
  NULL
}
