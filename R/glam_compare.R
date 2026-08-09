#' Compare Currie–Durbán–Eilers Gaussian GLAM vs TT-P-splines.
#'
#' Fits all methods on the **same training grid**, selects TT rank with
#' [tt_rank_select()] using **training data only**, then evaluates
#' \(\mathrm{RMSE}\) vs \(f_{\mathrm{true}}\) on an **independent test sample**.
#'
#' ## Timing (honest computational accounting)
#'
#' | Column | Meaning |
#' |--------|---------|
#' | `fit_time` | Cost of the **final** model fit (known rank / GLAM). |
#' | `rank_selection_time` | Cost of [tt_rank_select()] CV (0 for GLAM / oracle). |
#' | `total_procedure_time` | Practical cost = `fit_time + rank_selection_time`. |
#' | `oracle_probe_time` | Extra cost of probing all ranks for `TT-oracle` (**not**
#'   included in practical totals). |
#'
#' A referee can therefore see both “TT is fast once \(r\) is known” and
#' “selecting \(r\) by \(K\)-fold CV costs \(K\times|\mathcal R|\) fits.”
#'
#' @param d Dimension.
#' @param n_grid Training grid size per margin (scalar or length-`d`).
#' @param k Basis size per margin.
#' @param ranks Candidate uniform TT ranks.
#' @param lambda Fixed smoothing (preferred) or `"cGCV"` for TT
#'   (GLAM then uses fixed \(\lambda=1\) with a note — no GLAM-cGCV yet).
#' @param sigma Gaussian noise sd on the training grid.
#' @param max_sweeps ALS sweeps.
#' @param max_glam_npar Skip GLAM when `k^d` exceeds this.
#' @param degree B-spline degree.
#' @param seed RNG seed (data, test draw, CV folds).
#' @param backend TT backend.
#' @param rank_selection `"cv"` (default) or `"none"` (rank grid only; no
#'   minCV/1SE rows).
#' @param rank_rule `"both"` (emit minCV and 1SE), `"1se"`, or `"min"`.
#' @param folds,select_folds Number of CV folds (`folds` preferred;
#'   `select_folds` kept for backward compatibility).
#' @param include_oracle Include simulation-only `TT-oracle` row.
#' @param detail_ranks Append per-rank `TT-r*` rows (for RMSE-vs-rank plots).
#' @param n_test Independent test sample size (default `min(2000, 4 * n_train)`).
#' @param penalty_order Difference penalty order.
#'
#' @return A `data.frame` of class `"glam_tt_compare"`. Attributes:
#'   `selection`, `oracle_probe_time`, `test` (list with `X`, `truth`).
#'
#' @seealso [summarize_glam_tt_compare()], [plot.glam_tt_compare()],
#'   [tt_rank_select()].
#' @export
#' @examples
#' cmp <- compare_glam_tt_gaussian(
#'   d = 3, n_grid = c(8, 7, 6), k = 5, ranks = 1:3,
#'   lambda = 1, folds = 3, max_sweeps = 4, seed = 1
#' )
#' cmp[, c("method", "rank", "fit_time", "rank_selection_time",
#'         "total_procedure_time", "rmse_truth")]
compare_glam_tt_gaussian <- function(d = 3L,
                                     n_grid = NULL,
                                     k = NULL,
                                     ranks = 1:3,
                                     lambda = 1,
                                     sigma = 0.25,
                                     max_sweeps = 10L,
                                     max_glam_npar = 5000L,
                                     degree = 3L,
                                     seed = 51L,
                                     backend = "R",
                                     rank_selection = c("cv", "none"),
                                     rank_rule = c("both", "1se", "min"),
                                     folds = NULL,
                                     select_folds = 4L,
                                     include_oracle = TRUE,
                                     detail_ranks = FALSE,
                                     n_test = NULL,
                                     penalty_order = 2L) {
  rank_selection <- match.arg(rank_selection)
  rank_rule <- match.arg(rank_rule)
  d <- as.integer(d)
  if (d < 2L) stop("`d` must be >= 2.", call. = FALSE)
  degree <- as.integer(degree)[1L]
  penalty_order <- as.integer(penalty_order)[1L]

  if (is.null(n_grid)) {
    n_grid <- if (d == 3L) {
      c(14L, 12L, 10L)
    } else if (d == 5L) {
      rep(6L, 5L)
    } else if (d == 7L) {
      rep(3L, 7L)
    } else {
      rep(max(3L, 10L - d), d)
    }
  }
  n_grid <- rep(as.integer(n_grid), length.out = d)
  if (is.null(k)) {
    k <- if (d <= 3L) 6L else if (d <= 5L) 5L else max(degree + 1L, 3L)
  }
  k <- as.integer(k)[1L]
  if (k <= degree) stop("`k` must be > degree (", degree, ").", call. = FALSE)

  ranks <- sort(unique(as.integer(ranks)))
  if (!length(ranks) || any(ranks < 1L)) {
    stop("`ranks` must be positive integers.", call. = FALSE)
  }

  folds_k <- as.integer(if (!is.null(folds)) folds else select_folds)[1L]
  if (identical(rank_selection, "cv") && (folds_k < 2L)) {
    stop("`folds` must be >= 2 when rank_selection = \"cv\".", call. = FALSE)
  }

  lambda_is_cgcv <- is.character(lambda) && identical(lambda, "cGCV")
  if (lambda_is_cgcv) {
    lambda_tt <- "cGCV"
    lambda_glam <- rep(1, d)
    lambda_note_glam <- "GLAM fixed lambda=1 (no GLAM-cGCV); TT used cGCV"
  } else {
    lambda_tt <- rep(as.numeric(lambda), length.out = d)
    lambda_glam <- lambda_tt
    lambda_note_glam <- "ok"
  }

  max_glam_npar <- as.integer(max_glam_npar)[1L]
  set.seed(as.integer(seed))

  # ---- training grid -------------------------------------------------------
  axes <- lapply(n_grid, function(nk) seq(0, 1, length.out = nk))
  names(axes) <- paste0("x", seq_len(d))
  idx <- expand.grid(lapply(n_grid, seq_len), KEEP.OUT.ATTRS = FALSE)
  X_train <- do.call(cbind, lapply(seq_len(d), function(j) axes[[j]][idx[[j]]]))
  colnames(X_train) <- names(axes)
  truth_train <- .grid_truth_surface(X_train)
  Y_vec <- truth_train + stats::rnorm(length(truth_train), 0, sigma)
  Y <- array(Y_vec, dim = n_grid)
  n_train <- length(Y_vec)
  n_grid_margin <- as.integer(round(exp(mean(log(as.numeric(n_grid))))))
  npar_dense <- as.integer(k)^d
  glam_ok <- npar_dense <= max_glam_npar

  # ---- independent test sample (evaluation only; never for selection) ------
  if (is.null(n_test)) n_test <- as.integer(min(2000L, max(200L, 4L * n_train)))
  n_test <- as.integer(n_test)[1L]
  X_test <- matrix(stats::runif(n_test * d), n_test, d)
  colnames(X_test) <- colnames(X_train)
  truth_test <- .grid_truth_surface(X_test)
  # noisy test responses (optional rmse_y)
  y_test <- truth_test + stats::rnorm(n_test, 0, sigma)

  bb <- glam_grid_bases(axes, k = k, degree = degree)
  rmse <- function(a, b) {
    if (is.null(a) || is.null(b)) return(NA_real_)
    sqrt(mean((as.numeric(a) - as.numeric(b))^2))
  }

  ctrl <- tt_control(
    max_sweeps = as.integer(max_sweeps),
    backend = backend,
    compute_edf = FALSE,
    seed = as.integer(seed) + 7L
  )

  # ---- GLAM ----------------------------------------------------------------
  if (isTRUE(glam_ok)) {
    t0 <- proc.time()[["elapsed"]]
    glam <- glam_fit_gaussian(Y, bb$B, lambda = lambda_glam,
                              penalty_order = penalty_order)
    fit_time_glam <- proc.time()[["elapsed"]] - t0
    mu_glam_test <- .glam_predict_scattered(
      glam$Theta, glam$intercept, bb$knots, degree, X_test
    )
    glam_npar <- as.integer(glam$npar)
    glam_status <- "ok"
    glam_reason <- if (lambda_is_cgcv) lambda_note_glam else NA_character_
    glam_method <- "GLAM"
    glam_opt <- glam$method
  } else {
    fit_time_glam <- NA_real_
    mu_glam_test <- NULL
    glam_npar <- npar_dense
    glam_status <- "not_run"
    glam_reason <- sprintf(
      "dense tensor exceeds benchmark threshold (k^d=%d > max_glam_npar=%d)",
      npar_dense, max_glam_npar
    )
    glam_method <- "GLAM-infeasible"
    glam_opt <- "GLAM-skipped"
  }

  .row <- function(method, role, rank, npar, intrinsic_dim, compression,
                   mu_test, fit_time, rank_selection_time = 0,
                   oracle_probe_time = NA_real_,
                   status = "ok", reason = NA_character_,
                   backend = NA_character_, optimizer = NA_character_,
                   note = "ok") {
    fit_time <- as.numeric(fit_time)
    rank_selection_time <- as.numeric(rank_selection_time)
    total <- if (is.finite(fit_time)) {
      fit_time + (if (is.finite(rank_selection_time)) rank_selection_time else 0)
    } else {
      NA_real_
    }
    # Practical total excludes oracle probe by construction
    data.frame(
      d = d,
      n_grid_margin = n_grid_margin,
      n_cells = n_train,
      n_test = n_test,
      k = k,
      degree = degree,
      method = method,
      role = role,
      rank = if (is.null(rank) || !length(rank) || all(is.na(rank))) {
        NA_integer_
      } else {
        as.integer(rank)[1L]
      },
      npar = as.integer(npar),
      intrinsic_dim = if (is.null(intrinsic_dim) || !length(intrinsic_dim)) {
        NA_integer_
      } else {
        as.integer(intrinsic_dim)[1L]
      },
      npar_dense = npar_dense,
      compression = as.numeric(compression),
      rmse_truth = rmse(mu_test, truth_test),
      rmse_y = rmse(mu_test, y_test),
      rmse_vs_glam = if (!is.null(mu_glam_test)) {
        rmse(mu_test, mu_glam_test)
      } else {
        NA_real_
      },
      fit_time = fit_time,
      rank_selection_time = rank_selection_time,
      total_procedure_time = total,
      oracle_probe_time = as.numeric(oracle_probe_time),
      # backward-compatible aliases
      time_s = fit_time,
      select_time_s = rank_selection_time,
      time_vs_glam = if (isTRUE(glam_ok) && is.finite(fit_time) &&
                          is.finite(fit_time_glam)) {
        fit_time / max(fit_time_glam, 1e-12)
      } else {
        NA_real_
      },
      status = status,
      reason = if (is.null(reason) || length(reason) == 0L) {
        NA_character_
      } else {
        as.character(reason)[1L]
      },
      backend = backend,
      optimizer = optimizer,
      note = note,
      stringsAsFactors = FALSE
    )
  }

  rows <- list(.row(
    method = glam_method,
    role = "baseline",
    rank = NA_integer_,
    npar = glam_npar,
    intrinsic_dim = glam_npar,
    compression = 1,
    mu_test = mu_glam_test,
    fit_time = fit_time_glam,
    rank_selection_time = 0,
    status = glam_status,
    reason = glam_reason,
    backend = "array",
    optimizer = glam_opt,
    note = if (identical(glam_status, "ok")) "ok" else glam_reason
  ))

  # ---- TT fits on training -------------------------------------------------
  fit_tt <- function(r) {
    t1 <- proc.time()[["elapsed"]]
    fit <- ttpspline(
      Y_vec, X_train,
      family = stats::gaussian(),
      rank = r,
      k = k,
      degree = degree,
      penalty_order = penalty_order,
      lambda = lambda_tt,
      knots = bb$knots,
      control = ctrl
    )
    list(fit = fit, fit_time = proc.time()[["elapsed"]] - t1)
  }

  predict_tt <- function(fit) {
    predict(fit, newdata = X_test, type = "response")
  }

  # ---- CV rank selection (training only; no truth) -------------------------
  sel <- NULL
  select_time <- 0
  if (identical(rank_selection, "cv")) {
    t_sel0 <- proc.time()[["elapsed"]]
    sel <- tt_rank_select(
      Y_vec, X_train,
      ranks = ranks,
      family = stats::gaussian(),
      k = k,
      degree = degree,
      penalty_order = penalty_order,
      lambda = lambda_tt,
      folds = folds_k,
      rule = "1se",
      seed = as.integer(seed) + 11L,
      knots = bb$knots,
      control = ctrl
    )
    select_time <- proc.time()[["elapsed"]] - t_sel0
  }

  # Cache training fits we need
  need <- ranks
  if (!is.null(sel)) need <- unique(c(need, sel$rank_min, sel$rank_1se))
  cache <- vector("list", length(need))
  names(cache) <- as.character(need)
  for (r in need) cache[[as.character(r)]] <- fit_tt(r)

  # Oracle probe: RMSE on TEST truth after fitting (truth never in selection)
  oracle_probe_time <- NA_real_
  r_oracle <- NA_integer_
  if (isTRUE(include_oracle)) {
    oracle_rmse <- vapply(ranks, function(r) {
      rmse(predict_tt(cache[[as.character(r)]]$fit), truth_test)
    }, numeric(1))
    r_oracle <- ranks[which.min(oracle_rmse)]
    # Cost of probing all candidate ranks (simulation accounting only)
    oracle_probe_time <- sum(vapply(ranks, function(r) {
      cache[[as.character(r)]]$fit_time
    }, numeric(1)))
  }

  add_tt <- function(method, role, r, rank_selection_time = 0,
                     oracle_probe_time = NA_real_, note = "ok") {
    z <- cache[[as.character(r)]]
    fit <- z$fit
    rows[[length(rows) + 1L]] <<- .row(
      method = method,
      role = role,
      rank = r,
      npar = fit$npar_tt,
      intrinsic_dim = fit$npar_tt_intrinsic %||% NA_integer_,
      compression = fit$compression_ratio,
      mu_test = predict_tt(fit),
      fit_time = z$fit_time,
      rank_selection_time = rank_selection_time,
      oracle_probe_time = oracle_probe_time,
      status = "ok",
      backend = fit$backend,
      optimizer = fit$optimizer_used,
      note = note
    )
  }

  if (isTRUE(include_oracle) && is.finite(r_oracle)) {
    add_tt(
      "TT-oracle", "oracle", r_oracle,
      rank_selection_time = 0,
      oracle_probe_time = oracle_probe_time,
      note = "simulation-only oracle reference; uses f_true on test after fitting"
    )
  }

  if (identical(rank_selection, "cv") && !is.null(sel)) {
    if (rank_rule %in% c("both", "min")) {
      add_tt("TT-minCV", "cv-min", sel$rank_min,
             rank_selection_time = select_time,
             note = "CV minimum; training only; no truth")
    }
    if (rank_rule %in% c("both", "1se")) {
      add_tt("TT-1SE", "1se", sel$rank_1se,
             rank_selection_time = select_time,
             note = "CV 1-SE; training only; no truth")
    }
  }

  if (isTRUE(detail_ranks)) {
    for (r in ranks) {
      add_tt(sprintf("TT-r%d", r), "grid", r, note = "rank grid (detail)")
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # Deltas vs oracle (prediction cost of data-driven selection)
  if (any(out$method == "TT-oracle")) {
    ora_rmse <- out$rmse_truth[out$method == "TT-oracle"][1L]
    out$delta_rmse_oracle <- out$rmse_truth - ora_rmse
    out$delta_rel_oracle <- ifelse(
      is.finite(ora_rmse) & abs(ora_rmse) > 1e-12,
      (out$rmse_truth - ora_rmse) / ora_rmse,
      NA_real_
    )
    out$delta_rmse_oracle[out$method %in% c("GLAM", "GLAM-infeasible", "TT-oracle")] <-
      NA_real_
    out$delta_rel_oracle[out$method %in% c("GLAM", "GLAM-infeasible", "TT-oracle")] <-
      NA_real_
  } else {
    out$delta_rmse_oracle <- NA_real_
    out$delta_rel_oracle <- NA_real_
  }

  class(out) <- c("glam_tt_compare", "data.frame")
  attr(out, "selection") <- sel
  attr(out, "oracle_probe_time") <- oracle_probe_time
  attr(out, "test") <- list(X = X_test, truth = truth_test, y = y_test)
  attr(out, "train") <- list(X = X_train, y = Y_vec, truth = truth_train)
  attr(out, "call") <- match.call()
  out
}

#' @export
print.glam_tt_compare <- function(x, digits = 4, ...) {
  cat("GLAM vs TT comparison (Gaussian)\n")
  cat("Timing: fit_time = final fit; rank_selection_time = CV cost;\n")
  cat("        total_procedure_time = fit_time + rank_selection_time.\n")
  cat("        oracle_probe_time is simulation-only (not in totals).\n\n")
  cols <- c("method", "rank", "npar", "intrinsic_dim", "compression",
            "rmse_truth", "fit_time", "rank_selection_time",
            "total_procedure_time", "status")
  cols <- cols[cols %in% names(x)]
  print(as.data.frame(x)[, cols, drop = FALSE], digits = digits, row.names = FALSE)
  invisible(x)
}

#' Compact summary: Oracle | minCV | 1-SE | GLAM
#'
#' @param tab Output of [compare_glam_tt_gaussian()] / [compare_glam_tt_scale()].
#' @return A `data.frame` summary with timing columns separated.
#' @export
summarize_glam_tt_compare <- function(tab) {
  stopifnot(is.data.frame(tab), "method" %in% names(tab))
  ft <- function(row, name) {
    if (is.null(row)) return(NA_real_)
    if (name %in% names(row)) return(as.numeric(row[[name]]))
    if (identical(name, "fit_time") && "time_s" %in% names(row)) {
      return(as.numeric(row$time_s))
    }
    if (identical(name, "rank_selection_time") && "select_time_s" %in% names(row)) {
      return(as.numeric(row$select_time_s))
    }
    NA_real_
  }
  keys <- paste(tab$d, tab$n_grid_margin, tab$k, sep = ":")
  blocks <- split(tab, keys)
  rows <- lapply(blocks, function(b) {
    pick <- function(meth) {
      w <- b[b$method == meth, , drop = FALSE]
      if (!nrow(w)) return(NULL)
      w[1, ]
    }
    glam <- pick("GLAM")
    if (is.null(glam)) glam <- pick("GLAM-infeasible")
    ora <- pick("TT-oracle")
    mcv <- pick("TT-minCV")
    se1 <- pick("TT-1SE")
    data.frame(
      d = glam$d,
      n_grid = glam$n_grid_margin,
      n_train = glam$n_cells,
      n_test = if ("n_test" %in% names(glam)) glam$n_test else NA_integer_,
      k = glam$k,
      glam = glam$method,
      glam_status = if ("status" %in% names(glam)) glam$status else "ok",
      glam_npar = glam$npar,
      glam_rmse = glam$rmse_truth,
      glam_fit_time = ft(glam, "fit_time"),
      oracle_rank = if (is.null(ora)) NA_integer_ else ora$rank,
      oracle_npar = if (is.null(ora)) NA_integer_ else ora$npar,
      oracle_CR = if (is.null(ora)) NA_real_ else round(ora$compression, 1),
      oracle_rmse = if (is.null(ora)) NA_real_ else ora$rmse_truth,
      oracle_fit_time = ft(ora, "fit_time"),
      oracle_probe_time = ft(ora, "oracle_probe_time"),
      minCV_rank = if (is.null(mcv)) NA_integer_ else mcv$rank,
      minCV_npar = if (is.null(mcv)) NA_integer_ else mcv$npar,
      minCV_CR = if (is.null(mcv)) NA_real_ else round(mcv$compression, 1),
      minCV_rmse = if (is.null(mcv)) NA_real_ else mcv$rmse_truth,
      minCV_fit_time = ft(mcv, "fit_time"),
      minCV_select_time = ft(mcv, "rank_selection_time"),
      minCV_total_time = ft(mcv, "total_procedure_time"),
      minCV_delta = if (is.null(mcv)) NA_real_ else mcv$delta_rmse_oracle,
      se1_rank = if (is.null(se1)) NA_integer_ else se1$rank,
      se1_npar = if (is.null(se1)) NA_integer_ else se1$npar,
      se1_CR = if (is.null(se1)) NA_real_ else round(se1$compression, 1),
      se1_rmse = if (is.null(se1)) NA_real_ else se1$rmse_truth,
      se1_fit_time = ft(se1, "fit_time"),
      se1_select_time = ft(se1, "rank_selection_time"),
      se1_total_time = ft(se1, "total_procedure_time"),
      se1_delta = if (is.null(se1)) NA_real_ else se1$delta_rmse_oracle,
      se1_delta_rel = if (is.null(se1)) NA_real_ else se1$delta_rel_oracle,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Analyze GLAM-vs-TT with data-driven ranks.
#' @inheritParams summarize_glam_tt_compare
#' @export
analyze_glam_tt_rank_select <- function(tab) {
  stopifnot(is.data.frame(tab))
  if (!any(tab$method %in% c("TT-minCV", "TT-1SE"))) {
    stop("tab must include TT-minCV / TT-1SE from compare_glam_tt_gaussian().",
         call. = FALSE)
  }
  tbl <- summarize_glam_tt_compare(tab)
  rel <- function(a, b) ifelse(is.finite(a) & is.finite(b) & abs(b) > 1e-12, a / b, NA_real_)
  pct <- function(a, b) ifelse(is.finite(a) & is.finite(b) & abs(b) > 1e-12, 100 * (a / b - 1), NA_real_)
  tbl$minCV_vs_oracle_rmse <- rel(tbl$minCV_rmse, tbl$oracle_rmse)
  tbl$se1_vs_oracle_rmse <- rel(tbl$se1_rmse, tbl$oracle_rmse)
  tbl$se1_rmse_penalty_pct <- pct(tbl$se1_rmse, tbl$oracle_rmse)
  tbl$minCV_rmse_penalty_pct <- pct(tbl$minCV_rmse, tbl$oracle_rmse)
  tbl$se1_vs_glam_rmse <- rel(tbl$se1_rmse, tbl$glam_rmse)
  tbl$glam_feasible <- grepl("^GLAM$", tbl$glam)
  # Selection overhead vs final fit
  tbl$se1_select_over_fit <- rel(tbl$se1_select_time, tbl$se1_fit_time)

  msgs <- c(
    sprintf("Settings: %d (%d GLAM-feasible).", nrow(tbl), sum(tbl$glam_feasible)),
    sprintf(
      "Mean 1-SE RMSE penalty vs oracle: %+.1f%%.",
      mean(tbl$se1_rmse_penalty_pct, na.rm = TRUE)
    ),
    sprintf(
      "Mean 1-SE select_time / fit_time: %.1fx (CV dominates if >>1).",
      mean(tbl$se1_select_over_fit, na.rm = TRUE)
    ),
    "Oracle is simulation-only; practical totals exclude oracle_probe_time."
  )
  structure(
    list(table = tbl, messages = msgs, call = match.call()),
    class = "glam_tt_rank_select_analysis"
  )
}

#' @export
print.glam_tt_rank_select_analysis <- function(x, digits = 4, ...) {
  cat("GLAM-vs-TT with CV rank selection\n")
  cat("Central: Oracle | minCV | 1-SE | GLAM  (timings separated)\n\n")
  show <- x$table[, intersect(c(
    "d", "k",
    "oracle_rank", "oracle_rmse", "oracle_fit_time", "oracle_probe_time",
    "minCV_rank", "minCV_rmse", "minCV_fit_time", "minCV_select_time", "minCV_total_time",
    "se1_rank", "se1_rmse", "se1_fit_time", "se1_select_time", "se1_total_time",
    "glam", "glam_rmse", "glam_fit_time",
    "se1_delta", "se1_delta_rel"
  ), names(x$table)), drop = FALSE]
  print(show, digits = digits, row.names = FALSE)
  cat("\n")
  for (m in x$messages) cat(" - ", m, "\n", sep = "")
  invisible(x)
}

#' Plot GLAM-vs-TT comparison.
#'
#' @param x A `"glam_tt_compare"` object (use `detail_ranks = TRUE` for
#'   `type = "rmse_rank"`).
#' @param type `"rmse_rank"` or `"tradeoff"`.
#' @param ... Passed to [graphics::plot()].
#' @export
plot.glam_tt_compare <- function(x, type = c("rmse_rank", "tradeoff"), ...) {
  type <- match.arg(type)
  if (identical(type, "rmse_rank")) {
    grid <- x[x$role == "grid" | grepl("^TT-r", x$method), , drop = FALSE]
    if (!nrow(grid)) {
      stop("No rank-grid rows; re-run with detail_ranks = TRUE.", call. = FALSE)
    }
    plot(grid$rank, grid$rmse_truth, type = "b", pch = 19,
         xlab = "TT rank r", ylab = expression(RMSE[test] ~ "vs truth"),
         main = "Test RMSE vs rank", ...)
    mark <- function(meth, col, lty) {
      w <- x[x$method == meth, , drop = FALSE]
      if (nrow(w)) graphics::abline(v = w$rank[1], col = col, lty = lty)
    }
    mark("TT-oracle", "darkgreen", 2)
    mark("TT-minCV", "blue", 3)
    mark("TT-1SE", "darkred", 4)
    glam <- x[x$method == "GLAM", , drop = FALSE]
    if (nrow(glam) && is.finite(glam$rmse_truth[1])) {
      graphics::abline(h = glam$rmse_truth[1], col = "gray40", lty = 1)
    }
    graphics::legend("topright",
                     legend = c("TT grid", "oracle", "minCV", "1SE", "GLAM"),
                     col = c("black", "darkgreen", "blue", "darkred", "gray40"),
                     lty = c(1, 2, 3, 4, 1), pch = c(19, NA, NA, NA, NA),
                     bty = "n", cex = 0.8)
  } else {
    keep <- x$method %in% c("GLAM", "TT-oracle", "TT-minCV", "TT-1SE")
    z <- x[keep, , drop = FALSE]
    z <- z[is.finite(z$rmse_truth) & is.finite(z$compression), , drop = FALSE]
    plot(z$compression, z$rmse_truth, log = "x", pch = 19,
         xlab = "Compression (dense / stored)", ylab = expression(RMSE[test]),
         main = "Accuracy vs compression", ...)
    labs <- z$method
    graphics::text(z$compression, z$rmse_truth, labels = labs, pos = 3, cex = 0.7)
  }
  invisible(x)
}

#' Scale wrapper over [compare_glam_tt_gaussian()].
#' @inheritParams compare_glam_tt_gaussian
#' @param n_grid_sizes,k_values Isotropic grids / basis sizes.
#' @export
compare_glam_tt_scale <- function(d = 7L,
                                  n_grid_sizes = 3L,
                                  k_values = c(4L, 5L),
                                  ranks = 1:2,
                                  lambda = 1,
                                  sigma = 0.25,
                                  max_sweeps = 6L,
                                  max_glam_npar = 5000L,
                                  degree = 3L,
                                  seed = 71L,
                                  backend = "R",
                                  rank_selection = "cv",
                                  rank_rule = "both",
                                  folds = NULL,
                                  select_folds = 4L,
                                  include_oracle = TRUE,
                                  detail_ranks = FALSE,
                                  n_test = NULL) {
  d <- as.integer(d)
  blocks <- list()
  i <- 0L
  for (ng in as.integer(n_grid_sizes)) {
    for (kk in as.integer(k_values)) {
      i <- i + 1L
      blocks[[length(blocks) + 1L]] <- compare_glam_tt_gaussian(
        d = d, n_grid = rep(ng, d), k = kk, ranks = ranks,
        lambda = lambda, sigma = sigma, max_sweeps = max_sweeps,
        max_glam_npar = max_glam_npar, degree = degree,
        seed = as.integer(seed) + i, backend = backend,
        rank_selection = rank_selection, rank_rule = rank_rule,
        folds = folds, select_folds = select_folds,
        include_oracle = include_oracle, detail_ranks = detail_ranks,
        n_test = n_test
      )
    }
  }
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  class(out) <- c("glam_tt_compare", "data.frame")
  out
}
