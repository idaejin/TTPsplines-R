#' Data-driven TT rank selection via K-fold CV (+ optional 1-SE rule).
#'
#' Chooses a uniform TT rank \(r\) (chain `(1,r,...,r,1)`) using out-of-sample
#' predictive loss only — never simulation truth. Fitting stays separate:
#' [ttps()] still requires an explicit `rank`; use [tt_rank_refit()] /
#' [refit()] after selection to fit on all data at the chosen rank.
#'
#' Inner smoothing (`lambda` fixed or `"cGCV"`) is estimated **only** on each
#' fold's training portion (no validation leakage into cGCV).
#'
#' Selected ranks are a **predictive / working structural capacity** conditional
#' on the covariate ordering in `X` — not an estimate of a latent “true TT rank”.
#'
#' @param y,X Response and covariate matrix (as in [ttps()]).
#' @param ranks Integer vector of candidate **uniform** ranks (default `1:5`).
#' @param family,k,degree,penalty_order,lambda,optimizer,backend,control,offset,weights,linear
#'   Passed through to [ttps()] on each fold (and on final refit).
#' @param knots Optional list of knot vectors (as in [ttps()]). Default `NULL`
#'   builds knots once from the **full** `X` (shared across folds and stored for
#'   [tt_rank_refit()]), avoiding B-spline extrapolation warnings when a test
#'   fold has extremes outside a train-only knot span. Set `fold_knots = TRUE`
#'   to restore per-fold knot construction from the training subset only.
#' @param fold_knots If `TRUE`, ignore the full-`X` default and let each fold
#'   rebuild knots from its training rows (`knots` must be `NULL`). Default
#'   `FALSE`.
#' @param folds Number of CV folds (default 5). Ignored when `foldid` is supplied.
#' @param foldid Optional integer vector of length `n` giving the fold of each
#'   observation (cv.glmnet-style). When supplied, folds are taken from
#'   `foldid` (remapped to `1:K` if needed) and shared across all ranks.
#' @param rule `"1se"` (default, parsimonious) or `"min"` (minimum mean CV).
#' @param metric `"auto"` (family-aware) or one of `"mse"`, `"rmse"`,
#'   `"deviance"`, `"poisson_deviance"`, `"binomial_deviance"`, `"logloss"`.
#'   Defaults: Gaussian → MSE; Poisson → mean Poisson deviance;
#'   Bernoulli → mean binomial deviance.
#' @param n_starts Number of random TT initializations per fold×rank (default `1`).
#'   When `>1`, each fold keeps the start with best **training** objective among
#'   converged fits (else best objective) for validation scoring. Recommended
#'   when ALS at low rank is init-sensitive (e.g. Ishigami at `r=2`).
#' @param seed Optional RNG seed for fold assignment (when `foldid` is `NULL`)
#'   and multi-start inits.
#' @param keep_fits If `TRUE`, store fold fits (large); default `FALSE`.
#' @param rank_chain Reserved for future non-uniform rank search; must be `NULL`
#'   in this version.
#' @param ... Currently unused in [tt_rank_select()] (reserved). For [cv.ttps()],
#'   further arguments are passed to [tt_rank_select()].
#'
#' @return An object of class `c("cv.ttps", "tt_rank_selection")` with
#'   glmnet-style aliases `cvm`, `cvsd`, `cvraw`, `foldid`, `rank.min`,
#'   `rank.1se`, plus `variable_order`. With `n_starts>1`, `cv_results` adds
#'   `objective_best`, `objective_median`, `objective_sd`, `start_gap`
#'   (median−best), and `start_convergence_rate`.
#'
#' @section Complexity layers:
#' Rank \(r\) = structural / interaction capacity; \(\lambda\) = directional
#' smoothness; EDF = effective fitted flexibility. These are not interchangeable:
#' \(r \neq \lambda \neq \mathrm{EDF}\). A modestly larger \(r\) can also
#' improve **optimization robustness**.
#'
#' @section Failed folds:
#' Non-finite fold losses are retained in the mean CV (`Inf` ⇒ `cvm = Inf`).
#' The SE uses only finite folds (`sd / sqrt(K_finite)` when `K_finite >= 2`);
#' it is never computed by silently dropping failures from the mean.
#'
#' @seealso [cv.ttps()], [tt_rank_refit()], [refit()], [ttps()],
#'   [tt_truncate_rank()], [tt_rank_profile()] (in-sample rank diagnostic, not CV).
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
                           foldid = NULL,
                           rule = c("1se", "min"),
                           metric = c("auto", "mse", "rmse", "deviance",
                                      "poisson_deviance", "binomial_deviance",
                                      "logloss"),
                           n_starts = 1L,
                           seed = NULL,
                           keep_fits = FALSE,
                           knots = NULL,
                           fold_knots = FALSE,
                           offset = NULL,
                           linear = NULL,
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
  n_starts <- as.integer(n_starts)[1L]
  if (!is.finite(n_starts) || n_starts < 1L) {
    stop("`n_starts` must be a positive integer.", call. = FALSE)
  }

  offset_full <- normalize_offset(offset, n)
  weights_full <- normalize_weights(weights, n)
  linear_full <- normalize_linear(linear, n)

  fold_knots <- isTRUE(fold_knots)
  if (fold_knots && !is.null(knots)) {
    stop("`fold_knots = TRUE` requires `knots = NULL`.", call. = FALSE)
  }
  # Shared full-X knots by default: train/test use the same B-spline domain so
  # fold extremes do not fall outside train-only knot spans (predict warnings).
  knots_source <- if (!is.null(knots)) {
    "user"
  } else if (fold_knots) {
    "fold"
  } else {
    knots <- build_marginal_bases(X, k = k, degree = degree)$knots
    "full_X"
  }

  if (!is.null(foldid)) {
    fold_id <- .tt_normalize_foldid(foldid, n)
    folds <- length(unique(fold_id))
  } else {
    folds <- as.integer(folds)[1L]
    if (folds < 2L || folds > n) {
      stop("`folds` must satisfy 2 <= folds <= n.", call. = FALSE)
    }
    fold_id <- .tt_make_fold_id(n, folds, seed)
  }
  seed0 <- if (is.null(seed)) 1L else as.integer(seed)[1L]

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
    knots = knots,
    linear = linear_full
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
  start_obj <- vector("list", length(ranks))
  start_conv <- vector("list", length(ranks))
  names(start_obj) <- names(start_conv) <- as.character(ranks)
  fold_fits <- if (isTRUE(keep_fits)) {
    vector("list", length(ranks))
  } else {
    NULL
  }
  if (isTRUE(keep_fits)) names(fold_fits) <- as.character(ranks)

  for (i in seq_along(ranks)) {
    r <- ranks[i]
    lambda_list[[i]] <- vector("list", folds)
    start_obj[[i]] <- matrix(NA_real_, folds, n_starts)
    start_conv[[i]] <- matrix(FALSE, folds, n_starts)
    if (isTRUE(keep_fits)) fold_fits[[i]] <- vector("list", folds)
    for (f in seq_len(folds)) {
      test <- fold_id == f
      train <- !test
      off_tr <- offset_full[train]
      off_te <- offset_full[test]
      w_tr <- weights_full[train]
      w_te <- weights_full[test]
      lin_tr <- if (is.null(linear_full)) NULL else linear_full[train, , drop = FALSE]
      lin_te <- if (is.null(linear_full)) NULL else linear_full[test, , drop = FALSE]
      y_tr <- y[train]
      X_tr <- X[train, , drop = FALSE]
      X_te <- X[test, , drop = FALSE]

      best_fit <- NULL
      best_obj <- Inf
      best_conv <- FALSE
      best_lam <- NA_real_
      t_fold <- 0
      for (s in seq_len(n_starts)) {
        init_seed <- as.integer(seed0 + 1000L * r + 100L * f + s)
        init <- tt_initialize(
          X_tr, rank = r, k = k, seed = init_seed,
          sd = ctrl$init_sd %||% 0.15
        )
        ctrl_s <- ctrl
        ctrl_s$seed <- init_seed
        t0 <- proc.time()[["elapsed"]]
        fit <- tryCatch(
          ttps(
            y_tr, X_tr,
            family = fam,
            rank = r,
            k = k,
            degree = degree,
            penalty_order = penalty_order,
            lambda = lambda,
            optimizer = optimizer,
            backend = backend,
            init = init,
            control = ctrl_s,
            knots = knots,
            offset = off_tr,
            linear = lin_tr,
            weights = w_tr
          ),
          error = function(e) e
        )
        t_fold <- t_fold + (proc.time()[["elapsed"]] - t0)
        if (inherits(fit, "error")) {
          start_obj[[i]][f, s] <- Inf
          start_conv[[i]][f, s] <- FALSE
          next
        }
        obj <- tryCatch({
          o <- tt_objective(fit, X_tr, y_tr)
          as.numeric(o$value)
        }, error = function(e) as.numeric(fit$deviance))
        if (!is.finite(obj)) obj <- Inf
        start_obj[[i]][f, s] <- obj
        start_conv[[i]][f, s] <- isTRUE(fit$converged)
        better <- FALSE
        if (isTRUE(fit$converged) && !best_conv) {
          better <- TRUE
        } else if ((isTRUE(fit$converged) == best_conv) && (obj < best_obj)) {
          better <- TRUE
        }
        if (better) {
          best_fit <- fit
          best_obj <- obj
          best_conv <- isTRUE(fit$converged)
          best_lam <- as.numeric(fit$lambda)
        }
      }
      time_mat[i, f] <- t_fold
      if (is.null(best_fit)) {
        loss_mat[i, f] <- Inf
        conv_mat[i, f] <- FALSE
        lambda_list[[i]][[f]] <- NA_real_
        next
      }
      mu <- tryCatch(
        predict(best_fit, newdata = X_te, type = "response",
                offset = off_te, linear = lin_te),
        error = function(e) NULL
      )
      if (is.null(mu) || anyNA(mu) || !all(is.finite(mu))) {
        loss_mat[i, f] <- Inf
        conv_mat[i, f] <- FALSE
        lambda_list[[i]][[f]] <- best_lam
        next
      }
      loss_mat[i, f] <- .tt_cv_loss(y[test], mu, metric = metric, family = fam,
                                   weights = w_te)
      conv_mat[i, f] <- best_conv
      lambda_list[[i]][[f]] <- best_lam
      if (isTRUE(keep_fits)) fold_fits[[i]][[f]] <- best_fit
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
    p = p,
    start_obj = start_obj,
    start_conv = start_conv,
    n_starts = n_starts
  )

  rank_min <- .tt_rank_min_cv(cv_results)
  rank_1se <- .tt_rank_1se(cv_results, rank_min)
  selected <- if (identical(rule, "min")) rank_min else rank_1se

  out <- structure(
    list(
      ranks = ranks,
      metric = metric,
      metric_requested = metric_arg,
      folds = folds,
      n_starts = n_starts,
      fold_id = fold_id,
      foldid = fold_id,
      cv_results = cv_results,
      loss_by_fold = loss_mat,
      time_by_fold = time_mat,
      converged_by_fold = conv_mat,
      start_objective = start_obj,
      start_converged = start_conv,
      lambda_by_fold = lambda_list,
      rank_min = rank_min,
      rank_1se = rank_1se,
      selected_rank = selected,
      selected = selected,
      rule = rule,
      family = fam,
      family_key = key,
      lambda = lambda,
      lambda_method = lambda_method,
      fit_args = fit_args,
      y = y,
      X = X,
      variable_order = colnames(X),
      offset = offset_full,
      linear = linear_full,
      weights = weights_full,
      knots_source = knots_source,
      fold_knots = fold_knots,
      call = cl,
      timings = list(
        total_s = sum(time_mat, na.rm = TRUE),
        by_rank = rowSums(time_mat, na.rm = TRUE)
      ),
      keep_fits = isTRUE(keep_fits),
      fits = fold_fits
    ),
    class = c("cv.ttps", "tt_rank_selection")
  )
  .tt_attach_cv_aliases(out)
}

#' cv.glmnet-style alias for [tt_rank_select()].
#'
#' @inheritParams tt_rank_select
#' @param type.measure Synonym of `metric` in [tt_rank_select()].
#' @param nfolds Synonym of `folds` (ignored when `foldid` is supplied).
#' @param ... Passed to [tt_rank_select()] (`k`, `lambda`, `control`, …).
#' @return Same as [tt_rank_select()].
#' @export
cv.ttps <- function(y,
                    X,
                    ranks = 1:5,
                    family = stats::gaussian(),
                    type.measure = c("auto", "mse", "rmse", "deviance",
                                     "poisson_deviance", "binomial_deviance",
                                     "logloss"),
                    nfolds = 5L,
                    foldid = NULL,
                    rule = c("1se", "min"),
                    n_starts = 1L,
                    seed = NULL,
                    ...) {
  type.measure <- match.arg(type.measure)
  rule <- match.arg(rule)
  tt_rank_select(
    y = y,
    X = X,
    ranks = ranks,
    family = family,
    folds = nfolds,
    foldid = foldid,
    rule = rule,
    metric = type.measure,
    n_starts = n_starts,
    seed = seed,
    ...
  )
}

#' Refit [ttps()] on all data at the selected TT rank.
#'
#' Uses the fitting arguments stored on a `"tt_rank_selection"` object.
#' Does **not** change [ttps()] defaults: this is the explicit final fit
#' after model selection. Lambda is re-estimated on the full dataset when
#' `lambda = "cGCV"`; fold-specific lambdas are never averaged.
#'
#' By default `compute_edf` is restored to `TRUE` for the final fit (CV
#' itself disables EDF for speed).
#'
#' @param object A `"tt_rank_selection"` / `"cv.ttps"` from [tt_rank_select()]
#'   or [cv.ttps()].
#' @param y,X Optional override data (default: data stored on `object`).
#' @param rank Rank to refit (default: `object$selected_rank`). For [refit()],
#'   may also be `"1se"`, `"min"`, or an integer.
#' @param control Optional [tt_control()] override. If `NULL`, uses the
#'   selection control with `compute_edf = TRUE`.
#' @param ... Passed to [ttps()] (override stored args).
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
  if (is.null(control)) {
    if (is.null(args$control)) args$control <- tt_control()
    args$control$compute_edf <- TRUE
  } else {
    args$control <- control
  }
  extra <- list(...)
  args$y <- y
  args$X <- X
  args$rank <- as.integer(rank)[1L]
  args$offset <- if ("offset" %in% names(extra)) {
    extra$offset
  } else {
    object$offset
  }
  extra$offset <- NULL
  args$linear <- if ("linear" %in% names(extra)) {
    extra$linear
  } else if (!is.null(object$linear)) {
    object$linear
  } else {
    args$linear
  }
  extra$linear <- NULL
  args$weights <- if ("weights" %in% names(extra)) {
    extra$weights
  } else {
    object$weights
  }
  extra$weights <- NULL
  if (length(extra)) args <- utils::modifyList(args, extra)
  do.call(ttps, args)
}

#' Generic full-data refit after CV (S3).
#' @param object A fitted CV / selection object.
#' @param ... Passed to methods.
#' @export
refit <- function(object, ...) {
  UseMethod("refit")
}

#' @rdname tt_rank_refit
#' @export
refit.tt_rank_selection <- function(object,
                                    rank = "1se",
                                    ...) {
  if (is.character(rank)) {
    rank <- match.arg(rank, c("1se", "min"))
    rank <- if (identical(rank, "1se")) object$rank_1se else object$rank_min
  } else {
    rank <- as.integer(rank)[1L]
  }
  tt_rank_refit(object, rank = rank, ...)
}

#' @rdname tt_rank_refit
#' @export
refit.cv.ttps <- function(object, rank = "1se", ...) {
  refit.tt_rank_selection(object, rank = rank, ...)
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

#' Normalize user foldid to consecutive integers 1:K.
#' @keywords internal
#' @noRd
.tt_normalize_foldid <- function(foldid, n) {
  foldid <- as.integer(foldid)
  if (length(foldid) != as.integer(n)) {
    stop("`foldid` must have length n = ", n, ".", call. = FALSE)
  }
  if (anyNA(foldid)) stop("`foldid` contains NA.", call. = FALSE)
  u <- sort(unique(foldid))
  if (length(u) < 2L) {
    stop("`foldid` must define at least 2 folds.", call. = FALSE)
  }
  if (!identical(u, seq_along(u))) {
    foldid <- match(foldid, u)
  }
  as.integer(foldid)
}

#' Attach glmnet-style aliases on a rank-selection object.
#' @keywords internal
#' @noRd
.tt_attach_cv_aliases <- function(obj) {
  obj$cvm <- obj$cv_results$mean_cv
  obj$cvsd <- obj$cv_results$se_cv
  obj$cvup <- obj$cv_results$mean_cv + obj$cv_results$se_cv
  obj$cvlo <- obj$cv_results$mean_cv - obj$cv_results$se_cv
  obj$cvraw <- obj$loss_by_fold
  obj$foldid <- obj$fold_id
  obj$`rank.min` <- obj$rank_min
  obj$`rank.1se` <- obj$rank_1se
  obj
}

#' @keywords internal
#' @noRd
.tt_resolve_cv_metric <- function(metric, family_key) {
  if (identical(metric, "deviance")) {
    if (identical(family_key, "poisson")) return("poisson_deviance")
    if (identical(family_key, "bernoulli")) return("binomial_deviance")
    stop("`metric = \"deviance\"` is not defined for family_key=",
         family_key, call. = FALSE)
  }
  if (!identical(metric, "auto")) return(metric)
  if (identical(family_key, "gaussian")) return("mse")
  if (identical(family_key, "poisson")) return("poisson_deviance")
  if (identical(family_key, "bernoulli")) return("binomial_deviance")
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
  if (identical(metric, "mse")) {
    return(sum(w * (y - mu)^2) / sw)
  }
  if (identical(metric, "rmse")) {
    return(sqrt(sum(w * (y - mu)^2) / sw))
  }
  if (identical(metric, "poisson_deviance") ||
      identical(metric, "binomial_deviance")) {
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
                                  lambda_list, lambda_method, d, p,
                                  start_obj = NULL, start_conv = NULL,
                                  n_starts = 1L) {
  rows <- lapply(seq_along(ranks), function(i) {
    r <- ranks[i]
    losses <- as.numeric(loss_mat[i, ])
    times <- as.numeric(time_mat[i, ])
    conv <- as.logical(conv_mat[i, ])
    # Retain Inf in the mean (failed folds ⇒ Inf CV). SE uses finite folds only.
    losses_mean <- losses
    losses_mean[is.na(losses_mean)] <- Inf
    mean_cv <- mean(losses_mean)
    finite <- is.finite(losses_mean)
    K <- sum(finite)
    sd_cv <- if (K > 1) stats::sd(losses_mean[finite]) else NA_real_
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
    # Multi-start train-objective stability (pooled over folds × starts)
    obj_best <- obj_med <- obj_sd <- start_gap <- start_conv_rate <- NA_real_
    if (!is.null(start_obj) && length(start_obj) >= i) {
      objs <- as.numeric(start_obj[[i]])
      convs <- as.logical(start_conv[[i]])
      finite_o <- is.finite(objs)
      if (any(finite_o)) {
        obj_best <- min(objs[finite_o])
        obj_med <- stats::median(objs[finite_o])
        obj_sd <- if (sum(finite_o) > 1) stats::sd(objs[finite_o]) else 0
        start_gap <- obj_med - obj_best
      }
      start_conv_rate <- mean(convs)
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
      n_finite_folds = as.integer(K),
      n_inf_folds = as.integer(sum(!finite)),
      lambda_mean = if (length(lam_mean) == 1L) lam_mean else mean(lam_mean, na.rm = TRUE),
      n_starts = as.integer(n_starts),
      objective_best = obj_best,
      objective_median = obj_med,
      objective_sd = obj_sd,
      start_gap = start_gap,
      start_convergence_rate = start_conv_rate,
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
#' If SE is NA (fewer than 2 finite folds at r_min), SE is treated as 0.
#' @keywords internal
#' @noRd
.tt_rank_1se <- function(cv_results, rank_min, tol = 1e-12) {
  row_min <- cv_results[cv_results$rank == rank_min, , drop = FALSE]
  se <- row_min$se_cv
  if (length(se) != 1L || !is.finite(se)) se <- 0
  thresh <- row_min$mean_cv + se
  if (!is.finite(thresh)) return(rank_min)
  ok <- is.finite(cv_results$mean_cv) & (cv_results$mean_cv <= thresh + tol)
  if (!any(ok)) return(rank_min)
  min(cv_results$rank[ok])
}

#' @export
print.tt_rank_selection <- function(x, digits = 3, ...) {
  cat("Tensor-Train rank selection\n\n")
  cat(sprintf("Family:             %s\n", x$family$family))
  cat(sprintf("Metric:             %s\n", x$metric))
  cat(sprintf("Folds:              %d\n", x$folds))
  cat(sprintf("Starts per fold:    %d\n", x$n_starts %||% 1L))
  if (!is.null(x$variable_order)) {
    cat(sprintf("Variable order:     %s\n",
                paste(x$variable_order, collapse = ", ")))
  }
  lam_lab <- if (identical(x$lambda_method, "cGCV")) {
    "cGCV"
  } else if (is.numeric(x$lambda)) {
    paste(sprintf("%.4g", x$lambda), collapse = ",")
  } else {
    as.character(x$lambda)
  }
  cat(sprintf("Lambda:             %s\n", lam_lab))
  ks <- x$knots_source %||% if (isTRUE(x$fold_knots)) "fold" else "full_X"
  cat(sprintf("Knots:              %s\n\n",
              switch(ks,
                     full_X = "shared (full X)",
                     fold = "per-fold (train only)",
                     user = "user-supplied",
                     ks)))

  tab <- x$cv_results
  show <- data.frame(
    rank = tab$rank,
    `CV error` = round(tab$mean_cv, digits),
    SE = round(tab$se_cv, digits),
    `TT params` = tab$npar_tt,
    Compression = sprintf("%.1fx", tab$compression),
    check.names = FALSE
  )
  if (!is.null(tab$start_gap) && (x$n_starts %||% 1L) > 1L) {
    show[["start_gap"]] <- round(tab$start_gap, digits)
    show[["start_conv"]] <- round(tab$start_convergence_rate, 2)
  }
  print(show, row.names = FALSE)
  cat("\n")
  cat(sprintf("Minimum-CV rank:    %d\n", x$rank_min))
  cat(sprintf("1-SE rank:          %d\n", x$rank_1se))
  cat(sprintf("Selected rank:      %d\n", x$selected_rank))
  cat(sprintf("Rule:               %s\n",
              if (identical(x$rule, "1se")) "1-SE" else "minimum-CV"))
  cat("Ranks are predictive / working structural capacity, not a true TT rank.\n")
  invisible(x)
}

#' @export
print.cv.ttps <- function(x, digits = 3, ...) {
  print.tt_rank_selection(x, digits = digits, ...)
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
    n_inf = tab$n_inf_folds %||% 0L,
    mean_t = round(tab$mean_time, 3),
    total_t = round(tab$total_time, 3),
    stringsAsFactors = FALSE
  )
  print(show, row.names = FALSE)
  if (!is.null(tab$objective_best) && (x$n_starts %||% 1L) > 1L) {
    cat("\nMulti-start train-objective stability\n")
    stab <- data.frame(
      rank = tab$rank,
      obj_best = round(tab$objective_best, digits),
      obj_median = round(tab$objective_median, digits),
      obj_sd = round(tab$objective_sd, digits),
      start_gap = round(tab$start_gap, digits),
      start_conv = round(tab$start_convergence_rate, 2),
      stringsAsFactors = FALSE
    )
    print(stab, row.names = FALSE)
    cat("start_gap = objective_median - objective_best (large => init-sensitive).\n")
  }
  cat(sprintf("\nTotal CV wall time: %.3fs\n", x$timings$total_s))
  cat("SE is fold-to-fold variability of the CV loss, not a formal rank test.\n")
  cat("Failed folds contribute Inf to mean CV; SE uses finite folds only.\n")
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
         main = "CV error vs rank (predictive structural capacity)",
         ylim = ylim, ...)
    if (any(is.finite(se) & se > 0)) {
      graphics::arrows(r, y - se, r, y + se, length = 0.05, angle = 90, code = 3)
    }
    graphics::abline(v = x$rank_min, col = "gray40", lty = 2)
    graphics::abline(v = x$rank_1se, col = "darkred", lty = 3)
    graphics::legend("topright",
                     legend = c("mean CV", "rank.min", "rank.1se"),
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

#' @export
plot.cv.ttps <- function(x, type = c("error", "compression", "tradeoff"), ...) {
  plot.tt_rank_selection(x, type = type, ...)
}
