# Bernoulli PIRLS stabilization gate (FIX 1 + FIX 2)
#
# Compares:
#   ALS_old  — pre-fix: inconsistent W/z + soft deviance abort (local replica)
#   ALS_new  — package: consistent W/z + true-objective step-halving
#   LBFGS    — direct penalized Bernoulli objective
#
# Usage:
#   Rscript inst/benchmarks/benchmark_bernoulli_stabilization.R --quick
#   Rscript inst/benchmarks/benchmark_bernoulli_stabilization.R --full
#
# Results -> inst/benchmarks/results/bernoulli_stabilization/
# Does NOT overwrite bernoulli_audit/.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (any(grepl("^--full$", args))) "full" else "quick"

local({
  bench_dir <- (function() {
    fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(fa)) return(dirname(normalizePath(sub("^--file=", "", fa[[1]]))))
    of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(of)) return(dirname(normalizePath(of)))
    if (file.exists("helpers.R")) return(normalizePath("."))
    if (file.exists("inst/benchmarks/helpers.R"))
      return(normalizePath("inst/benchmarks"))
    getwd()
  })()
  source(file.path(bench_dir, "helpers.R"), local = TRUE)
  pkg_root <- normalizePath(file.path(bench_dir, "..", ".."))
  if (!file.exists(file.path(pkg_root, "DESCRIPTION")))
    pkg_root <- normalizePath(file.path(bench_dir, ".."))
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    .ttps_bench_ensure_pkg()
  }

  out_dir <- file.path(.ttps_bench_outdir(), "bernoulli_stabilization")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== Bernoulli PIRLS stabilization (", mode, ") ===")
  message("Results -> ", out_dir)

  # ---- config (match audit; full = user §17 n=800) ----
  if (identical(mode, "quick")) {
    n <- 500L; n_te <- 1500L; k <- 6L
    rank_grid <- c(2L, 3L)
    rank_lam_grid <- c(2L, 3L)
    n_init <- 5L
    lambdas <- c(0.1, 1, 10, 100)
  } else {
    n <- 800L; n_te <- 2500L; k <- 6L
    rank_grid <- c(2L, 3L)
    rank_lam_grid <- c(2L, 3L, 4L)
    n_init <- 20L
    lambdas <- c(0.1, 1, 10, 100)
  }

  # Diagnostic-only catastrophic threshold (does NOT alter fits)
  ETA_CATASTROPHIC <- 20

  set.seed(3)
  dat <- simulate_bernoulli(n = n, d = 3L, seed = 3L)
  te <- holdout_bernoulli(n_te = n_te, d = 3L, shift = dat$shift, seed = 99L)
  y <- dat$y; X <- dat$X

  # ------------------------------------------------------------------
  # Legacy ALS (audit-era): floored W, unfloored z; soft deviance abort
  # ------------------------------------------------------------------
  .legacy_als_fit <- function(y, X, init, rank, k, lambda,
                              pirls_maxit = 40L, als_sweeps = 10L,
                              soft_abort_factor = 1.25, tol = 1e-8) {
    t0 <- proc.time()[["elapsed"]]
    d <- ncol(X)
    rank_chain <- tt_rank(rank, d = d)
    bs <- build_marginal_bases(X, k = k, degree = 3)
    basis <- bs$basis
    p <- bs$k
    cores <- init
    pens <- lapply(seq_len(d), function(kk)
      core_penalty(rank_chain[kk], p, rank_chain[kk + 1L], 2L))
    lam <- rep(as.numeric(lambda), d)
    intercept <- init_intercept(binomial(), y)
    eta <- intercept + tt_contraction(cores, basis)
    mu <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12)
    dev <- glm_deviance(binomial(), y, mu)
    n_pirls <- 0L
    n_als <- 0L
    reason <- "maxit"
    ok <- TRUE

    for (it in seq_len(pirls_maxit)) {
      cores_old <- cores
      intercept_old <- intercept
      eta_old <- eta
      dev_old <- dev

      # inconsistent W/z (pre-FIX1)
      mu_w <- pmin(pmax(plogis(eta), 1e-5), 1 - 1e-5)
      var <- mu_w * (1 - mu_w)
      w <- pmax(var, 1e-4)
      z <- eta + (y - mu_w) / var

      for (sw in seq_len(als_sweeps)) {
        zc <- z - intercept
        for (kk in seq_len(d)) {
          L <- left_interfaces(cores, basis)
          R <- right_interfaces(cores, basis)
          Xk <- tt_design_core(L[[kk]], R[[kk]], basis[[kk]])
          sws <- sqrt(pmax(w, 0))
          Xw <- Xk * sws
          yw <- zc * sws
          A <- crossprod(Xw) + lam[kk] * pens[[kk]]
          b <- as.numeric(crossprod(Xw, yw))
          g <- tryCatch(
            as.numeric(solve(A, b)),
            error = function(e) as.numeric(qr.coef(qr(A), b))
          )
          cores[[kk]] <- array(g, c(rank_chain[kk], p, rank_chain[kk + 1L]))
        }
        f <- tt_contraction(cores, basis)
        intercept <- sum(w * (z - f)) / max(sum(w), 1e-12)
        n_als <- n_als + 1L
      }

      eta <- intercept + tt_contraction(cores, basis)
      mu <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12)
      dev <- glm_deviance(binomial(), y, mu)
      n_pirls <- it

      # soft deviance abort (NOT step-halving)
      if (is.finite(dev_old) && is.finite(dev) &&
          dev > soft_abort_factor * max(dev_old, 1e-12)) {
        cores <- cores_old
        intercept <- intercept_old
        eta <- eta_old
        mu <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12)
        dev <- glm_deviance(binomial(), y, mu)
        reason <- "soft deviance abort"
        ok <- TRUE # legacy marked "ok" while stopping early
        break
      }

      rel <- abs(dev_old - dev) / max(1, abs(dev_old))
      if (!is.finite(dev)) {
        ok <- FALSE
        reason <- "non-finite deviance"
        break
      }
      if (it > 2L && rel < tol) {
        reason <- "relative deviance change"
        break
      }
    }

    obj <- tt_glm_penalized_objective(
      y, cores, intercept, basis, pens, lam, binomial()
    )
    list(
      cores = cores, intercept = intercept, lambda = lam,
      eta = eta, mu = mu, deviance = dev,
      n_pirls = n_pirls, n_als_sweeps = n_als,
      n_step_halvings = 0L,
      elapsed = proc.time()[["elapsed"]] - t0,
      converged = isTRUE(ok) && is.finite(obj$value),
      convergence = list(
        overall = isTRUE(ok) && is.finite(obj$value),
        pirls = isTRUE(ok), als = TRUE,
        reason = reason,
        n_pirls = n_pirls, n_als_sweeps = n_als, n_step_halvings = 0L
      ),
      objective = obj$value,
      knots = bs$knots, degree = bs$degree, penalties = pens
    )
  }

  .predict_legacy_link <- function(fit, Xnew) {
    basis <- eval_marginal_bases(Xnew, fit$knots, fit$degree)
    fit$intercept + tt_contraction(fit$cores, basis)
  }

  .is_catastrophic <- function(obj, max_abs_eta, converged, reason = "") {
    early_fail <- grepl("line search failed", reason, fixed = TRUE)
    (!is.finite(obj)) ||
      (is.finite(max_abs_eta) && max_abs_eta > ETA_CATASTROPHIC) ||
      isTRUE(early_fail)
  }

  .metrics <- function(fit, tag, r, lam, seed = NA_integer_, legacy = FALSE) {
    if (legacy) {
      o_val <- fit$objective
      nll <- NA_real_
      pen <- NA_real_
      eta_tr <- fit$eta
      eta_te <- .predict_legacy_link(fit, te$X)
      p_tr <- plogis(eta_tr)
      n_halve <- 0L
      conv <- fit$convergence
      time_s <- fit$elapsed
      n_pirls <- fit$n_pirls
      n_als <- fit$n_als_sweeps
      reason <- conv$reason
      overall <- conv$overall
    } else {
      o <- tt_objective(fit, X, y)
      o_val <- o$value
      nll <- o$nll_or_sse
      pen <- o$penalty
      eta_tr <- fit$linear.predictors
      eta_te <- predict(fit, te$X, type = "link")
      p_tr <- predict(fit, type = "response")
      n_halve <- fit$convergence$n_step_halvings %||% 0L
      time_s <- fit$timing
      n_pirls <- fit$n_pirls %||% NA_real_
      n_als <- fit$convergence$n_als_sweeps %||% NA_real_
      reason <- fit$convergence$reason %||% NA_character_
      overall <- fit$convergence$overall %||% fit$converged
    }
    max_ae <- max(abs(eta_tr))
    data.frame(
      tag = tag, rank = r, lambda = lam, seed = seed,
      objective = o_val, nll = nll, penalty = pen,
      deviance = fit$deviance,
      rmse_eta_train = rmse(eta_tr, dat$truth_eta),
      rmse_eta_test = rmse(eta_te, te$truth_eta),
      max_abs_eta = max_ae,
      min_p = min(p_tr), max_p = max(p_tr),
      time_s = time_s,
      n_pirls = n_pirls,
      n_als_sweeps = n_als,
      n_step_halvings = n_halve,
      converged = isTRUE(overall),
      reason = as.character(reason),
      catastrophic = .is_catastrophic(o_val, max_ae, overall, as.character(reason)),
      stringsAsFactors = FALSE
    )
  }

  .fit_new_als <- function(init, r, lam) {
    ttps(
      y, X, family = binomial(), rank = r, k = k, lambda = lam,
      optimizer = "ALS", init = init,
      control = tt_control(
        backend = "R", pirls_maxit = 40L, als_sweeps_per_pirls = 10L,
        pirls_step_halving = TRUE, compute_edf = FALSE, seed = 123L
      )
    )
  }

  .fit_lbfgs <- function(init, r, lam) {
    ttps(
      y, X, family = binomial(), rank = r, k = k, lambda = lam,
      optimizer = "LBFGS", init = init,
      control = tt_control(
        backend = "R", lbfgs_maxit = 300L, compute_edf = FALSE
      )
    )
  }

  .fit_trio <- function(init, r, lam, seed = NA_integer_) {
    rows <- list()
    leg <- .legacy_als_fit(y, X, init, r, k, lam)
    rows[[1]] <- .metrics(leg, "ALS_old", r, lam, seed = seed, legacy = TRUE)
    fit_n <- .fit_new_als(init, r, lam)
    rows[[2]] <- .metrics(fit_n, "ALS_new", r, lam, seed = seed)
    fit_l <- .fit_lbfgs(init, r, lam)
    rows[[3]] <- .metrics(fit_l, "LBFGS", r, lam, seed = seed)
    list(rows = do.call(rbind, rows), fit_new = fit_n)
  }

  # ==================================================================
  # Exp 1 — original fixed-λ benchmark (shared init)
  # ==================================================================
  message("\n--- Exp1: original benchmark (lambda=1, shared init) ---")
  exp1_rows <- list()
  hist_saved <- list()
  for (r in rank_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    trio <- .fit_trio(init, r, 1)
    exp1_rows[[length(exp1_rows) + 1L]] <- trio$rows
    # save new ALS history for trajectory diagnostics
    if (is.data.frame(trio$fit_new$history) && nrow(trio$fit_new$history)) {
      h <- trio$fit_new$history
      h$rank <- r
      hist_saved[[length(hist_saved) + 1L]] <- h
    }
    for (i in seq_len(nrow(trio$rows))) {
      m <- trio$rows[i, ]
      message(sprintf(
        "  r=%d %-7s | L=%.4g | RMSE_te=%.3f | max|eta|=%.2f | halvings=%s | cat=%s | %s",
        r, m$tag, m$objective, m$rmse_eta_test, m$max_abs_eta,
        m$n_step_halvings, m$catastrophic, m$reason
      ))
    }
  }
  exp1 <- do.call(rbind, exp1_rows)
  write.csv(exp1, file.path(out_dir, "exp1_original_benchmark.csv"), row.names = FALSE)
  if (length(hist_saved)) {
    write.csv(do.call(rbind, hist_saved),
              file.path(out_dir, "exp1_als_new_history.csv"), row.names = FALSE)
  }

  # ==================================================================
  # Exp 2 — multi-init
  # ==================================================================
  message(sprintf("\n--- Exp2: multi-init (%d seeds) ---", n_init))
  multi_rows <- list()
  for (r in rank_grid) {
    for (s in seq_len(n_init)) {
      init <- tt_initialize(X, rank = r, k = k, seed = 2000L + s, sd = 0.05)
      multi_rows[[length(multi_rows) + 1L]] <- .fit_trio(init, r, 1, seed = s)$rows
    }
  }
  multi <- do.call(rbind, multi_rows)
  write.csv(multi, file.path(out_dir, "exp2_multiinit.csv"), row.names = FALSE)

  .summ <- function(sub, col) {
    x <- sub[[col]]
    x <- x[is.finite(x)]
    if (!length(x)) return(c(med = NA, q1 = NA, q3 = NA, mn = NA, mx = NA))
    c(
      med = stats::median(x),
      q1 = as.numeric(stats::quantile(x, 0.25)),
      q3 = as.numeric(stats::quantile(x, 0.75)),
      mn = min(x), mx = max(x)
    )
  }

  summ_rows <- list()
  for (r in rank_grid) {
    for (tag in c("ALS_old", "ALS_new", "LBFGS")) {
      sub <- multi[multi$rank == r & multi$tag == tag, ]
      s_rmse <- .summ(sub, "rmse_eta_test")
      s_obj <- .summ(sub, "objective")
      s_eta <- .summ(sub, "max_abs_eta")
      s_t <- .summ(sub, "time_s")
      cat_rate <- mean(sub$catastrophic)
      summ_rows[[length(summ_rows) + 1L]] <- data.frame(
        rank = r, tag = tag, n = nrow(sub),
        catastrophic_rate = cat_rate,
        rmse_med = s_rmse["med"], rmse_q1 = s_rmse["q1"], rmse_q3 = s_rmse["q3"],
        rmse_min = s_rmse["mn"], rmse_max = s_rmse["mx"],
        obj_med = s_obj["med"], obj_q1 = s_obj["q1"], obj_q3 = s_obj["q3"],
        obj_min = s_obj["mn"], obj_max = s_obj["mx"],
        maxeta_med = s_eta["med"], maxeta_q1 = s_eta["q1"], maxeta_q3 = s_eta["q3"],
        maxeta_min = s_eta["mn"], maxeta_max = s_eta["mx"],
        time_med = s_t["med"], time_q1 = s_t["q1"], time_q3 = s_t["q3"],
        stringsAsFactors = FALSE
      )
      message(sprintf(
        "  r=%d %-7s | cat=%.0f%% | med RMSE=%.3f IQR=[%.3f,%.3f] | med max|eta|=%.2f | med t=%.2fs",
        r, tag, 100 * cat_rate, s_rmse["med"], s_rmse["q1"], s_rmse["q3"],
        s_eta["med"], s_t["med"]
      ))
    }
  }
  write.csv(do.call(rbind, summ_rows),
            file.path(out_dir, "exp2_multiinit_summary.csv"), row.names = FALSE)

  # ==================================================================
  # Exp 3 — rank × lambda
  # ==================================================================
  message("\n--- Exp3: rank × lambda ---")
  lam_rows <- list()
  for (r in rank_lam_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    for (lam in lambdas) {
      lam_rows[[length(lam_rows) + 1L]] <- .fit_trio(init, r, lam)$rows
    }
  }
  lam <- do.call(rbind, lam_rows)
  write.csv(lam, file.path(out_dir, "exp3_rank_lambda.csv"), row.names = FALSE)

  for (r in rank_lam_grid) {
    for (tag in c("ALS_old", "ALS_new", "LBFGS")) {
      sub <- lam[lam$rank == r & lam$tag == tag, ]
      message(sprintf(
        "  r=%d %-7s | cat=%.0f%% | RMSE range [%.3f, %.3f] | max|eta| range [%.2f, %.2f]",
        r, tag, 100 * mean(sub$catastrophic),
        min(sub$rmse_eta_test), max(sub$rmse_eta_test),
        min(sub$max_abs_eta), max(sub$max_abs_eta)
      ))
    }
  }

  # Append hybrid comparison (Outcome B follow-up) onto existing results
  message("\n--- Exp4: hybrid vs ALS_new vs LBFGS (multi-init) ---")
  hy_rows <- list()
  for (r in rank_grid) {
    for (s in seq_len(n_init)) {
      init <- tt_initialize(X, rank = r, k = k, seed = 2000L + s, sd = 0.05)
      for (opt in c("ALS", "hybrid", "LBFGS")) {
        ctrl <- if (identical(opt, "ALS")) {
          tt_control(backend = "R", pirls_maxit = 40L, als_sweeps_per_pirls = 10L,
                     pirls_step_halving = TRUE, compute_edf = FALSE)
        } else if (identical(opt, "hybrid")) {
          tt_control(backend = "R", pirls_maxit = 40L, als_sweeps_per_pirls = 10L,
                     pirls_step_halving = TRUE, hybrid_lbfgs_maxit = 50L,
                     compute_edf = FALSE)
        } else {
          tt_control(backend = "R", lbfgs_maxit = 300L, compute_edf = FALSE)
        }
        fit <- ttps(
          y, X, family = binomial(), rank = r, k = k, lambda = 1,
          optimizer = opt, init = init, control = ctrl
        )
        tag <- if (identical(opt, "ALS")) "ALS_new" else opt
        hy_rows[[length(hy_rows) + 1L]] <- .metrics(fit, tag, r, 1, seed = s)
      }
    }
  }
  hy <- do.call(rbind, hy_rows)
  write.csv(hy, file.path(out_dir, "exp4_hybrid.csv"), row.names = FALSE)
  hy_summ <- list()
  for (r in rank_grid) {
    for (tag in c("ALS_new", "hybrid", "LBFGS")) {
      sub <- hy[hy$rank == r & hy$tag == tag, ]
      s_rmse <- .summ(sub, "rmse_eta_test")
      s_eta <- .summ(sub, "max_abs_eta")
      s_t <- .summ(sub, "time_s")
      s_obj <- .summ(sub, "objective")
      cat_rate <- mean(sub$catastrophic)
      hy_summ[[length(hy_summ) + 1L]] <- data.frame(
        rank = r, tag = tag, n = nrow(sub),
        catastrophic_rate = cat_rate,
        rmse_med = s_rmse["med"], rmse_q1 = s_rmse["q1"], rmse_q3 = s_rmse["q3"],
        obj_med = s_obj["med"],
        maxeta_med = s_eta["med"],
        time_med = s_t["med"],
        stringsAsFactors = FALSE
      )
      message(sprintf(
        "  r=%d %-7s | cat=%.0f%% | med RMSE=%.3f | med max|eta|=%.2f | med t=%.2fs | med L=%.4g",
        r, tag, 100 * cat_rate, s_rmse["med"], s_eta["med"], s_t["med"], s_obj["med"]
      ))
    }
  }
  write.csv(do.call(rbind, hy_summ),
            file.path(out_dir, "exp4_hybrid_summary.csv"), row.names = FALSE)

  # Gate helpers
  gate <- list(
    mode = mode, n = n, k = k, n_init = n_init,
    eta_catastrophic_threshold = ETA_CATASTROPHIC,
    exp1 = exp1,
    multi_summary = do.call(rbind, summ_rows),
    hybrid_summary = do.call(rbind, hy_summ),
    timestamp = as.character(Sys.time())
  )
  saveRDS(gate, file.path(out_dir, "gate_snapshot.rds"))
  writeLines(c(
    paste0("mode=", mode),
    paste0("n=", n),
    paste0("n_init=", n_init),
    paste0("ETA_CATASTROPHIC=", ETA_CATASTROPHIC),
    paste0("timestamp=", gate$timestamp)
  ), file.path(out_dir, "README.txt"))

  message("\n=== Stabilization benchmarks complete ===")
})
