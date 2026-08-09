# Bernoulli: PIRLS-ALS vs GD vs LBFGS (fixed λ)
#
# Scientific contrast (docs/OPTIMIZER_GD.md):
#   GD ≈ LBFGS ≫ PIRLS-ALS  →  PIRLS/alternating path is the issue
#   LBFGS ≫ GD ≈ PIRLS-ALS →  quasi-Newton curvature matters
#
# Usage:
#   Rscript inst/benchmarks/benchmark_bernoulli_gd.R --quick
#   Rscript inst/benchmarks/benchmark_bernoulli_gd.R --full
#
# Results -> inst/benchmarks/results/bernoulli_gd/

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

  out_dir <- file.path(.ttps_bench_outdir(), "bernoulli_gd")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== Bernoulli GD vs LBFGS vs PIRLS-ALS (", mode, ") ===")
  message("Results -> ", out_dir)

  if (identical(mode, "quick")) {
    n <- 500L; n_te <- 1500L; k <- 6L
    rank_grid <- c(2L, 3L)
    rank_lam <- c(2L, 3L)
    n_init <- 5L
    lambdas <- c(0.1, 1, 10, 100)
    gd_maxit <- 1500L
  } else {
    n <- 800L; n_te <- 2500L; k <- 6L
    rank_grid <- c(2L, 3L)
    rank_lam <- c(2L, 3L, 4L)
    n_init <- 20L
    lambdas <- c(0.1, 1, 10, 100)
    gd_maxit <- 3000L
  }

  ETA_CATASTROPHIC <- 20
  opts <- c("PIRLS-ALS", "GD", "LBFGS")

  set.seed(3)
  dat <- simulate_bernoulli(n = n, d = 3L, seed = 3L)
  te <- holdout_bernoulli(n_te = n_te, d = 3L, shift = dat$shift, seed = 99L)
  y <- dat$y; X <- dat$X

  .ctrl <- function(opt) {
    if (identical(opt, "PIRLS-ALS")) {
      tt_control(
        backend = "R", pirls_maxit = 40L, als_sweeps_per_pirls = 10L,
        pirls_step_halving = TRUE, compute_edf = FALSE
      )
    } else if (identical(opt, "GD")) {
      tt_control(
        backend = "R",
        gd_lr = 1, gd_maxit = gd_maxit, gd_tol = 1e-7,
        gd_linesearch = TRUE, compute_edf = FALSE
      )
    } else {
      tt_control(backend = "R", lbfgs_maxit = 300L, compute_edf = FALSE)
    }
  }

  .is_cat <- function(obj, max_ae) {
    (!is.finite(obj)) || (is.finite(max_ae) && max_ae > ETA_CATASTROPHIC)
  }

  .metrics <- function(fit, tag, r, lam, seed = NA_integer_) {
    o <- tt_objective(fit, X, y)
    eta_te <- predict(fit, te$X, type = "link")
    p_tr <- predict(fit, type = "response")
    max_ae <- max(abs(fit$linear.predictors))
    data.frame(
      tag = tag, rank = r, lambda = lam, seed = seed,
      objective = o$value, nll = o$nll_or_sse, penalty = o$penalty,
      deviance = fit$deviance,
      rmse_eta_train = rmse(fit$linear.predictors, dat$truth_eta),
      rmse_eta_test = rmse(eta_te, te$truth_eta),
      max_abs_eta = max_ae,
      min_p = min(p_tr), max_p = max(p_tr),
      time_s = fit$timing,
      n_opt_iter = fit$n_opt_iter %||% NA_real_,
      n_pirls = fit$n_pirls %||% NA_real_,
      converged = isTRUE(fit$convergence$overall %||% fit$converged),
      reason = as.character(fit$convergence$reason %||% NA_character_),
      catastrophic = .is_cat(o$value, max_ae),
      stringsAsFactors = FALSE
    )
  }

  .fit_one <- function(init, opt, r, lam, seed = NA_integer_) {
    fit <- ttpspline(
      y, X, family = binomial(), rank = r, k = k, lambda = lam,
      optimizer = opt, init = init, control = .ctrl(opt)
    )
    .metrics(fit, opt, r, lam, seed = seed)
  }

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

  # ---- Exp1: shared init, λ=1 ----
  message("\n--- Exp1: shared init, lambda=1 ---")
  exp1_rows <- list()
  for (r in rank_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    for (opt in opts) {
      m <- .fit_one(init, opt, r, 1)
      exp1_rows[[length(exp1_rows) + 1L]] <- m
      message(sprintf(
        "  r=%d %-10s | L=%.4g | RMSE_te=%.3f | max|eta|=%.2f | t=%.2fs | cat=%s",
        r, opt, m$objective, m$rmse_eta_test, m$max_abs_eta, m$time_s, m$catastrophic
      ))
    }
  }
  exp1 <- do.call(rbind, exp1_rows)
  write.csv(exp1, file.path(out_dir, "exp1_shared_init.csv"), row.names = FALSE)

  # ---- Exp2: multi-init ----
  message(sprintf("\n--- Exp2: multi-init (%d seeds) ---", n_init))
  multi_rows <- list()
  for (r in rank_grid) {
    for (s in seq_len(n_init)) {
      init <- tt_initialize(X, rank = r, k = k, seed = 2000L + s, sd = 0.05)
      for (opt in opts) {
        multi_rows[[length(multi_rows) + 1L]] <- .fit_one(init, opt, r, 1, seed = s)
      }
    }
  }
  multi <- do.call(rbind, multi_rows)
  write.csv(multi, file.path(out_dir, "exp2_multiinit.csv"), row.names = FALSE)

  summ_rows <- list()
  for (r in rank_grid) {
    for (opt in opts) {
      sub <- multi[multi$rank == r & multi$tag == opt, ]
      s_rmse <- .summ(sub, "rmse_eta_test")
      s_obj <- .summ(sub, "objective")
      s_eta <- .summ(sub, "max_abs_eta")
      s_t <- .summ(sub, "time_s")
      cat_rate <- mean(sub$catastrophic)
      summ_rows[[length(summ_rows) + 1L]] <- data.frame(
        rank = r, tag = opt, n = nrow(sub),
        catastrophic_rate = cat_rate,
        rmse_med = s_rmse["med"], rmse_q1 = s_rmse["q1"], rmse_q3 = s_rmse["q3"],
        rmse_min = s_rmse["mn"], rmse_max = s_rmse["mx"],
        obj_med = s_obj["med"], obj_q1 = s_obj["q1"], obj_q3 = s_obj["q3"],
        maxeta_med = s_eta["med"], maxeta_q1 = s_eta["q1"], maxeta_q3 = s_eta["q3"],
        maxeta_max = s_eta["mx"],
        time_med = s_t["med"], time_q1 = s_t["q1"], time_q3 = s_t["q3"],
        stringsAsFactors = FALSE
      )
      message(sprintf(
        "  r=%d %-10s | cat=%.0f%% | med RMSE=%.3f IQR=[%.3f,%.3f] | med max|eta|=%.2f | med t=%.2fs | med L=%.4g",
        r, opt, 100 * cat_rate, s_rmse["med"], s_rmse["q1"], s_rmse["q3"],
        s_eta["med"], s_t["med"], s_obj["med"]
      ))
    }
  }
  summ <- do.call(rbind, summ_rows)
  write.csv(summ, file.path(out_dir, "exp2_multiinit_summary.csv"), row.names = FALSE)

  # ---- Exp3: rank × λ ----
  message("\n--- Exp3: rank × lambda ---")
  lam_rows <- list()
  for (r in rank_lam) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    for (lam in lambdas) {
      for (opt in opts) {
        lam_rows[[length(lam_rows) + 1L]] <- .fit_one(init, opt, r, lam)
      }
    }
  }
  lam <- do.call(rbind, lam_rows)
  write.csv(lam, file.path(out_dir, "exp3_rank_lambda.csv"), row.names = FALSE)

  for (r in rank_lam) {
    for (opt in opts) {
      sub <- lam[lam$rank == r & lam$tag == opt, ]
      message(sprintf(
        "  r=%d %-10s | cat=%.0f%% | RMSE [%.3f, %.3f] | max|eta| [%.2f, %.2f]",
        r, opt, 100 * mean(sub$catastrophic),
        min(sub$rmse_eta_test, na.rm = TRUE), max(sub$rmse_eta_test, na.rm = TRUE),
        min(sub$max_abs_eta, na.rm = TRUE), max(sub$max_abs_eta, na.rm = TRUE)
      ))
    }
  }

  # ---- Decision snapshot ----
  # Compare median test RMSE: GD vs LBFGS vs PIRLS at each rank
  decide <- function(r) {
    sub <- summ[summ$rank == r, ]
    getv <- function(tag, col) sub[sub$tag == tag, col]
    pirls <- getv("PIRLS-ALS", "rmse_med")
    gd <- getv("GD", "rmse_med")
    lb <- getv("LBFGS", "rmse_med")
    # relative gaps
    gap_gd_lb <- (gd - lb) / max(lb, 1e-8)
    gap_pirls_lb <- (pirls - lb) / max(lb, 1e-8)
    if (isTRUE(gap_gd_lb < 0.25) && isTRUE(gap_pirls_lb > 0.5)) {
      "GD ≈ LBFGS ≫ PIRLS-ALS (path issue)"
    } else if (isTRUE(gap_gd_lb > 0.5) && isTRUE(abs(pirls - gd) / max(gd, 1e-8) < 0.35)) {
      "LBFGS ≫ GD ≈ PIRLS-ALS (quasi-Newton matters)"
    } else {
      "mixed / inconclusive at this rank"
    }
  }
  decisions <- data.frame(
    rank = rank_grid,
    interpretation = vapply(rank_grid, decide, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(decisions, file.path(out_dir, "decision_snapshot.csv"), row.names = FALSE)
  message("\n--- Decision snapshot ---")
  print(decisions, row.names = FALSE)

  saveRDS(
    list(mode = mode, n = n, n_init = n_init, gd_maxit = gd_maxit,
         ETA_CATASTROPHIC = ETA_CATASTROPHIC, exp1 = exp1, summ = summ,
         decisions = decisions, timestamp = as.character(Sys.time())),
    file.path(out_dir, "gate_snapshot.rds")
  )
  writeLines(c(
    paste0("mode=", mode),
    paste0("n=", n),
    paste0("n_init=", n_init),
    paste0("gd_maxit=", gd_maxit),
    paste0("ETA_CATASTROPHIC=", ETA_CATASTROPHIC),
    paste0("timestamp=", Sys.time())
  ), file.path(out_dir, "README.txt"))

  message("\n=== Bernoulli GD benchmark complete ===")
})
