# Diagnostic: ALS vs L-BFGS (Gaussian) — fixed λ then cGCV
#
# Priority checks (not a large benchmark):
#  1) Training penalized objective (same definition)
#  2) Long ALS with strict tol vs LBFGS
#  3) Multi-init (shared seeds)
#  4) Cross warm-start ALS <-> LBFGS
#  5) cGCV λ trajectories
#
#   Rscript inst/benchmarks/diagnose_als_vs_lbfgs.R

local({
  bench_dir <- (function() {
    of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(of)) return(dirname(normalizePath(of)))
    if (file.exists("helpers.R")) return(normalizePath("."))
    if (file.exists("inst/benchmarks/helpers.R"))
      return(normalizePath("inst/benchmarks"))
    system.file("benchmarks", package = "TTPsplines")
  })()
  source(file.path(bench_dir, "helpers.R"), local = TRUE)
  # Always prefer development tree for diagnostics
  pkg_root <- normalizePath(file.path(bench_dir, "..", ".."))
  if (file.exists(file.path(pkg_root, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    .ttps_bench_ensure_pkg()
  }
  out_dir <- .ttps_bench_outdir()

  message("=== Diagnostic: ALS vs LBFGS (Gaussian) ===")
  n <- 800L
  d <- 3L
  k <- 6L
  ranks <- c(2L, 3L)
  lambda_fixed <- 1
  n_init <- as.integer(Sys.getenv("TTPSPLINES_DIAG_NINIT", "12"))
  dat <- simulate_gaussian(n = n, d = d, sigma = 0.3, seed = 21L)
  te <- holdout_gaussian(n_te = 2500L, d = d, seed = 211L)

  .fit_als <- function(init, r, sweeps, tol, lambda) {
    ttpspline(
      dat$y, dat$X, family = gaussian(), rank = r, k = k,
      lambda = lambda, optimizer = "ALS", init = init,
      control = tt_control(
        backend = "R", max_sweeps = sweeps, tol = tol,
        seed = 1L, compute_edf = FALSE, trace = FALSE
      )
    )
  }
  .fit_lbfgs <- function(init, r, maxit, lambda) {
    ttpspline(
      dat$y, dat$X, family = gaussian(), rank = r, k = k,
      lambda = lambda, optimizer = "LBFGS", init = init,
      control = tt_control(
        backend = "R", lbfgs_maxit = maxit,
        seed = 1L, compute_edf = FALSE, trace = FALSE
      )
    )
  }
  .metrics <- function(fit, tag, r, init_id = NA_integer_) {
    obj <- tt_objective(fit, dat$X)
    eta_te <- predict(fit, te$X, type = "link")
    data.frame(
      tag = tag,
      rank = r,
      init_id = init_id,
      objective = obj$value,
      rss = obj$rss,
      penalty = obj$penalty,
      rmse_train_truth = rmse(fit$linear.predictors, dat$truth),
      rmse_test_truth = rmse(eta_te, te$truth),
      time_s = fit$timing,
      n_sweeps = fit$n_sweeps %||% NA_real_,
      n_opt_iter = fit$n_opt_iter %||% NA_real_,
      lambda_geo = exp(mean(log(pmax(fit$lambda, 1e-12)))),
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  path_rows <- list()

  # ------------------------------------------------------------------
  # A. Single shared init: short ALS vs long ALS vs LBFGS (fixed λ)
  # ------------------------------------------------------------------
  message("\n--- A. Fixed λ=1: short / long ALS vs LBFGS ---")
  for (r in ranks) {
    init <- tt_initialize(dat$X, rank = r, k = k, seed = 123L, sd = 0.1)
    als_short <- .fit_als(init, r, sweeps = 20L, tol = 1e-6, lambda = lambda_fixed)
    als_long <- .fit_als(init, r, sweeps = 200L, tol = 1e-12, lambda = lambda_fixed)
    lbfgs <- .fit_lbfgs(init, r, maxit = 500L, lambda = lambda_fixed)

    for (fit in list(als_short, als_long)) {
      if (is.list(fit$history) && length(fit$history)) {
        for (h in fit$history) {
          path_rows[[length(path_rows) + 1L]] <- data.frame(
            phase = "fixed_path",
            optimizer = "ALS",
            rank = r,
            sweep = h$sweep,
            objective = h$objective %||% NA_real_,
            rss = h$rss,
            d_eta = h$d_eta %||% NA_real_,
            lambda_1 = h$lambda[1],
            lambda_2 = h$lambda[2],
            lambda_3 = h$lambda[3],
            stringsAsFactors = FALSE
          )
        }
      }
    }

    rows[[length(rows) + 1L]] <- .metrics(als_short, "ALS_short20", r)
    rows[[length(rows) + 1L]] <- .metrics(als_long, "ALS_long200", r)
    rows[[length(rows) + 1L]] <- .metrics(lbfgs, "LBFGS", r)

    # Cross warm-start
    message(sprintf("  r=%d cross warm-start ...", r))
    als_from_lbfgs <- .fit_als(lbfgs$cores, r, sweeps = 100L, tol = 1e-12, lambda = lambda_fixed)
    lbfgs_from_als <- .fit_lbfgs(als_long$cores, r, maxit = 500L, lambda = lambda_fixed)
    rows[[length(rows) + 1L]] <- .metrics(als_from_lbfgs, "ALS_from_LBFGS", r)
    rows[[length(rows) + 1L]] <- .metrics(lbfgs_from_als, "LBFGS_from_ALS", r)

    o_als <- tt_objective(als_long, dat$X)$value
    o_lbf <- tt_objective(lbfgs, dat$X)$value
    o_l_from_a <- tt_objective(lbfgs_from_als, dat$X)$value
    o_a_from_l <- tt_objective(als_from_lbfgs, dat$X)$value
    message(sprintf(
      "  r=%d | obj ALS_long=%.6g | LBFGS=%.6g | LBFGS<-ALS=%.6g | ALS<-LBFGS=%.6g",
      r, o_als, o_lbf, o_l_from_a, o_a_from_l
    ))
    message(sprintf(
      "  r=%d | Δobj(LBFGS - ALS_long)=%.4g | Δobj(LBFGS_from_ALS - ALS_long)=%.4g",
      r, o_lbf - o_als, o_l_from_a - o_als
    ))
  }

  # ------------------------------------------------------------------
  # B. Multi-init fixed λ
  # ------------------------------------------------------------------
  message(sprintf("\n--- B. Multi-init fixed λ (%d seeds) ---", n_init))
  multi <- list()
  for (r in ranks) {
    for (s in seq_len(n_init)) {
      init <- tt_initialize(dat$X, rank = r, k = k, seed = 1000L + s, sd = 0.1)
      als <- .fit_als(init, r, sweeps = 80L, tol = 1e-10, lambda = lambda_fixed)
      lbf <- .fit_lbfgs(init, r, maxit = 400L, lambda = lambda_fixed)
      ma <- .metrics(als, "ALS", r, init_id = s)
      ml <- .metrics(lbf, "LBFGS", r, init_id = s)
      multi[[length(multi) + 1L]] <- ma
      multi[[length(multi) + 1L]] <- ml
      rows[[length(rows) + 1L]] <- transform(ma, tag = paste0("multi_ALS"))
      rows[[length(rows) + 1L]] <- transform(ml, tag = paste0("multi_LBFGS"))
    }
  }
  multi_tab <- do.call(rbind, multi)
  for (r in ranks) {
    sub <- multi_tab[multi_tab$rank == r, ]
    als_o <- sub$objective[sub$tag == "ALS"]
    lbf_o <- sub$objective[sub$tag == "LBFGS"]
    als_te <- sub$rmse_test_truth[sub$tag == "ALS"]
    lbf_te <- sub$rmse_test_truth[sub$tag == "LBFGS"]
    message(sprintf(
      "  r=%d | obj ALS median=%.4g (min=%.4g) | LBFGS median=%.4g (min=%.4g)",
      r, median(als_o), min(als_o), median(lbf_o), min(lbf_o)
    ))
    message(sprintf(
      "  r=%d | RMSE_te ALS median=%.4f | LBFGS median=%.4f | ALS<LBFGS obj in %d/%d inits",
      r, median(als_te), median(lbf_te), sum(als_o < lbf_o), n_init
    ))
    message(sprintf(
      "  r=%d | mean(obj_LBFGS - obj_ALS)=%.4g",
      r, mean(lbf_o - als_o)
    ))
  }

  # ------------------------------------------------------------------
  # C. cGCV paths (shared init)
  # ------------------------------------------------------------------
  message("\n--- C. cGCV λ trajectories (shared init) ---")
  for (r in ranks) {
    init <- tt_initialize(dat$X, rank = r, k = k, seed = 123L, sd = 0.1)
    als_c <- ttpspline(
      dat$y, dat$X, family = gaussian(), rank = r, k = k,
      lambda = "cGCV", optimizer = "ALS", init = init,
      control = tt_control(
        backend = "R", max_sweeps = 30L, tol_lambda = 1e-4,
        lambda_bounds = c(1e-2, 1e2), compute_edf = FALSE
      )
    )
    lbf_c <- ttpspline(
      dat$y, dat$X, family = gaussian(), rank = r, k = k,
      lambda = "cGCV", optimizer = "LBFGS", init = init,
      control = tt_control(
        backend = "R", lbfgs_maxit = 150L, outer_maxit = 12L,
        outer_tol = 1e-4, lambda_bounds = c(1e-2, 1e2), compute_edf = FALSE
      )
    )
    rows[[length(rows) + 1L]] <- .metrics(als_c, "ALS_cGCV", r)
    rows[[length(rows) + 1L]] <- .metrics(lbf_c, "LBFGS_cGCV", r)

    if (is.list(als_c$history) && length(als_c$history)) {
      for (h in als_c$history) {
        path_rows[[length(path_rows) + 1L]] <- data.frame(
          phase = "cGCV_ALS",
          optimizer = "ALS",
          rank = r,
          sweep = h$sweep,
          objective = h$objective %||% NA_real_,
          rss = h$rss,
          d_eta = h$d_eta %||% NA_real_,
          lambda_1 = h$lambda[1],
          lambda_2 = h$lambda[2],
          lambda_3 = h$lambda[3],
          stringsAsFactors = FALSE
        )
      }
    }
    if (is.list(lbf_c$history) && length(lbf_c$history)) {
      for (h in lbf_c$history) {
        path_rows[[length(path_rows) + 1L]] <- data.frame(
          phase = "cGCV_LBFGS",
          optimizer = "LBFGS",
          rank = r,
          sweep = h$outer %||% NA_integer_,
          objective = h$value %||% NA_real_,
          rss = NA_real_,
          d_eta = NA_real_,
          lambda_1 = h$lambda[1],
          lambda_2 = h$lambda[2],
          lambda_3 = h$lambda[3],
          stringsAsFactors = FALSE
        )
      }
    }
    message(sprintf(
      "  r=%d | ALS-cGCV λ=(%s) RMSE_te=%.4f | LBFGS-cGCV λ=(%s) RMSE_te=%.4f",
      r,
      paste(sprintf("%.3g", als_c$lambda), collapse = ","),
      rmse(predict(als_c, te$X), te$truth),
      paste(sprintf("%.3g", lbf_c$lambda), collapse = ","),
      rmse(predict(lbf_c, te$X), te$truth)
    ))
  }

  tab <- do.call(rbind, rows)
  paths <- if (length(path_rows)) do.call(rbind, path_rows) else NULL
  .write_bench_csv(tab, "diagnose_als_vs_lbfgs.csv", out_dir)
  if (!is.null(paths)) {
    .write_bench_csv(paths, "diagnose_als_vs_lbfgs_paths.csv", out_dir)
  }
  .write_bench_csv(multi_tab, "diagnose_als_vs_lbfgs_multiinit.csv", out_dir)

  # Compact summary printed at end
  message("\n=== Summary table (fixed λ single-init + cross) ===")
  print(subset(tab, grepl("^(ALS_|LBFGS)", tag) & !grepl("cGCV|multi", tag))[
    , c("tag", "rank", "objective", "rss", "rmse_test_truth", "time_s")
  ], row.names = FALSE, digits = 5)

  message("\n=== cGCV summary ===")
  print(subset(tab, grepl("cGCV", tag))[
    , c("tag", "rank", "objective", "rss", "rmse_test_truth", "lambda_geo", "time_s")
  ], row.names = FALSE, digits = 5)

  invisible(list(summary = tab, paths = paths, multi = multi_tab))
})
