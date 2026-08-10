# Gaussian scattered-data benchmark
#
# Usage (from package root after load_all / install):
#   source("inst/benchmarks/benchmark_gaussian.R")
# Or:
#   Sys.setenv(TTPSPLINES_BENCH_OUT = "path/to/out")
#   source(system.file("benchmarks/benchmark_gaussian.R", package = "TTPsplines"))
#
# NOT run during R CMD check.

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
  .ttps_bench_ensure_pkg()
  out_dir <- .ttps_bench_outdir()

  message("=== TTPsplines benchmark: Gaussian (scattered) ===")
  dat <- simulate_gaussian(n = 1500L, d = 3L, sigma = 0.3, seed = 11L)
  te <- holdout_gaussian(n_te = 4000L, d = 3L, seed = 111L)

  ranks <- c(1L, 2L, 3L, 4L, 6L)
  lambdas <- list(fixed = 1, cGCV = "cGCV")
  k <- 8L
  rows <- list()

  for (lam_name in names(lambdas)) {
    for (r in ranks) {
      ctrl <- .default_control(backend = "auto", max_sweeps = 12L, seed = 7L)
      t0 <- proc.time()[["elapsed"]]
      fit <- ttps(
        dat$y, dat$X,
        family = gaussian(),
        rank = r,
        k = k,
        lambda = lambdas[[lam_name]],
        control = ctrl
      )
      elapsed <- proc.time()[["elapsed"]] - t0
      eta_te <- predict(fit, te$X, type = "link")
      rows[[length(rows) + 1L]] <- data.frame(
        family = "gaussian",
        lambda_mode = lam_name,
        rank = r,
        k = k,
        n = length(dat$y),
        d = ncol(dat$X),
        npar_tt = fit$npar_tt,
        npar_dense = fit$npar_dense,
        compression = fit$compression_ratio,
        rss_train = fit$deviance,
        rmse_train_truth = rmse(fitted(fit), dat$truth),
        rmse_test_truth = rmse(eta_te, te$truth),
        cor_test = suppressWarnings(cor(eta_te, te$truth)),
        time_s = elapsed,
        backend = fit$backend,
        lambda_geo = exp(mean(log(pmax(fit$lambda, 1e-12)))),
        stringsAsFactors = FALSE
      )
      message(sprintf(
        "  λ=%-5s r=%d | RMSE_te=%.4f | npar=%d | %.2fs",
        lam_name, r, rmse(eta_te, te$truth), fit$npar_tt, elapsed
      ))
    }
  }

  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE)
  .write_bench_csv(tab, "benchmark_gaussian.csv", out_dir)

  # Quick figure
  fig <- file.path(out_dir, "benchmark_gaussian_rank.png")
  png(fig, width = 900, height = 400)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  for (lm in c("fixed", "cGCV")) {
    sub <- tab[tab$lambda_mode == lm, ]
    plot(sub$rank, sub$rmse_test_truth, type = "b", pch = 19,
         xlab = "rank r", ylab = "RMSE vs truth (test)",
         main = paste("Gaussian —", lm))
  }
  par(op)
  dev.off()
  message("Wrote ", fig)
})
