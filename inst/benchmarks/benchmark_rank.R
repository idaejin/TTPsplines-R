# Paper 1 — rank profile diagnostic (tt_rank_profile)
# source("inst/benchmarks/benchmark_rank.R")

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

  message("=== TTPsplines benchmark: rank profile ===")
  dat <- simulate_gaussian(n = 1200L, d = 3L, sigma = 0.3, seed = 41L)
  te <- holdout_gaussian(n_te = 3000L, d = 3L, seed = 411L)

  ranks <- c(1L, 2L, 3L, 4L, 6L)
  ctrl <- .default_control(backend = "auto", max_sweeps = 12L, seed = 5L)

  # Use package helper + add test RMSE
  prof <- tt_rank_profile(
    dat$y, dat$X,
    ranks = ranks,
    family = gaussian(),
    lambda = "cGCV",
    k = 8L,
    control = ctrl
  )
  tab <- prof$table
  tab$rmse_test_truth <- NA_real_
  for (i in seq_len(nrow(tab))) {
    fit <- ttpspline(
      dat$y, dat$X, family = gaussian(), rank = tab$rank[i], k = 8L,
      lambda = "cGCV", control = ctrl
    )
    tab$rmse_test_truth[i] <- rmse(predict(fit, te$X), te$truth)
  }

  print(tab, row.names = FALSE)
  .write_bench_csv(tab, "benchmark_rank.csv", out_dir)

  fig <- file.path(out_dir, "benchmark_rank.png")
  png(fig, width = 900, height = 400)
  plot(prof)
  # overlay test RMSE in a second device panel already used by plot.tt_rank_profile
  # rewrite custom
  dev.off()
  png(fig, width = 900, height = 400)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  plot(tab$rank, tab$deviance, type = "b", pch = 19,
       xlab = "rank", ylab = "train RSS", main = "Train fit")
  plot(tab$rank, tab$rmse_test_truth, type = "b", pch = 19,
       xlab = "rank", ylab = "test RMSE", main = "Test vs truth")
  par(op)
  dev.off()
  message("Wrote ", fig)
})
