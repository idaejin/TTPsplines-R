# Paper 1 — Bernoulli / binomial scattered-data benchmark
# source("inst/benchmarks/benchmark_bernoulli.R")

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

  message("=== TTPsplines benchmark: Bernoulli (scattered) ===")
  dat <- simulate_bernoulli(n = 2000L, d = 3L, seed = 31L)
  te <- holdout_bernoulli(n_te = 4000L, d = 3L, shift = dat$shift, seed = 311L)

  ranks <- c(1L, 2L, 3L, 4L)
  # Prefer slightly stronger fixed λ for Bernoulli stability
  settings <- list(
    fixed = list(lambda = 5),
    cGCV = list(lambda = "cGCV")
  )
  k <- 8L
  rows <- list()

  for (nm in names(settings)) {
    for (r in ranks) {
      ctrl <- .default_control(
        backend = "auto", pirls_maxit = 15L, als_sweeps_per_pirls = 3L,
        damping = TRUE, seed = 9L,
        lambda_bounds = c(1e-1, 1e2)
      )
      t0 <- proc.time()[["elapsed"]]
      fit <- ttpspline(
        dat$y, dat$X,
        family = binomial(),
        rank = r,
        k = k,
        lambda = settings[[nm]]$lambda,
        control = ctrl
      )
      elapsed <- proc.time()[["elapsed"]] - t0
      eta_te <- predict(fit, te$X, type = "link")
      p_te <- predict(fit, te$X, type = "response")
      rows[[length(rows) + 1L]] <- data.frame(
        family = "bernoulli",
        lambda_mode = nm,
        rank = r,
        k = k,
        n = length(dat$y),
        npar_tt = fit$npar_tt,
        npar_dense = fit$npar_dense,
        compression = fit$compression_ratio,
        deviance = fit$deviance,
        rmse_p_train = rmse(fitted(fit), dat$truth_p),
        rmse_p_test = rmse(p_te, te$truth_p),
        brier_test = mean((p_te - te$truth_p)^2),
        logloss_vs_truth = {
          p <- pmin(pmax(p_te, 1e-12), 1 - 1e-12)
          -mean(te$truth_p * log(p) + (1 - te$truth_p) * log(1 - p))
        },
        cor_eta_test = suppressWarnings(cor(eta_te, te$truth_eta)),
        time_s = elapsed,
        backend = fit$backend,
        n_pirls = fit$n_pirls,
        max_abs_eta = max(abs(fit$linear.predictors)),
        stringsAsFactors = FALSE
      )
      message(sprintf(
        "  λ=%-5s r=%d | cor_η=%.3f | RMSE_p=%.3f | %.2fs",
        nm, r, suppressWarnings(cor(eta_te, te$truth_eta)),
        rmse(p_te, te$truth_p), elapsed
      ))
    }
  }

  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE)
  .write_bench_csv(tab, "benchmark_bernoulli.csv", out_dir)

  fig <- file.path(out_dir, "benchmark_bernoulli_rank.png")
  png(fig, width = 900, height = 400)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  for (lm in c("fixed", "cGCV")) {
    sub <- tab[tab$lambda_mode == lm, ]
    plot(sub$rank, sub$rmse_p_test, type = "b", pch = 19,
         xlab = "rank r", ylab = "RMSE(p̂)", main = paste("Bernoulli —", lm))
  }
  par(op)
  dev.off()
  message("Wrote ", fig)
})
