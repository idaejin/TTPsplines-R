# Poisson scattered-data benchmark
# source("inst/benchmarks/benchmark_poisson.R")

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

  message("=== TTPsplines benchmark: Poisson (scattered) ===")
  dat <- simulate_poisson(n = 1500L, d = 3L, seed = 21L)
  te <- holdout_poisson(n_te = 4000L, d = 3L, shift = dat$shift, seed = 211L)

  ranks <- c(1L, 2L, 3L, 4L)
  settings <- list(
    fixed = list(lambda = 1),
    cGCV = list(lambda = "cGCV")
  )
  k <- 8L
  rows <- list()

  for (nm in names(settings)) {
    for (r in ranks) {
      ctrl <- .default_control(
        backend = "auto", pirls_maxit = 15L, als_sweeps_per_pirls = 3L, seed = 8L
      )
      t0 <- proc.time()[["elapsed"]]
      fit <- ttpspline(
        dat$y, dat$X,
        family = poisson(),
        rank = r,
        k = k,
        lambda = settings[[nm]]$lambda,
        control = ctrl
      )
      elapsed <- proc.time()[["elapsed"]] - t0
      eta_te <- predict(fit, te$X, type = "link")
      mu_te <- predict(fit, te$X, type = "response")
      rows[[length(rows) + 1L]] <- data.frame(
        family = "poisson",
        lambda_mode = nm,
        rank = r,
        k = k,
        n = length(dat$y),
        npar_tt = fit$npar_tt,
        npar_dense = fit$npar_dense,
        compression = fit$compression_ratio,
        deviance = fit$deviance,
        rmse_mu_train = rmse(fitted(fit), dat$truth_mu),
        rmse_mu_test = rmse(mu_te, te$truth_mu),
        cor_eta_test = suppressWarnings(cor(eta_te, te$truth_eta)),
        time_s = elapsed,
        backend = fit$backend,
        n_pirls = fit$n_pirls,
        stringsAsFactors = FALSE
      )
      message(sprintf(
        "  λ=%-5s r=%d | cor_η=%.3f | RMSE_μ=%.3f | %.2fs",
        nm, r, suppressWarnings(cor(eta_te, te$truth_eta)),
        rmse(mu_te, te$truth_mu), elapsed
      ))
    }
  }

  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE)
  .write_bench_csv(tab, "benchmark_poisson.csv", out_dir)

  fig <- file.path(out_dir, "benchmark_poisson_rank.png")
  png(fig, width = 900, height = 400)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  for (lm in c("fixed", "cGCV")) {
    sub <- tab[tab$lambda_mode == lm, ]
    plot(sub$rank, sub$cor_eta_test, type = "b", pch = 19, ylim = c(0, 1),
         xlab = "rank r", ylab = "cor(η̂, η)", main = paste("Poisson —", lm))
  }
  par(op)
  dev.off()
  message("Wrote ", fig)
})
