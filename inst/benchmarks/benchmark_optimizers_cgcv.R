# Optimizer comparison — cGCV (same conditional GCV definition)
#
# ALS: cGCV on each core visit.
# LBFGS: outer alternation (optimize cores → conditional cGCV → repeat).
# Adam: stub / unavailable.
#
#   Rscript inst/benchmarks/benchmark_optimizers_cgcv.R

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

  message("=== TTPsplines benchmark: optimizers + cGCV ===")

  n_tr <- as.integer(Sys.getenv("TTPSPLINES_BENCH_N", "800"))
  k <- as.integer(Sys.getenv("TTPSPLINES_BENCH_K", "6"))
  ranks <- as.integer(strsplit(Sys.getenv("TTPSPLINES_BENCH_RANKS", "2,3"), ",")[[1]])
  optimizers <- c("ALS", "LBFGS", "Adam")
  seed_data <- 31L
  seed_init <- 123L

  # Focus on Gaussian + Poisson for cGCV outer cost; Bernoulli optional via env
  fams <- strsplit(Sys.getenv("TTPSPLINES_BENCH_FAMS", "gaussian,poisson"), ",")[[1]]
  fams <- trimws(fams)

  sim_fun <- list(
    gaussian = function() simulate_gaussian(n = n_tr, d = 3L, sigma = 0.3, seed = seed_data),
    poisson = function() simulate_poisson(n = n_tr, d = 3L, seed = seed_data + 1L),
    bernoulli = function() simulate_bernoulli(n = n_tr, d = 3L, seed = seed_data + 2L)
  )
  hold_fun <- list(
    gaussian = function(dat) holdout_gaussian(n_te = 2500L, d = 3L, seed = 311L),
    poisson = function(dat) holdout_poisson(n_te = 2500L, d = 3L, shift = dat$shift, seed = 312L),
    bernoulli = function(dat) holdout_bernoulli(n_te = 2500L, d = 3L, shift = dat$shift, seed = 313L)
  )
  fam_obj <- list(gaussian = gaussian(), poisson = poisson(), bernoulli = binomial())

  rows <- list()
  for (fam_name in fams) {
    dat <- sim_fun[[fam_name]]()
    te <- hold_fun[[fam_name]](dat)
    for (r in ranks) {
      init <- tt_initialize(dat$X, rank = r, k = k, seed = seed_init, sd = 0.1)
      for (opt in optimizers) {
        ctrl <- .default_control(
          backend = "R",
          max_sweeps = 12L,
          pirls_maxit = 12L,
          als_sweeps_per_pirls = 2L,
          lbfgs_maxit = 120L,
          outer_maxit = 8L,
          outer_tol = 1e-4,
          lambda_bounds = c(1e-2, 1e2),
          seed = seed_init,
          init_sd = 0.1,
          compute_edf = TRUE
        )
        t0 <- proc.time()[["elapsed"]]
        fit <- tryCatch(
          ttps(
            dat$y, dat$X,
            family = fam_obj[[fam_name]],
            rank = r,
            k = k,
            lambda = "cGCV",
            optimizer = opt,
            init = init,
            control = ctrl
          ),
          error = function(e) e
        )
        elapsed <- proc.time()[["elapsed"]] - t0

        if (inherits(fit, "error")) {
          rows[[length(rows) + 1L]] <- data.frame(
            setting = "cGCV",
            family = fam_name,
            optimizer = opt,
            rank = r,
            k = k,
            n = length(dat$y),
            status = "error",
            message = conditionMessage(fit),
            lambda_geo = NA_real_,
            lambda_1 = NA_real_,
            lambda_2 = NA_real_,
            lambda_3 = NA_real_,
            npar_tt = NA_real_,
            edf = NA_real_,
            deviance = NA_real_,
            rmse_eta_test = NA_real_,
            cor_eta_test = NA_real_,
            time_s = elapsed,
            n_criterion_evals = NA_real_,
            n_outer = NA_real_,
            n_opt_iter = NA_real_,
            backend = NA_character_,
            converged = FALSE,
            stringsAsFactors = FALSE
          )
          message(sprintf("  %-9s %-7s r=%d | ERROR", fam_name, opt, r))
          next
        }

        eta_te <- predict(fit, te$X, type = "link")
        truth_eta <- if (identical(fam_name, "gaussian")) te$truth else te$truth_eta
        rows[[length(rows) + 1L]] <- data.frame(
          setting = "cGCV",
          family = fam_name,
          optimizer = opt,
          rank = r,
          k = k,
          n = length(dat$y),
          status = "ok",
          message = "",
          lambda_geo = exp(mean(log(pmax(fit$lambda, 1e-12)))),
          lambda_1 = fit$lambda[1],
          lambda_2 = fit$lambda[2],
          lambda_3 = fit$lambda[3],
          npar_tt = fit$npar_tt,
          edf = fit$edf %||% NA_real_,
          deviance = fit$deviance,
          rmse_eta_test = rmse(eta_te, truth_eta),
          cor_eta_test = suppressWarnings(cor(eta_te, truth_eta)),
          time_s = elapsed,
          n_criterion_evals = fit$n_criterion_evals %||% NA_real_,
          n_outer = fit$n_outer %||% NA_real_,
          n_opt_iter = fit$n_opt_iter %||% NA_real_,
          backend = fit$backend,
          converged = isTRUE(fit$converged),
          stringsAsFactors = FALSE
        )
        message(sprintf(
          "  %-9s %-7s r=%d | RMSE_te=%.4f | λ_geo=%.3g | GCV_evals=%s | %.2fs",
          fam_name, opt, r, rmse(eta_te, truth_eta),
          exp(mean(log(pmax(fit$lambda, 1e-12)))),
          as.character(fit$n_criterion_evals %||% NA),
          elapsed
        ))
      }
    }
  }

  tab <- do.call(rbind, rows)
  print(tab[, c("family", "optimizer", "rank", "status", "rmse_eta_test",
                "lambda_geo", "edf", "n_criterion_evals", "time_s")],
        row.names = FALSE, digits = 4)
  .write_bench_csv(tab, "benchmark_optimizers_cgcv.csv", out_dir)
  invisible(tab)
})
