# Optimizer comparison — FIXED lambda (fair init / bases / data)
#
# Compares ALS vs LBFGS (Adam stub recorded as unavailable).
# Metrics: eta/mu prediction quality, deviance/RSS, wall time, iterations.
# Do NOT compare TT cores (gauge non-identifiability).
#
#   Rscript inst/benchmarks/benchmark_optimizers_fixed.R
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

  message("=== TTPsplines benchmark: optimizers (fixed lambda) ===")

  # Slightly lighter defaults so the full grid finishes in a few minutes
  n_tr <- as.integer(Sys.getenv("TTPSPLINES_BENCH_N", "1000"))
  k <- as.integer(Sys.getenv("TTPSPLINES_BENCH_K", "6"))
  ranks <- as.integer(strsplit(Sys.getenv("TTPSPLINES_BENCH_RANKS", "2,3"), ",")[[1]])
  optimizers <- c("ALS", "LBFGS", "Adam")
  seed_data <- 21L
  seed_init <- 123L
  lambda_fixed <- 1

  sim_list <- list(
    gaussian = simulate_gaussian(n = n_tr, d = 3L, sigma = 0.3, seed = seed_data),
    poisson = simulate_poisson(n = n_tr, d = 3L, seed = seed_data + 1L),
    bernoulli = simulate_bernoulli(n = n_tr, d = 3L, seed = seed_data + 2L)
  )
  hold_list <- list(
    gaussian = holdout_gaussian(n_te = 3000L, d = 3L, seed = 211L),
    poisson = holdout_poisson(n_te = 3000L, d = 3L, shift = sim_list$poisson$shift, seed = 212L),
    bernoulli = holdout_bernoulli(n_te = 3000L, d = 3L, shift = sim_list$bernoulli$shift, seed = 213L)
  )
  fam_obj <- list(
    gaussian = gaussian(),
    poisson = poisson(),
    bernoulli = binomial()
  )

  rows <- list()

  for (fam_name in names(sim_list)) {
    dat <- sim_list[[fam_name]]
    te <- hold_list[[fam_name]]
    for (r in ranks) {
      init <- tt_initialize(dat$X, rank = r, k = k, seed = seed_init, sd = 0.1)
      for (opt in optimizers) {
        ctrl <- .default_control(
          backend = "R",
          max_sweeps = 20L,
          pirls_maxit = 15L,
          als_sweeps_per_pirls = 3L,
          lbfgs_maxit = 200L,
          seed = seed_init,
          init_sd = 0.1,
          compute_edf = TRUE,
          damping = TRUE
        )
        t0 <- proc.time()[["elapsed"]]
        fit <- tryCatch(
          ttpspline(
            dat$y, dat$X,
            family = fam_obj[[fam_name]],
            rank = r,
            k = k,
            lambda = lambda_fixed,
            optimizer = opt,
            init = init,
            control = ctrl
          ),
          error = function(e) e
        )
        elapsed <- proc.time()[["elapsed"]] - t0

        if (inherits(fit, "error")) {
          rows[[length(rows) + 1L]] <- data.frame(
            setting = "fixed_lambda",
            family = fam_name,
            optimizer = opt,
            rank = r,
            k = k,
            n = length(dat$y),
            d = 3L,
            lambda_mode = "fixed",
            status = "error",
            message = conditionMessage(fit),
            npar_tt = NA_real_,
            npar_tt_intrinsic = NA_real_,
            compression = NA_real_,
            edf = NA_real_,
            deviance = NA_real_,
            rmse_eta_train = NA_real_,
            rmse_eta_test = NA_real_,
            rmse_mu_test = NA_real_,
            cor_eta_test = NA_real_,
            time_s = elapsed,
            n_sweeps = NA_real_,
            n_pirls = NA_real_,
            n_opt_iter = NA_real_,
            n_criterion_evals = NA_real_,
            backend = NA_character_,
            converged = FALSE,
            stringsAsFactors = FALSE
          )
          message(sprintf(
            "  %-9s %-7s r=%d | ERROR: %s",
            fam_name, opt, r, conditionMessage(fit)
          ))
          next
        }

        eta_tr <- fit$linear.predictors
        eta_te <- predict(fit, te$X, type = "link")
        mu_te <- predict(fit, te$X, type = "response")

        if (identical(fam_name, "gaussian")) {
          rmse_tr <- rmse(eta_tr, dat$truth)
          rmse_te <- rmse(eta_te, te$truth)
          rmse_mu <- rmse_te
          cor_te <- suppressWarnings(cor(eta_te, te$truth))
        } else if (identical(fam_name, "poisson")) {
          rmse_tr <- rmse(eta_tr, dat$truth_eta)
          rmse_te <- rmse(eta_te, te$truth_eta)
          rmse_mu <- rmse(mu_te, te$truth_mu)
          cor_te <- suppressWarnings(cor(eta_te, te$truth_eta))
        } else {
          rmse_tr <- rmse(eta_tr, dat$truth_eta)
          rmse_te <- rmse(eta_te, te$truth_eta)
          rmse_mu <- rmse(mu_te, te$truth_p)
          cor_te <- suppressWarnings(cor(eta_te, te$truth_eta))
        }

        rows[[length(rows) + 1L]] <- data.frame(
          setting = "fixed_lambda",
          family = fam_name,
          optimizer = opt,
          rank = r,
          k = k,
          n = length(dat$y),
          d = 3L,
          lambda_mode = "fixed",
          status = "ok",
          message = "",
          npar_tt = fit$npar_tt,
          npar_tt_intrinsic = fit$npar_tt_intrinsic %||% NA_real_,
          compression = fit$compression_ratio,
          edf = fit$edf %||% NA_real_,
          deviance = fit$deviance,
          rmse_eta_train = rmse_tr,
          rmse_eta_test = rmse_te,
          rmse_mu_test = rmse_mu,
          cor_eta_test = cor_te,
          time_s = elapsed,
          n_sweeps = fit$n_sweeps %||% NA_real_,
          n_pirls = fit$n_pirls %||% NA_real_,
          n_opt_iter = fit$n_opt_iter %||% NA_real_,
          n_criterion_evals = fit$n_criterion_evals %||% NA_real_,
          backend = fit$backend,
          converged = isTRUE(fit$converged),
          stringsAsFactors = FALSE
        )
        message(sprintf(
          "  %-9s %-7s r=%d | RMSE_te=%.4f | dev=%.3g | %.2fs | backend=%s",
          fam_name, opt, r, rmse_te, fit$deviance, elapsed, fit$backend
        ))
      }
    }
  }

  tab <- do.call(rbind, rows)
  print(tab[, c("family", "optimizer", "rank", "status", "rmse_eta_test",
                "deviance", "edf", "time_s", "converged")],
        row.names = FALSE, digits = 4)
  .write_bench_csv(tab, "benchmark_optimizers_fixed.csv", out_dir)

  # ALS vs LBFGS eta agreement on training (same init) for gaussian r=2 if both ok
  ok_g <- subset(tab, family == "gaussian" & rank == ranks[1] & status == "ok")
  if (nrow(ok_g) >= 2L) {
    fig <- file.path(out_dir, "benchmark_optimizers_fixed.png")
    png(fig, width = 1000, height = 420)
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
    on.exit({ par(op); dev.off() }, add = TRUE)
    for (fam in c("gaussian", "poisson")) {
      sub <- subset(tab, family == fam & status == "ok" & optimizer %in% c("ALS", "LBFGS"))
      if (!nrow(sub)) next
      plot(sub$rank, sub$rmse_eta_test, type = "n",
           xlab = "rank", ylab = "RMSE eta (test)",
           main = paste0(fam, " — fixed λ=", lambda_fixed))
      for (opt in unique(sub$optimizer)) {
        s2 <- sub[sub$optimizer == opt, ]
        lines(s2$rank, s2$rmse_eta_test, type = "b", pch = 19,
              col = if (opt == "ALS") "steelblue" else "darkorange")
      }
      legend("topright", legend = unique(sub$optimizer),
             col = c("steelblue", "darkorange")[seq_along(unique(sub$optimizer))],
             lty = 1, pch = 19, bty = "n")
    }
    message("Wrote ", fig)
  }

  invisible(tab)
})
