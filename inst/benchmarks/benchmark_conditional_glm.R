# Conditional GLM solvers gate: PIRLS-ALS vs Damped-Newton-ALS vs LBFGS-ALS vs LBFGS
#
# Scientific goal (docs/OPTIMIZER_CONDITIONAL_GLM.md):
#   Separate the PIRLS path from alternating / core-wise TT optimization.
#   Fixed λ only; no cGCV; no GD in this gate.
#
# Usage:
#   Rscript inst/benchmarks/benchmark_conditional_glm.R --quick
#   Rscript inst/benchmarks/benchmark_conditional_glm.R --full
#
# Results -> inst/benchmarks/results/conditional_glm/

args <- commandArgs(trailingOnly = TRUE)
mode <- if (any(grepl("^--full$", args))) "full" else "quick"
family_arg <- if (any(grepl("^--poisson-only$", args))) {
  "poisson"
} else if (any(grepl("^--bernoulli-only$", args))) {
  "bernoulli"
} else {
  "both"
}

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

  out_dir <- file.path(.ttps_bench_outdir(), "conditional_glm")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== Conditional GLM gate (", mode, ", ", family_arg, ") ===")
  message("Results -> ", out_dir)

  if (identical(mode, "quick")) {
    n <- 400L; n_te <- 1200L; k <- 6L
    rank_grid <- c(2L, 3L)
    n_init <- 3L
  } else {
    n <- 800L; n_te <- 2500L; k <- 6L
    rank_grid <- c(2L, 3L)
    n_init <- 10L
  }

  ETA_CATASTROPHIC <- 20
  opts <- c("PIRLS-ALS", "Damped-Newton-ALS", "LBFGS-ALS", "LBFGS")
  lambda_fixed <- 1

  .ctrl <- function(opt) {
    if (identical(opt, "PIRLS-ALS")) {
      tt_control(
        backend = "R", pirls_maxit = 40L, als_sweeps_per_pirls = 10L,
        pirls_step_halving = TRUE, compute_edf = FALSE
      )
    } else if (identical(opt, "Damped-Newton-ALS")) {
      tt_control(
        backend = "R", dn_max_sweeps = 40L, dn_armijo_c = 1e-4,
        dn_step_factor = 0.5, dn_ridge = 0, compute_edf = FALSE, tol = 1e-8
      )
    } else if (identical(opt, "LBFGS-ALS")) {
      tt_control(
        backend = "R", block_lbfgs_sweeps = 40L, block_lbfgs_maxit = 50L,
        compute_edf = FALSE, tol = 1e-8
      )
    } else {
      tt_control(backend = "R", lbfgs_maxit = 300L, compute_edf = FALSE)
    }
  }

  .is_cat <- function(obj, max_ae) {
    (!is.finite(obj)) || (is.finite(max_ae) && max_ae > ETA_CATASTROPHIC)
  }

  .bernoulli_logloss <- function(y, p) {
    p <- pmin(pmax(as.numeric(p), 1e-12), 1 - 1e-12)
    -mean(y * log(p) + (1 - y) * log(1 - p))
  }

  .bernoulli_deviance <- function(y, p) {
    p <- pmin(pmax(as.numeric(p), 1e-12), 1 - 1e-12)
    2 * sum(
      ifelse(y > 0, y * log(y / p), 0) +
        ifelse(y < 1, (1 - y) * log((1 - y) / (1 - p)), 0)
    )
  }

  .poisson_deviance <- function(y, mu) {
    mu <- pmax(as.numeric(mu), 1e-12)
    2 * sum(ifelse(y > 0, y * log(y / mu) - (y - mu), mu))
  }

  .min_weight_from_eta <- function(eta, family_name) {
    if (identical(family_name, "bernoulli")) {
      mu <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12)
      return(min(mu * (1 - mu)))
    }
    if (identical(family_name, "poisson")) {
      return(min(pmax(exp(pmin(pmax(eta, -20), 20)), 1e-12)))
    }
    NA_real_
  }

  .metrics_bern <- function(fit, tag, r, lam, seed, dat, te) {
    o <- tt_objective(fit, dat$X, dat$y)
    eta_te <- predict(fit, te$X, type = "link")
    p_tr <- predict(fit, type = "response")
    p_te <- predict(fit, te$X, type = "response")
    max_ae_tr <- max(abs(fit$linear.predictors))
    max_ae_te <- max(abs(eta_te))
    # Shared hard test labels across methods (seed from DGP shift only).
    set.seed(99000L + as.integer(round(1000 * abs(dat$shift))))
    y_te <- rbinom(length(te$truth_p), 1L, te$truth_p)
    data.frame(
      family = "bernoulli", method = tag, rank = r, lambda = lam, seed = seed,
      train_objective = o$value, train_nll = o$nll_or_sse, train_penalty = o$penalty,
      train_deviance = fit$deviance,
      rmse_eta_train = rmse(fit$linear.predictors, dat$truth_eta),
      rmse_eta_test = rmse(eta_te, te$truth_eta),
      logloss_test = .bernoulli_logloss(y_te, p_te),
      deviance_test = .bernoulli_deviance(y_te, p_te),
      max_abs_eta_train = max_ae_tr,
      max_abs_eta_test = max_ae_te,
      min_p = min(p_tr), max_p = max(p_tr),
      min_weight = .min_weight_from_eta(fit$linear.predictors, "bernoulli"),
      time_s = fit$timing,
      n_sweeps = fit$n_sweeps %||% fit$convergence$n_als_sweeps %||% NA_real_,
      n_opt_iter = fit$n_opt_iter %||% NA_real_,
      n_backtrack = fit$convergence$n_backtrack %||%
        fit$convergence$n_step_halvings %||% NA_real_,
      median_alpha = fit$convergence$median_alpha %||% NA_real_,
      converged = isTRUE(fit$convergence$overall %||% fit$converged),
      reason = as.character(fit$convergence$reason %||% NA_character_),
      catastrophic = .is_cat(o$value, max_ae_tr),
      stringsAsFactors = FALSE
    )
  }

  .metrics_pois <- function(fit, tag, r, lam, seed, dat, te) {
    o <- tt_objective(fit, dat$X, dat$y)
    eta_te <- predict(fit, te$X, type = "link")
    mu_te <- predict(fit, te$X, type = "response")
    max_ae_tr <- max(abs(fit$linear.predictors))
    max_ae_te <- max(abs(eta_te))
    set.seed(88000L + as.integer(round(1000 * abs(dat$shift))))
    y_te <- rpois(length(te$truth_mu), te$truth_mu)
    data.frame(
      family = "poisson", method = tag, rank = r, lambda = lam, seed = seed,
      train_objective = o$value, train_nll = o$nll_or_sse, train_penalty = o$penalty,
      train_deviance = fit$deviance,
      rmse_eta_train = rmse(fit$linear.predictors, dat$truth_eta),
      rmse_eta_test = rmse(eta_te, te$truth_eta),
      logloss_test = NA_real_,
      deviance_test = .poisson_deviance(y_te, mu_te),
      max_abs_eta_train = max_ae_tr,
      max_abs_eta_test = max_ae_te,
      min_p = NA_real_, max_p = NA_real_,
      min_weight = .min_weight_from_eta(fit$linear.predictors, "poisson"),
      time_s = fit$timing,
      n_sweeps = fit$n_sweeps %||% fit$convergence$n_als_sweeps %||% NA_real_,
      n_opt_iter = fit$n_opt_iter %||% NA_real_,
      n_backtrack = fit$convergence$n_backtrack %||%
        fit$convergence$n_step_halvings %||% NA_real_,
      median_alpha = fit$convergence$median_alpha %||% NA_real_,
      converged = isTRUE(fit$convergence$overall %||% fit$converged),
      reason = as.character(fit$convergence$reason %||% NA_character_),
      catastrophic = .is_cat(o$value, max_ae_tr),
      stringsAsFactors = FALSE
    )
  }

  .fit_bern <- function(init, opt, r, lam, seed, dat, te) {
    fit <- ttpspline(
      dat$y, dat$X, family = binomial(), rank = r, k = k, lambda = lam,
      optimizer = opt, init = init, control = .ctrl(opt)
    )
    .metrics_bern(fit, opt, r, lam, seed, dat, te)
  }

  .fit_pois <- function(init, opt, r, lam, seed, dat, te) {
    fit <- ttpspline(
      dat$y, dat$X, family = poisson(), rank = r, k = k, lambda = lam,
      optimizer = opt, init = init, control = .ctrl(opt)
    )
    .metrics_pois(fit, opt, r, lam, seed, dat, te)
  }

  .summ <- function(sub, col) {
    x <- sub[[col]]
    x <- x[is.finite(x)]
    if (!length(x)) return(c(med = NA_real_, q1 = NA_real_, q3 = NA_real_))
    c(
      med = stats::median(x),
      q1 = as.numeric(stats::quantile(x, 0.25)),
      q3 = as.numeric(stats::quantile(x, 0.75))
    )
  }

  .summarize <- function(multi, ranks, methods) {
    rows <- list()
    for (r in ranks) {
      for (opt in methods) {
        sub <- multi[multi$rank == r & multi$method == opt, , drop = FALSE]
        if (!nrow(sub)) next
        s_rmse <- .summ(sub, "rmse_eta_test")
        s_ll <- .summ(sub, "logloss_test")
        s_eta <- .summ(sub, "max_abs_eta_train")
        s_t <- .summ(sub, "time_s")
        rows[[length(rows) + 1L]] <- data.frame(
          family = sub$family[[1]],
          method = opt, rank = r, n = nrow(sub),
          median_rmse = unname(s_rmse["med"]),
          median_logloss = unname(s_ll["med"]),
          median_max_abs_eta = unname(s_eta["med"]),
          median_runtime = unname(s_t["med"]),
          failure_rate = mean(sub$catastrophic),
          stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows)
  }

  run_family <- function(family_name) {
    message("\n######## ", toupper(family_name), " ########")
    if (identical(family_name, "bernoulli")) {
      set.seed(3)
      dat <- simulate_bernoulli(n = n, d = 3L, seed = 3L)
      te <- holdout_bernoulli(n_te = n_te, d = 3L, shift = dat$shift, seed = 99L)
      .fit <- .fit_bern
    } else {
      set.seed(2)
      dat <- simulate_poisson(n = n, d = 3L, seed = 2L)
      te <- holdout_poisson(n_te = n_te, d = 3L, shift = dat$shift, seed = 99L)
      .fit <- .fit_pois
    }

    message("\n--- Exp1: shared init, lambda=1 ---")
    exp1_rows <- list()
    for (r in rank_grid) {
      init <- tt_initialize(dat$X, rank = r, k = k, seed = 123L, sd = 0.05)
      for (opt in opts) {
        m <- .fit(init, opt, r, lambda_fixed, seed = NA_integer_, dat, te)
        exp1_rows[[length(exp1_rows) + 1L]] <- m
        message(sprintf(
          "  r=%d %-20s | L=%.4g | RMSE_te=%.3f | max|eta|=%.2f | t=%.2fs | cat=%s",
          r, opt, m$train_objective, m$rmse_eta_test, m$max_abs_eta_train,
          m$time_s, m$catastrophic
        ))
      }
    }
    exp1 <- do.call(rbind, exp1_rows)
    write.csv(exp1, file.path(out_dir, paste0(family_name, "_exp1_shared_init.csv")),
              row.names = FALSE)

    message(sprintf("\n--- Exp2: multi-init (%d seeds) ---", n_init))
    multi_rows <- list()
    for (r in rank_grid) {
      for (s in seq_len(n_init)) {
        init <- tt_initialize(dat$X, rank = r, k = k, seed = 2000L + s, sd = 0.05)
        for (opt in opts) {
          multi_rows[[length(multi_rows) + 1L]] <-
            .fit(init, opt, r, lambda_fixed, seed = s, dat, te)
        }
      }
    }
    multi <- do.call(rbind, multi_rows)
    write.csv(multi, file.path(out_dir, paste0(family_name, "_exp2_multiinit.csv")),
              row.names = FALSE)
    summ <- .summarize(multi, rank_grid, opts)
    write.csv(summ, file.path(out_dir, paste0(family_name, "_exp2_summary.csv")),
              row.names = FALSE)
    message("\n--- Summary ---")
    print(summ, row.names = FALSE)
    list(exp1 = exp1, multi = multi, summ = summ)
  }

  results <- list()
  if (family_arg %in% c("both", "bernoulli")) {
    results$bernoulli <- run_family("bernoulli")
  }
  if (family_arg %in% c("both", "poisson")) {
    results$poisson <- run_family("poisson")
  }

  saveRDS(
    list(mode = mode, n = n, n_init = n_init, ETA_CATASTROPHIC = ETA_CATASTROPHIC,
         opts = opts, results = results, timestamp = as.character(Sys.time())),
    file.path(out_dir, "gate_snapshot.rds")
  )
  writeLines(c(
    paste0("mode=", mode),
    paste0("n=", n),
    paste0("n_te=", n_te),
    paste0("n_init=", n_init),
    paste0("lambda=", lambda_fixed),
    paste0("opts=", paste(opts, collapse = ",")),
    paste0("ETA_CATASTROPHIC=", ETA_CATASTROPHIC),
    paste0("timestamp=", Sys.time())
  ), file.path(out_dir, "README.txt"))

  message("\n=== Conditional GLM gate complete ===")
})
