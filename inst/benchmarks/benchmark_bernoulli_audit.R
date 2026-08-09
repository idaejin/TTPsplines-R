# Bernoulli / PIRLS audit (fixed λ) — diagnose ALS vs L-BFGS gap
#
# Usage:
#   Rscript inst/benchmarks/benchmark_bernoulli_audit.R --quick
#   Rscript inst/benchmarks/benchmark_bernoulli_audit.R --full
#
# Does NOT change the main fitting algorithm. Diagnostics only.

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

  out_dir <- file.path(.ttps_bench_outdir(), "bernoulli_audit")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== Bernoulli PIRLS audit (", mode, ") ===")
  message("Results -> ", out_dir)

  # ---- config ----
  if (identical(mode, "quick")) {
    n <- 500L; n_te <- 1500L; k <- 6L; rank_grid <- c(2L, 3L)
    n_init <- 5L; lambdas <- c(0.1, 1, 10, 100)
    als_sweep_grid <- c(1L, 2L, 5L, 10L, 20L)
    pirls_grid <- c(10L, 25L, 50L)
  } else {
    n <- 800L; n_te <- 2500L; k <- 6L; rank_grid <- c(2L, 3L)
    n_init <- 20L; lambdas <- c(0.1, 1, 10, 100)
    als_sweep_grid <- c(1L, 2L, 5L, 10L, 20L, 40L)
    pirls_grid <- c(10L, 25L, 50L, 100L)
  }

  set.seed(3)
  dat <- simulate_bernoulli(n = n, d = 3L, seed = 3L)
  te <- holdout_bernoulli(n_te = n_te, d = 3L, shift = dat$shift, seed = 99L)
  y <- dat$y; X <- dat$X

  .obj <- function(fit) tt_objective(fit, X, y)
  .metrics <- function(fit, tag, r, lam, seed = NA_integer_) {
    o <- .obj(fit)
    eta_te <- predict(fit, te$X, type = "link")
    p_tr <- predict(fit, type = "response")
    data.frame(
      tag = tag, rank = r, lambda = lam, seed = seed,
      objective = o$value, nll = o$nll_or_sse, penalty = o$penalty,
      deviance = fit$deviance,
      rmse_eta_train = rmse(fit$linear.predictors, dat$truth_eta),
      rmse_eta_test = rmse(eta_te, te$truth_eta),
      max_abs_eta = max(abs(fit$linear.predictors)),
      min_p = min(p_tr), max_p = max(p_tr),
      prop_p_lo = mean(p_tr < 1e-6), prop_p_hi = mean(p_tr > 1 - 1e-6),
      time_s = fit$timing, n_pirls = fit$n_pirls %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # ==================================================================
  # Q1 — Objective parity (same cores / intercept / λ / bases)
  # ==================================================================
  message("\n--- Q1: objective parity ---")
  parity_rows <- list()
  for (r in rank_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    # Fit briefly only to get a nontrivial point, then evaluate both paths
    fit_a <- ttpspline(
      y, X, family = binomial(), rank = r, k = k, lambda = 1,
      optimizer = "ALS", init = init,
      control = tt_control(backend = "R", pirls_maxit = 3L,
                           als_sweeps_per_pirls = 2L, compute_edf = FALSE,
                           damping = FALSE, seed = 123L)
    )
    # Evaluate objective via public helper
    o1 <- tt_objective(fit_a, X, y)
    # Rebuild via LBFGS internal objective on same packed cores
    basis <- eval_marginal_bases(X, fit_a$knots, fit_a$degree)
    d <- length(fit_a$cores)
    p <- ncol(basis[[1]])
    ranks_full <- fit_a$rank
    pens <- lapply(seq_len(d), function(kk)
      core_penalty(ranks_full[kk], p, ranks_full[kk + 1L], 2L))
    th <- TTPsplines:::.tt_pack_cores(fit_a$cores)
    o2 <- TTPsplines:::.tt_glm_objective(
      th, y, fit_a$intercept, basis, fit_a$cores, pens, fit_a$lambda,
      binomial()
    )
    eta_a <- fit_a$linear.predictors
    eta_b <- o2$eta
    rel_obj <- abs(o1$value - o2$value) / max(1, abs(o1$value))
    rel_eta <- sqrt(mean((eta_a - eta_b)^2)) / max(1, sqrt(mean(eta_a^2)))
    parity_rows[[length(parity_rows) + 1L]] <- data.frame(
      rank = r,
      obj_tt_objective = o1$value,
      obj_lbfgs_internal = o2$value,
      rel_obj = rel_obj,
      rmse_eta = sqrt(mean((eta_a - eta_b)^2)),
      rel_eta = rel_eta,
      nll = o1$nll_or_sse,
      penalty = o1$penalty,
      stringsAsFactors = FALSE
    )
    message(sprintf(
      "  r=%d | rel_obj=%.3e | rmse_eta=%.3e | PASS=%s",
      r, rel_obj, sqrt(mean((eta_a - eta_b)^2)),
      rel_obj < 1e-8 && sqrt(mean((eta_a - eta_b)^2)) < 1e-8
    ))
  }
  parity <- do.call(rbind, parity_rows)
  write.csv(parity, file.path(out_dir, "q1_objective_parity.csv"), row.names = FALSE)
  if (any(parity$rel_obj >= 1e-8)) {
    message("Q1 FAILED: objective mismatch — aborting further interpretation.")
    quit(save = "no", status = 1)
  }
  message("Q1 PASS: same Bernoulli penalized objective.")

  # ==================================================================
  # Instrumented PIRLS (local copy — does not modify package algorithm)
  # ==================================================================
  .glm_working_diag <- function(y, eta) {
    mu <- plogis(eta)
    mu_c <- pmin(pmax(mu, 1e-5), 1 - 1e-5)
    var <- mu_c * (1 - mu_c)
    w <- pmax(var, 1e-4)
    z <- eta + (y - mu_c) / var # NOTE: uses var, not w (package behavior)
    list(mu = mu_c, var = var, weight = w, z = z,
         min_p = min(mu_c), max_p = max(mu_c),
         min_w = min(w), med_w = stats::median(w), max_w = max(w),
         min_z = min(z), max_z = max(z))
  }

  .bern_obj_raw <- function(y, eta, cores, intercept, basis, pens, lambda) {
    mu <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12)
    nll <- -sum(y * log(mu) + (1 - y) * log(1 - mu))
    pen <- TTPsplines:::.tt_penalty_value_grad(cores, pens, lambda)$value
    list(value = nll + pen, nll = nll, penalty = pen, mu = mu)
  }

  .traced_pirls <- function(init, r, lambda, maxit, als_sweeps,
                            step_halving = FALSE, tol = 1e-8) {
    rank_chain <- tt_rank(r, d = 3L)
    bs <- build_marginal_bases(X, k = k, degree = 3)
    basis <- bs$basis
    p <- bs$k
    cores <- init
    pens <- lapply(seq_len(3), function(kk)
      core_penalty(rank_chain[kk], p, rank_chain[kk + 1L], 2L))
    lam <- rep(lambda, 3)
    intercept <- init_intercept(binomial(), y)
    eta <- intercept + tt_contraction(cores, basis)
    hist <- list()
    n_reject <- 0L
    alphas <- numeric()
    cond_hist <- list()

    for (it in seq_len(maxit)) {
      work <- .glm_working_diag(y, eta)
      w <- work$weight; z <- work$z
      obj_before <- .bern_obj_raw(y, eta, cores, intercept, basis, pens, lam)

      # inner ALS on working LS
      for (sw in seq_len(als_sweeps)) {
        zc <- z - intercept
        for (kk in 1:3) {
          L <- left_interfaces(cores, basis)
          R <- right_interfaces(cores, basis)
          Xk <- tt_design_core(L[[kk]], R[[kk]], basis[[kk]])
          sws <- sqrt(pmax(w, 0))
          Xw <- Xk * sws
          yw <- zc * sws
          S <- crossprod(Xw)
          b <- as.numeric(crossprod(Xw, yw))
          A <- S + lam[kk] * pens[[kk]]
          # conditioning diagnostics (first core, first sweep only to save time)
          if (sw == 1L && kk == 1L) {
            ev <- tryCatch(eigen(A, symmetric = TRUE, only.values = TRUE)$values,
                           error = function(e) NA_real_)
            cond <- if (all(is.finite(ev))) max(ev) / max(min(ev), 1e-16) else NA_real_
            cond_hist[[length(cond_hist) + 1L]] <- data.frame(
              pirls = it, core = kk, cond = cond,
              min_ev = if (all(is.finite(ev))) min(ev) else NA_real_,
              max_ev = if (all(is.finite(ev))) max(ev) else NA_real_
            )
          }
          g <- tryCatch(
            solve_spd_ridge(A, b),
            error = function(e) rep(0, nrow(A))
          )
          cores[[kk]] <- array(g, c(rank_chain[kk], p, rank_chain[kk + 1L]))
        }
        f <- tt_contraction(cores, basis)
        intercept <- sum(w * (z - f)) / max(sum(w), 1e-12)
      }

      eta_new <- intercept + tt_contraction(cores, basis)
      obj_full <- .bern_obj_raw(y, eta_new, cores, intercept, basis, pens, lam)

      alpha <- 1
      if (step_halving) {
        eta_old <- eta
        cores_new <- cores
        intercept_new <- intercept
        # step-halving on eta only for objective check; if needed, shrink toward old
        # Full model step: interpolate cores+intercept in eta space by damping cores
        for (a in c(1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64)) {
          # linear blend of eta, then least-commitment: accept if obj decreases
          eta_try <- eta_old + a * (eta_new - eta_old)
          # Approximate: scale core update by a from previous cores snapshot
          # We don't have previous cores easily after overwrite — store before ALS
          obj_try <- .bern_obj_raw(y, eta_try, cores_new, intercept_new, basis, pens, lam)
          # For true core blend we'd need cores_old; use eta-based NLL + current penalty as proxy
          # Better: keep cores_old
          if (obj_try$value <= obj_before$value * (1 + 1e-12) || a < 1/64 + 1e-12) {
            alpha <- a
            if (a < 1) n_reject <- n_reject + 1L
            # If alpha < 1, shrink cores toward previous by blending eta via intercept/cores
            if (a < 1) {
              # Re-run is expensive; accept eta_try with cores_new scaled update:
              # cores := cores_old + a (cores_new - cores_old) — need cores_old
            }
            break
          }
        }
        alphas <- c(alphas, alpha)
        # Simplified accept: if alpha==1 use eta_new else keep old (reject) unless a tiny
        if (alpha < 1 && obj_full$value > obj_before$value) {
          # reject full step (keep previous) — matches cautious PIRLS
          hist[[length(hist) + 1L]] <- data.frame(
            pirls = it, objective = obj_before$value, nll = obj_before$nll,
            penalty = obj_before$penalty,
            deviance = glm_deviance(binomial(), y, obj_before$mu),
            max_abs_eta = max(abs(eta)),
            d_eta = 0,
            alpha = alpha, rejected = TRUE,
            min_p = work$min_p, max_p = work$max_p,
            min_w = work$min_w, med_w = work$med_w, max_w = work$max_w,
            min_z = work$min_z, max_z = work$max_z
          )
          next
        }
      } else {
        alphas <- c(alphas, 1)
      }

      # Soft package-like guard (record only; optional)
      d_eta <- sqrt(mean((eta_new - eta)^2))
      hist[[length(hist) + 1L]] <- data.frame(
        pirls = it,
        objective = obj_full$value,
        nll = obj_full$nll,
        penalty = obj_full$penalty,
        deviance = glm_deviance(binomial(), y, obj_full$mu),
        max_abs_eta = max(abs(eta_new)),
        d_eta = d_eta,
        alpha = alpha,
        rejected = FALSE,
        min_p = work$min_p, max_p = work$max_p,
        min_w = work$min_w, med_w = work$med_w, max_w = work$max_w,
        min_z = work$min_z, max_z = work$max_z,
        obj_increased = obj_full$value > obj_before$value + 1e-10
      )
      eta <- eta_new
      if (it > 2L && d_eta < tol) break
    }

    list(
      cores = cores, intercept = intercept, lambda = lam, ranks = rank_chain,
      eta = eta, mu = plogis(eta), history = do.call(rbind, hist),
      cond = if (length(cond_hist)) do.call(rbind, cond_hist) else NULL,
      n_reject = n_reject, alphas = alphas,
      basis = basis, knots = bs$knots, degree = bs$degree, k = p
    )
  }

  # ==================================================================
  # Q3–Q5: trajectories + separation (package ALS vs LBFGS + traced)
  # ==================================================================
  message("\n--- Q3–Q5: trajectories / separation (λ=1) ---")
  traj_rows <- list()
  sep_rows <- list()
  for (r in rank_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    fit_als <- ttpspline(
      y, X, family = binomial(), rank = r, k = k, lambda = 1,
      optimizer = "ALS", init = init,
      control = tt_control(backend = "R", pirls_maxit = 50L,
                           als_sweeps_per_pirls = 3L, tol = 1e-8,
                           compute_edf = FALSE, damping = TRUE, seed = 123)
    )
    fit_lbf <- ttpspline(
      y, X, family = binomial(), rank = r, k = k, lambda = 1,
      optimizer = "LBFGS", init = init,
      control = tt_control(backend = "R", lbfgs_maxit = 400L,
                           compute_edf = FALSE, seed = 123)
    )
    tr <- .traced_pirls(init, r, 1, maxit = 40L, als_sweeps = 3L, step_halving = FALSE)
    tr$history$rank <- r
    tr$history$method <- "traced_PIRLS"
    traj_rows[[length(traj_rows) + 1L]] <- tr$history
    if (!is.null(tr$cond)) {
      tr$cond$rank <- r
      write.csv(tr$cond, file.path(out_dir, sprintf("q6_cond_r%d.csv", r)),
                row.names = FALSE)
    }
    for (nm in c("ALS", "LBFGS")) {
      fit <- if (nm == "ALS") fit_als else fit_lbf
      o <- .obj(fit)
      p <- predict(fit, type = "response")
      sep_rows[[length(sep_rows) + 1L]] <- data.frame(
        rank = r, method = nm,
        objective = o$value, rmse_test = rmse(predict(fit, te$X), te$truth_eta),
        max_abs_eta = max(abs(fit$linear.predictors)),
        min_p = min(p), max_p = max(p),
        prop_lo = mean(p < 1e-6), prop_hi = mean(p > 1 - 1e-6),
        n_obj_increases = if (nm == "ALS") sum(tr$history$obj_increased) else NA_integer_,
        stringsAsFactors = FALSE
      )
    }
    message(sprintf(
      "  r=%d | ALS obj=%.4g RMSE=%.3f max|eta|=%.2f | LBFGS obj=%.4g RMSE=%.3f max|eta|=%.2f | PIRLS obj↑=%d/%d",
      r, .obj(fit_als)$value, rmse(predict(fit_als, te$X), te$truth_eta),
      max(abs(fit_als$linear.predictors)),
      .obj(fit_lbf)$value, rmse(predict(fit_lbf, te$X), te$truth_eta),
      max(abs(fit_lbf$linear.predictors)),
      sum(tr$history$obj_increased), nrow(tr$history)
    ))

    # plots
    png(file.path(out_dir, sprintf("q3_traj_r%d.png", r)), width = 900, height = 700)
    op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    h <- tr$history
    plot(h$pirls, h$objective, type = "b", main = paste0("obj r=", r),
         xlab = "PIRLS", ylab = "penalized NLL")
    plot(h$pirls, h$deviance, type = "b", main = "deviance", xlab = "PIRLS", ylab = "dev")
    plot(h$pirls, h$max_abs_eta, type = "b", main = "max|eta|", xlab = "PIRLS", ylab = "|eta|")
    plot(h$pirls, h$med_w, type = "b", main = "median W", xlab = "PIRLS", ylab = "W")
    par(op); dev.off()
  }
  write.csv(do.call(rbind, traj_rows), file.path(out_dir, "q3_trajectories.csv"),
            row.names = FALSE)
  write.csv(do.call(rbind, sep_rows), file.path(out_dir, "q5_separation.csv"),
            row.names = FALSE)

  # ==================================================================
  # Q4 — step-halving experiment (traced)
  # ==================================================================
  message("\n--- Q4: step-halving ---")
  sh_rows <- list()
  for (r in rank_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    raw <- .traced_pirls(init, r, 1, 30L, 3L, step_halving = FALSE)
    # Proper step-halving with core blend
    rank_chain <- tt_rank(r, d = 3L)
    bs <- build_marginal_bases(X, k = k, degree = 3)
    basis <- bs$basis; p <- bs$k
    pens <- lapply(seq_len(3), function(kk)
      core_penalty(rank_chain[kk], p, rank_chain[kk + 1L], 2L))
    cores <- init
    lam <- rep(1, 3)
    intercept <- init_intercept(binomial(), y)
    eta <- intercept + tt_contraction(cores, basis)
    n_reject <- 0L; alphas <- c()
    hist_sh <- list()
    for (it in seq_len(30)) {
      work <- .glm_working_diag(y, eta)
      w <- work$weight; z <- work$z
      cores_old <- cores; intercept_old <- intercept; eta_old <- eta
      obj0 <- .bern_obj_raw(y, eta, cores, intercept, basis, pens, lam)
      for (sw in 1:3) {
        zc <- z - intercept
        for (kk in 1:3) {
          L <- left_interfaces(cores, basis); R <- right_interfaces(cores, basis)
          Xk <- tt_design_core(L[[kk]], R[[kk]], basis[[kk]])
          sws <- sqrt(pmax(w, 0)); Xw <- Xk * sws; yw <- zc * sws
          g <- solve_spd_ridge(crossprod(Xw) + lam[kk] * pens[[kk]],
                               as.numeric(crossprod(Xw, yw)))
          cores[[kk]] <- array(g, c(rank_chain[kk], p, rank_chain[kk + 1L]))
        }
        f <- tt_contraction(cores, basis)
        intercept <- sum(w * (z - f)) / max(sum(w), 1e-12)
      }
      cores_new <- cores; intercept_new <- intercept
      alpha_acc <- NA
      for (a in c(1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64)) {
        cores_try <- lapply(seq_along(cores_old), function(kk)
          cores_old[[kk]] + a * (cores_new[[kk]] - cores_old[[kk]]))
        intercept_try <- intercept_old + a * (intercept_new - intercept_old)
        eta_try <- intercept_try + tt_contraction(cores_try, basis)
        obj_try <- .bern_obj_raw(y, eta_try, cores_try, intercept_try, basis, pens, lam)
        if (obj_try$value <= obj0$value + 1e-10) {
          alpha_acc <- a
          cores <- cores_try; intercept <- intercept_try; eta <- eta_try
          if (a < 1) n_reject <- n_reject + 1L
          break
        }
      }
      if (is.na(alpha_acc)) {
        # keep old
        cores <- cores_old; intercept <- intercept_old; eta <- eta_old
        alpha_acc <- 0; n_reject <- n_reject + 1L
      }
      alphas <- c(alphas, alpha_acc)
      obj1 <- .bern_obj_raw(y, eta, cores, intercept, basis, pens, lam)
      hist_sh[[it]] <- data.frame(pirls = it, objective = obj1$value, alpha = alpha_acc)
    }
    sh_rows[[length(sh_rows) + 1L]] <- data.frame(
      rank = r,
      raw_final_obj = tail(raw$history$objective, 1),
      raw_obj_increases = sum(raw$history$obj_increased),
      sh_final_obj = tail(do.call(rbind, hist_sh)$objective, 1),
      sh_n_reject = n_reject,
      sh_median_alpha = stats::median(alphas),
      sh_min_alpha = min(alphas),
      stringsAsFactors = FALSE
    )
    message(sprintf(
      "  r=%d | raw obj↑=%d final=%.4g | step-halving rejects=%d med_alpha=%.3g final=%.4g",
      r, sum(raw$history$obj_increased), tail(raw$history$objective, 1),
      n_reject, stats::median(alphas), tail(do.call(rbind, hist_sh)$objective, 1)
    ))
  }
  write.csv(do.call(rbind, sh_rows), file.path(out_dir, "q4_step_halving.csv"),
            row.names = FALSE)

  # ==================================================================
  # Q7 — ALS sweeps per PIRLS
  # ==================================================================
  message("\n--- Q7: als_sweeps_per_pirls ---")
  sweep_rows <- list()
  for (r in rank_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    for (nsw in als_sweep_grid) {
      t0 <- proc.time()[["elapsed"]]
      fit <- ttpspline(
        y, X, family = binomial(), rank = r, k = k, lambda = 1,
        optimizer = "ALS", init = init,
        control = tt_control(backend = "R", pirls_maxit = 40L,
                             als_sweeps_per_pirls = nsw, tol = 1e-8,
                             compute_edf = FALSE, damping = TRUE)
      )
      m <- .metrics(fit, paste0("sweeps_", nsw), r, 1)
      m$als_sweeps_per_pirls <- nsw
      m$time_s <- proc.time()[["elapsed"]] - t0
      sweep_rows[[length(sweep_rows) + 1L]] <- m
      message(sprintf("  r=%d sweeps=%2d | obj=%.4g RMSE_te=%.3f max|eta|=%.2f",
                      r, nsw, m$objective, m$rmse_eta_test, m$max_abs_eta))
    }
    # LBFGS reference
    fit <- ttpspline(
      y, X, family = binomial(), rank = r, k = k, lambda = 1,
      optimizer = "LBFGS", init = init,
      control = tt_control(backend = "R", lbfgs_maxit = 400L, compute_edf = FALSE)
    )
    m <- .metrics(fit, "LBFGS", r, 1)
    m$als_sweeps_per_pirls <- NA
    sweep_rows[[length(sweep_rows) + 1L]] <- m
  }
  write.csv(do.call(rbind, sweep_rows), file.path(out_dir, "q7_inner_als.csv"),
            row.names = FALSE)

  # ==================================================================
  # Q8 — outer PIRLS maxit
  # ==================================================================
  message("\n--- Q8: pirls_maxit ---")
  outer_rows <- list()
  for (r in rank_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    for (pm in pirls_grid) {
      fit <- ttpspline(
        y, X, family = binomial(), rank = r, k = k, lambda = 1,
        optimizer = "ALS", init = init,
        control = tt_control(backend = "R", pirls_maxit = pm,
                             als_sweeps_per_pirls = 10L, tol = 1e-8,
                             compute_edf = FALSE, damping = TRUE)
      )
      m <- .metrics(fit, paste0("pirls_", pm), r, 1)
      m$pirls_maxit <- pm
      outer_rows[[length(outer_rows) + 1L]] <- m
    }
  }
  write.csv(do.call(rbind, outer_rows), file.path(out_dir, "q8_outer_pirls.csv"),
            row.names = FALSE)

  # ==================================================================
  # Q9 — cross warm-start (preserve intercept when refining with LBFGS)
  # ==================================================================
  message("\n--- Q9: cross warm-start ---")
  cross_rows <- list()
  for (r in rank_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    ctrl_als <- tt_control(backend = "R", pirls_maxit = 50L,
                           als_sweeps_per_pirls = 10L, compute_edf = FALSE,
                           damping = TRUE)
    ctrl_lbf <- tt_control(backend = "R", lbfgs_maxit = 400L, compute_edf = FALSE)
    A <- ttpspline(y, X, family = binomial(), rank = r, k = k, lambda = 1,
                   optimizer = "ALS", init = init, control = ctrl_als)
    # LBFGS from A with SAME intercept (public API would reset intercept)
    basis_A <- eval_marginal_bases(X, A$knots, A$degree)
    lam_spec <- parse_lambda_spec(1, d = 3)
    rank_chain <- tt_rank(r, d = 3L)
    A2_raw <- TTPsplines:::.tt_lbfgs_optimize_cores(
      y, basis_A, rank_chain, rep(1, 3), ctrl_lbf, 2L, A$cores,
      family = binomial(), intercept0 = A$intercept
    )
    A2_fit <- A
    A2_fit$cores <- A2_raw$cores
    A2_fit$intercept <- A2_raw$intercept
    A2_fit$linear.predictors <- A2_raw$intercept + tt_contraction(A2_raw$cores, basis_A)
    A2_fit$fitted.values <- plogis(A2_fit$linear.predictors)
    A2_fit$deviance <- glm_deviance(binomial(), y, A2_fit$fitted.values)

    B <- ttpspline(y, X, family = binomial(), rank = r, k = k, lambda = 1,
                   optimizer = "LBFGS", init = init, control = ctrl_lbf)
    B2 <- ttpspline(y, X, family = binomial(), rank = r, k = k, lambda = 1,
                    optimizer = "ALS", init = B$cores, control = ctrl_als)
    oa <- .obj(A); oa2 <- tt_objective(A2_fit, X, y)
    ob <- .obj(B); ob2 <- .obj(B2)
    cross_rows[[length(cross_rows) + 1L]] <- data.frame(
      rank = r,
      obj_A = oa$value, obj_A2 = oa2$value, delta_A = oa2$value - oa$value,
      rmse_A = rmse(predict(A, te$X), te$truth_eta),
      rmse_A2 = rmse(A2_fit$intercept +
                       tt_contraction(A2_fit$cores, eval_marginal_bases(te$X, A$knots, A$degree)),
                     te$truth_eta),
      obj_B = ob$value, obj_B2 = ob2$value, delta_B = ob2$value - ob$value,
      rmse_B = rmse(predict(B, te$X), te$truth_eta),
      rmse_B2 = rmse(predict(B2, te$X), te$truth_eta),
      max_eta_A = max(abs(A$linear.predictors)),
      max_eta_A2 = max(abs(A2_fit$linear.predictors)),
      max_eta_B = max(abs(B$linear.predictors)),
      max_eta_B2 = max(abs(B2$linear.predictors)),
      stringsAsFactors = FALSE
    )
    message(sprintf(
      "  r=%d | A→A2 Δobj=%.4g RMSE %.3f→%.3f | B→B2 Δobj=%.4g RMSE %.3f→%.3f",
      r, oa2$value - oa$value,
      rmse(predict(A, te$X), te$truth_eta),
      rmse(A2_fit$intercept +
             tt_contraction(A2_fit$cores, eval_marginal_bases(te$X, A$knots, A$degree)),
           te$truth_eta),
      ob2$value - ob$value,
      rmse(predict(B, te$X), te$truth_eta),
      rmse(predict(B2, te$X), te$truth_eta)
    ))
  }
  write.csv(do.call(rbind, cross_rows), file.path(out_dir, "q9_cross_warmstart.csv"),
            row.names = FALSE)

  # ==================================================================
  # Q10 — multi-init
  # ==================================================================
  message(sprintf("\n--- Q10: multi-init (%d seeds) ---", n_init))
  multi_rows <- list()
  for (r in rank_grid) {
    for (s in seq_len(n_init)) {
      init <- tt_initialize(X, rank = r, k = k, seed = 2000L + s, sd = 0.05)
      for (opt in c("ALS", "LBFGS")) {
        fit <- ttpspline(
          y, X, family = binomial(), rank = r, k = k, lambda = 1,
          optimizer = opt, init = init,
          control = if (opt == "ALS")
            tt_control(backend = "R", pirls_maxit = 40L, als_sweeps_per_pirls = 10L,
                       compute_edf = FALSE, damping = TRUE)
          else
            tt_control(backend = "R", lbfgs_maxit = 300L, compute_edf = FALSE)
        )
        multi_rows[[length(multi_rows) + 1L]] <- .metrics(fit, opt, r, 1, seed = s)
      }
    }
  }
  multi <- do.call(rbind, multi_rows)
  write.csv(multi, file.path(out_dir, "q10_multiinit.csv"), row.names = FALSE)
  for (r in rank_grid) {
    for (opt in c("ALS", "LBFGS")) {
      sub <- multi[multi$rank == r & multi$tag == opt, ]
      message(sprintf(
        "  r=%d %-5s | med obj=%.4g IQR=[%.4g,%.4g] | med RMSE=%.3f | med max|eta|=%.2f",
        r, opt, median(sub$objective),
        quantile(sub$objective, 0.25), quantile(sub$objective, 0.75),
        median(sub$rmse_eta_test), median(sub$max_abs_eta)
      ))
    }
  }

  # ==================================================================
  # Q11 — lambda sensitivity
  # ==================================================================
  message("\n--- Q11: lambda sensitivity ---")
  lam_rows <- list()
  for (r in rank_grid) {
    init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
    for (lam in lambdas) {
      for (opt in c("ALS", "LBFGS")) {
        fit <- ttpspline(
          y, X, family = binomial(), rank = r, k = k, lambda = lam,
          optimizer = opt, init = init,
          control = if (opt == "ALS")
            tt_control(backend = "R", pirls_maxit = 40L, als_sweeps_per_pirls = 10L,
                       compute_edf = FALSE, damping = TRUE)
          else
            tt_control(backend = "R", lbfgs_maxit = 300L, compute_edf = FALSE)
        )
        lam_rows[[length(lam_rows) + 1L]] <- .metrics(fit, opt, r, lam)
      }
    }
  }
  write.csv(do.call(rbind, lam_rows), file.path(out_dir, "q11_lambda.csv"),
            row.names = FALSE)

  # ==================================================================
  # Q12 — tiny dense comparison (d=2, k=5) if quick allows
  # ==================================================================
  message("\n--- Q12: small dense-ish check (d=2) ---")
  set.seed(7)
  n2 <- 300L
  X2 <- matrix(runif(n2 * 2), n2, 2)
  eta2 <- 1.2 * sin(2 * pi * X2[, 1]) * cos(2 * pi * X2[, 2])
  y2 <- rbinom(n2, 1, plogis(eta2))
  dense_rows <- list()
  for (r in c(2L, 4L)) {
    init <- tt_initialize(X2, rank = r, k = 5, seed = 1, sd = 0.05)
    for (opt in c("ALS", "LBFGS")) {
      fit <- ttpspline(
        y2, X2, family = binomial(), rank = r, k = 5, lambda = 1,
        optimizer = opt, init = init,
        control = if (opt == "ALS")
          tt_control(backend = "R", pirls_maxit = 40L, als_sweeps_per_pirls = 10L,
                     compute_edf = FALSE, damping = TRUE)
        else
          tt_control(backend = "R", lbfgs_maxit = 300L, compute_edf = FALSE)
      )
      o <- tt_objective(fit, X2, y2)
      dense_rows[[length(dense_rows) + 1L]] <- data.frame(
        rank = r, optimizer = opt, objective = o$value,
        max_abs_eta = max(abs(fit$linear.predictors)),
        rmse_eta = rmse(fit$linear.predictors, eta2),
        stringsAsFactors = FALSE
      )
    }
  }
  write.csv(do.call(rbind, dense_rows), file.path(out_dir, "q12_d2_check.csv"),
            row.names = FALSE)

  # ==================================================================
  # Q13 — conditional EDF sanity (final ALS fit)
  # ==================================================================
  message("\n--- Q13: conditional EDF vs lambda ---")
  edf_rows <- list()
  r <- 2L
  init <- tt_initialize(X, rank = r, k = k, seed = 123L, sd = 0.05)
  fit <- ttpspline(
    y, X, family = binomial(), rank = r, k = k, lambda = 1,
    optimizer = "ALS", init = init,
    control = tt_control(backend = "R", pirls_maxit = 20L, als_sweeps_per_pirls = 5L,
                         compute_edf = FALSE, damping = TRUE)
  )
  basis <- eval_marginal_bases(X, fit$knots, fit$degree)
  work <- .glm_working_diag(y, fit$linear.predictors)
  L <- left_interfaces(fit$cores, basis)
  R <- right_interfaces(fit$cores, basis)
  for (lam in c(0.1, 1, 10, 100)) {
    for (kk in 1:3) {
      Xk <- tt_design_core(L[[kk]], R[[kk]], basis[[kk]])
      sws <- sqrt(pmax(work$weight, 0))
      S <- crossprod(Xk * sws)
      Pk <- fit$penalties[[kk]]
      if (is.null(Pk)) {
        rk <- fit$rank
        Pk <- core_penalty(rk[kk], ncol(basis[[1]]), rk[kk + 1], 2)
      }
      ed <- tryCatch(TTPsplines:::.ed_S(S, Pk, lam), error = function(e) NA_real_)
      edf_rows[[length(edf_rows) + 1L]] <- data.frame(
        core = kk, lambda = lam, edf = ed, stringsAsFactors = FALSE
      )
    }
  }
  write.csv(do.call(rbind, edf_rows), file.path(out_dir, "q13_edf.csv"),
            row.names = FALSE)

  message("\n=== Bernoulli audit (", mode, ") complete ===")
  invisible(TRUE)
})
