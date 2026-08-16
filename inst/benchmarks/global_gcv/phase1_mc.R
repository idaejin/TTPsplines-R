# Phase 1 — compact Monte Carlo across replicates (prediction / GCV / stability).
# Usage: TT_GGCV_QUICK=true Rscript inst/benchmarks/global_gcv/phase1_mc.R

source(file.path("inst", "benchmarks", "global_gcv", "phase1_config.R"))
.phase1_load_pkg()
cfg <- .phase1_config()
.phase1_ensure_dirs(cfg)

# Methods compared per replicate:
#   full_grid (full-TP GCV min on coarse), tt_grid, cgcv_default, tt_optim (if available)

rows <- list()
ri <- 0L
t_all0 <- proc.time()[["elapsed"]]

for (rep in seq_len(cfg$n_rep_mc)) {
  seed <- cfg$seed0 + (rep - 1L) * 97L
  for (sc in cfg$scenarios) {
    des <- TTPsplines:::.tt_lab_phase1_make_design(
      scenario = sc, n = cfg$n, n_test = cfg$n_test, k = cfg$k, seed = seed
    )
    probes <- .phase1_probe_bank(des$n, cfg$M_bank, seed + 99L)
    ranks <- .phase1_ranks(des$k)
    rk <- ranks$sufficient
    ctrl <- TTPsplines:::.tt_lab_phase1_control(seed = seed, max_sweeps = cfg$max_sweeps)
    common_init <- TTPsplines:::.tt_with_preserved_seed({
      tt_initialize(d = 2L, rank = rk, k = des$k, seed = seed)
    })
    cache <- TTPsplines:::.tt_lab_new_cache()
    coarse <- seq(cfg$coarse_lo, cfg$coarse_hi, by = cfg$coarse_by)

    # --- full-TP grid min ---
    best_full <- list(gcv = Inf, th = c(0, 0))
    for (a in coarse) for (b in coarse) {
      full <- TTPsplines:::.tt_lab_eval_full_gcv(des, 10^c(a, b))
      if (is.finite(full$gcv) && full$gcv < best_full$gcv) {
        best_full <- list(gcv = full$gcv, th = c(a, b), rss = full$rss, edf = full$edf)
      }
    }

    # --- TT grid min (cold) ---
    tab <- TTPsplines:::.tt_lab_eval_lambda_grid(
      design = des, rank = rk, log10_grid = coarse, probes = probes,
      M = cfg$M_search, init_policy = "cold_common", traverse = "row_asc",
      common_init = common_init, mc_bank_id = "bank0", cache = cache,
      control = ctrl, epsilon_rel = cfg$epsilon_rel
    )
    ok <- is.finite(tab$global_gcv_tt)
    i_tt <- which(ok)[which.min(tab$global_gcv_tt[ok])]
    ev_tt <- TTPsplines:::.tt_lab_eval_ggcv(
      lambda = c(tab$lambda1[i_tt], tab$lambda2[i_tt]), design = des, rank = rk,
      probes = probes, M = cfg$M_final, init_cores = common_init,
      init_policy = "cold_common", mc_bank_id = "bank0", cache = cache,
      control = ctrl, epsilon_rel = cfg$epsilon_rel
    )

    # --- cGCV default ---
    ctrl_c <- tt_control(
      max_sweeps = cfg$max_sweeps, outer_maxit = if (cfg$quick) 6L else 10L,
      compute_edf = FALSE, seed = seed, backend = "R",
      cgcv_update = "outer_simultaneous", lambda_start = 1,
      warn_lambda_boundary = FALSE
    )
    fit_c <- TTPsplines:::.tt_with_preserved_seed({
      ttps(
        des$y, des$X, family = gaussian(), rank = rk, k = des$k,
        lambda = "cGCV", optimizer = "ALS", backend = "R",
        init = TTPsplines:::.tt_clone_cores(common_init), control = ctrl_c
      )
    })
    ev_c <- TTPsplines:::.tt_lab_eval_ggcv(
      lambda = as.numeric(fit_c$lambda), design = des, rank = rk,
      probes = probes, M = cfg$M_final, init_cores = common_init,
      init_policy = "cold_common", mc_bank_id = "bank0", cache = cache,
      control = ctrl, epsilon_rel = cfg$epsilon_rel
    )
    mu_c <- as.numeric(predict(fit_c, newdata = des$X_test, type = "response"))

    # Score full-TP solution under TT-gGCV for comparable delta
    ev_full_as_tt <- TTPsplines:::.tt_lab_eval_ggcv(
      lambda = 10^best_full$th, design = des, rank = rk,
      probes = probes, M = cfg$M_final, init_cores = common_init,
      init_policy = "cold_common", mc_bank_id = "bank0", cache = cache,
      control = ctrl, epsilon_rel = cfg$epsilon_rel
    )

    oracle_gcv <- ev_tt$global_gcv

    add_row <- function(method, th, ev, rmse_test, ise_true, extra = list()) {
      ri <<- ri + 1L
      rows[[ri]] <<- data.frame(
        rep = rep, seed = seed, scenario = sc, method = method,
        log10_l1 = th[1], log10_l2 = th[2],
        lambda1 = 10^th[1], lambda2 = 10^th[2],
        global_gcv = ev$global_gcv,
        delta_gcv_vs_oracle = ev$global_gcv - oracle_gcv,
        dist_log10 = sqrt(sum((th - c(tab$log10_l1[i_tt], tab$log10_l2[i_tt]))^2)),
        gdf = ev$gdf, gdf_se = ev$gdf_mc_se, rss = ev$rss,
        rmse_test = rmse_test, ise_true = ise_true,
        flat_1pct_n = sum(tab$flat_1pct),
        stringsAsFactors = FALSE
      )
    }

    add_row(
      "tt_grid_oracle", c(tab$log10_l1[i_tt], tab$log10_l2[i_tt]), ev_tt,
      ev_tt$rmse_test, ev_tt$ise_true
    )
    add_row(
      "full_tp_gcv", best_full$th, ev_full_as_tt,
      ev_full_as_tt$rmse_test, ev_full_as_tt$ise_true
    )
    add_row(
      "cgcv", log10(as.numeric(fit_c$lambda)), ev_c,
      sqrt(mean((mu_c - des$y_test)^2)), mean((mu_c - des$f_test)^2)
    )

    # Optional: one Nelder-Mead from center (costly)
    if (!cfg$quick || sc %in% c("smooth_smooth", "strong_aniso")) {
      opt <- TTPsplines:::.tt_lab_optimize_ggcv(
        design = des, rank = rk, probes = probes,
        M_search = cfg$M_search, M_final = cfg$M_final,
        starts = list(c(0, 0)), box = c(cfg$coarse_lo, cfg$coarse_hi),
        common_init = common_init, mc_bank_id = "bank0", cache = cache,
        control = ctrl, epsilon_rel = cfg$epsilon_rel,
        optim_control = list(maxit = if (cfg$quick) 20L else 40L, reltol = 1e-3)
      )
      if (isTRUE(opt$ok[1]) || isTRUE(as.logical(opt$ok[1]))) {
        ev_o <- TTPsplines:::.tt_lab_eval_ggcv(
          lambda = c(opt$lambda1[1], opt$lambda2[1]), design = des, rank = rk,
          probes = probes, M = cfg$M_final, init_cores = common_init,
          init_policy = "cold_common", mc_bank_id = "bank0", cache = cache,
          control = ctrl, epsilon_rel = cfg$epsilon_rel
        )
        add_row(
          "tt_nelder", c(opt$log10_l1[1], opt$log10_l2[1]), ev_o,
          ev_o$rmse_test, ev_o$ise_true
        )
      }
    }

    message(sprintf(
      "rep=%d %s oracle_gcv=%.4g cgcv_delta=%.4g",
      rep, sc, oracle_gcv, ev_c$global_gcv - oracle_gcv
    ))
    # incremental checkpoint
    .phase1_checkpoint_write(
      do.call(rbind, rows),
      file.path(cfg$out_dir, "phase1_mc.csv")
    )
  }
}

tab <- do.call(rbind, rows)
.phase1_checkpoint_write(tab, file.path(cfg$out_dir, "phase1_mc.csv"))

# Summary by method
summ <- aggregate(
  cbind(delta_gcv_vs_oracle, dist_log10, rmse_test, ise_true, global_gcv) ~ method + scenario,
  data = tab, FUN = mean, na.rm = TRUE
)
.phase1_checkpoint_write(summ, file.path(cfg$out_dir, "phase1_mc_summary.csv"))

png(file.path(cfg$fig_dir, "phase1_prediction_by_method.png"), width = 800, height = 480)
op <- par(mar = c(8, 4, 3, 1))
boxplot(ise_true ~ method + scenario, data = tab, las = 2, cex.axis = 0.7,
        main = "ISE vs truth (test)", ylab = "ISE")
par(op)
dev.off()

png(file.path(cfg$fig_dir, "phase1_delta_gcv_by_method.png"), width = 800, height = 480)
op <- par(mar = c(8, 4, 3, 1))
boxplot(delta_gcv_vs_oracle ~ method + scenario, data = tab, las = 2, cex.axis = 0.7,
        main = "Delta TT-gGCV vs grid oracle", ylab = "delta")
abline(h = 0, lty = 2)
par(op)
dev.off()

elapsed <- proc.time()[["elapsed"]] - t_all0
message(sprintf("MC compact done in %.1fs; rows=%d", elapsed, nrow(tab)))
print(summ)
