# Phase 1B FULL selective orchestrator.
#
# Smoke:
#   TT_GGCV_MODE=SMOKE TT_GGCV_R=2 Rscript inst/benchmarks/global_gcv/phase1b_run_full.R
# Full:
#   TT_GGCV_MODE=FULL TT_GGCV_R=50 \
#   TT_GGCV_SCENARIOS=smooth_smooth,strong_aniso \
#   Rscript inst/benchmarks/global_gcv/phase1b_run_full.R

t_wall0 <- proc.time()[["elapsed"]]
source(file.path("inst", "benchmarks", "global_gcv", "phase1b_config.R"))
source(file.path("inst", "benchmarks", "global_gcv", "phase1b_helpers.R"))
.phase1_load_pkg()
`%||%` <- function(a, b) if (!is.null(a)) a else b

cfg <- .phase1b_config()
.phase1b_ensure_dirs(cfg)

message("=== Phase 1B ", cfg$mode, " ===")
message("scenarios: ", paste(cfg$scenarios, collapse = ", "))
message("R=", cfg$n_rep, " n=", cfg$n, " n_test=", cfg$n_test, " k=", cfg$k)
message("M_search=", cfg$M_search, " M_final=", cfg$M_final)
message("out_dir: ", cfg$out_dir)

path_rep <- file.path(cfg$out_dir, "phase1b_replication_results.csv")
path_stab <- file.path(cfg$out_dir, "phase1b_initialization_stability.csv")
path_aniso <- file.path(cfg$out_dir, "phase1b_anisotropy.csv")
path_fail <- file.path(cfg$out_dir, "phase1b_failures.csv")
path_time <- file.path(cfg$out_dir, "phase1b_timing.csv")
path_surf <- file.path(cfg$out_dir, "phase1b_surfaces_sample.csv")

done <- unique(.phase1b_done_keys(path_rep, c("scenario", "rep")))
message("Already completed rep-keys: ", length(done))

timing_rows <- list()
fail_rows <- list()
n_done_now <- 0L
n_skip <- 0L
n_fail <- 0L
n_total <- length(cfg$scenarios) * cfg$n_rep
t_jobs <- numeric(0)

for (sc in cfg$scenarios) {
  for (rep in seq_len(cfg$n_rep)) {
    key <- paste(sc, rep, sep = "|")
    if (key %in% done) {
      n_skip <- n_skip + 1L
      next
    }
    message(sprintf(
      "[%s] %s rep %d/%d ...",
      format(Sys.time(), "%H:%M:%S"), sc, rep, cfg$n_rep
    ))
    t0 <- proc.time()[["elapsed"]]
    out <- tryCatch(
      .phase1b_one_rep(sc, rep, cfg),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    if (inherits(out, "error")) {
      n_fail <- n_fail + 1L
      fr <- data.frame(
        mode = cfg$mode, scenario = sc, rep = rep,
        error = conditionMessage(out), elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      )
      fail_rows[[length(fail_rows) + 1L]] <- fr
      .phase1b_append_csv(fr, path_fail)
      message("  FAIL: ", conditionMessage(out))
      next
    }
    .phase1b_append_csv(out$replication, path_rep)
    .phase1b_append_csv(out$stability, path_stab)
    .phase1b_append_csv(out$anisotropy, path_aniso)
    # Keep a light surface sample (rep 1 only) for figures
    if (rep == 1L && !is.null(out$surface)) {
      .phase1b_append_csv(out$surface, path_surf)
    }
    n_done_now <- n_done_now + 1L
    t_jobs <- c(t_jobs, elapsed)
    eta <- if (length(t_jobs)) {
      remaining <- n_total - (length(done) + n_skip + n_done_now + n_fail)
      # approximate: keys left among incomplete
      remaining <- max(0, n_total - length(unique(c(done, key))) - n_skip)
      mean(t_jobs) * max(0, n_total - length(done) - n_done_now - n_skip - n_fail)
    } else NA_real_
    message(sprintf(
      "  ok in %.1fs | flat=%.0f%% edge=%s | ETA~%.0fs",
      elapsed, 100 * (out$replication$flat_frac[1] %||% NA),
      out$replication$oracle_on_edge[1], eta
    ))
    timing_rows[[length(timing_rows) + 1L]] <- data.frame(
      mode = cfg$mode, scenario = sc, rep = rep,
      elapsed_sec = elapsed, status = "ok", stringsAsFactors = FALSE
    )
  }
}

if (length(timing_rows)) {
  .phase1b_append_csv(do.call(rbind, timing_rows), path_time)
}

# -------------------- Aggregate paired comparisons --------------------
if (!file.exists(path_rep)) {
  stop("No replication results written — nothing to aggregate.", call. = FALSE)
}
rep_tab <- utils::read.csv(path_rep, stringsAsFactors = FALSE)

# Prefer primary methods for paired analysis
pick <- function(method) {
  rep_tab[rep_tab$method == method, , drop = FALSE]
}
oracle <- pick("tt_grid_oracle")
cg_def <- pick("cGCV_default")
cg_best <- pick("cGCV_best_global")
full <- pick("full_tp_gcv")
nelder <- pick("tt_nelder")

merge_pair <- function(a, b, name_a, name_b) {
  m <- merge(
    a[, c("scenario", "rep", "test_mse", "ise", "global_gcv",
          "log10_l1", "log10_l2", "log_ratio", "gdf", "gdf_se", "dist_to_oracle")],
    b[, c("scenario", "rep", "test_mse", "ise", "global_gcv",
          "log10_l1", "log10_l2", "log_ratio", "gdf", "gdf_se")],
    by = c("scenario", "rep"),
    suffixes = c("_a", "_b")
  )
  if (!nrow(m)) return(NULL)
  m$pair <- paste(name_a, "vs", name_b)
  m$delta_test_mse <- m$test_mse_a - m$test_mse_b
  m$delta_ise <- m$ise_a - m$ise_b
  m$delta_gcv <- m$global_gcv_a - m$global_gcv_b
  m$dist_log10 <- sqrt((m$log10_l1_a - m$log10_l1_b)^2 +
                         (m$log10_l2_a - m$log10_l2_b)^2)
  # MC-aware: flag tiny GCV diffs relative to oracle SE if available
  m$gcv_diff_within_mc <- abs(m$delta_gcv) <
    pmax(2 * (m$gdf_se_b %||% 0), 1e-4) # crude placeholder; refined below
  m
}

pairs <- list(
  merge_pair(cg_def, oracle, "cGCV_default", "tt_grid_oracle"),
  merge_pair(cg_best, oracle, "cGCV_best_global", "tt_grid_oracle"),
  merge_pair(cg_best, cg_def, "cGCV_best_global", "cGCV_default")
)
if (nrow(nelder)) {
  pairs[[length(pairs) + 1L]] <- merge_pair(nelder, oracle, "tt_nelder", "tt_grid_oracle")
}
if (nrow(full)) {
  pairs[[length(pairs) + 1L]] <- merge_pair(full, oracle, "full_tp_gcv", "tt_grid_oracle")
}
paired <- do.call(rbind, Filter(Negate(is.null), pairs))
utils::write.csv(
  paired,
  file.path(cfg$out_dir, "phase1b_paired_comparisons.csv"),
  row.names = FALSE
)

# Summary stats
summ_list <- list()
for (pr in unique(paired$pair)) {
  for (sc in unique(paired$scenario)) {
    sub <- paired[paired$pair == pr & paired$scenario == sc, ]
    if (!nrow(sub)) next
    for (met in c("delta_test_mse", "delta_ise", "delta_gcv", "dist_log10")) {
      st <- .phase1b_paired_stats(sub[[met]], label = paste(pr, sc, met, sep = "|"))
      st$pair <- pr
      st$scenario <- sc
      st$metric <- met
      summ_list[[length(summ_list) + 1L]] <- st
    }
  }
}
summ <- do.call(rbind, summ_list)
utils::write.csv(summ, file.path(cfg$out_dir, "phase1b_paired_summary.csv"), row.names = FALSE)

# -------------------- Figures --------------------
fig <- function(name, expr) {
  png(file.path(cfg$fig_dir, name), width = 800, height = 520)
  tryCatch(expr, error = function(e) {
    plot.new(); title(paste("figure failed:", conditionMessage(e)))
  })
  dev.off()
}

fig("phase1b_paired_test_mse.png", {
  sub <- paired[paired$pair == "cGCV_default vs tt_grid_oracle", ]
  boxplot(delta_test_mse ~ scenario, data = sub,
          main = "Paired test-MSE: cGCV_default - oracle",
          ylab = "delta test-MSE (<0 => cGCV better)")
  abline(h = 0, lty = 2, col = "red")
})
fig("phase1b_paired_ise.png", {
  sub <- paired[paired$pair == "cGCV_default vs tt_grid_oracle", ]
  boxplot(delta_ise ~ scenario, data = sub,
          main = "Paired ISE: cGCV_default - oracle",
          ylab = "delta ISE (<0 => cGCV better)")
  abline(h = 0, lty = 2, col = "red")
})
fig("phase1b_delta_gcv_vs_delta_mse.png", {
  sub <- paired[paired$pair == "cGCV_default vs tt_grid_oracle", ]
  plot(sub$delta_gcv, sub$delta_test_mse, col = as.factor(sub$scenario),
       pch = 16, xlab = "delta TT-gGCV (cGCV - oracle)",
       ylab = "delta test-MSE (cGCV - oracle)",
       main = "GCV gap vs predictive gap")
  abline(h = 0, v = 0, lty = 2)
  legend("topleft", legend = levels(factor(sub$scenario)),
         col = seq_along(levels(factor(sub$scenario))), pch = 16, bty = "n")
})
fig("phase1b_lambda_cloud.png", {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (sc in cfg$scenarios) {
    plot(NA, xlim = c(-5, 5), ylim = c(-5, 5),
         xlab = "log10(l1)", ylab = "log10(l2)", main = sc)
    o <- oracle[oracle$scenario == sc, ]
    c0 <- cg_def[cg_def$scenario == sc, ]
    points(o$log10_l1, o$log10_l2, pch = 8, col = "red")
    points(c0$log10_l1, c0$log10_l2, pch = 16, col = adjustcolor("darkgreen", 0.5))
    legend("topright", c("oracle", "cGCV_default"), pch = c(8, 16),
           col = c("red", "darkgreen"), bty = "n", cex = 0.8)
  }
  par(op)
})
fig("phase1b_anisotropy.png", {
  if (file.exists(path_aniso)) {
    an <- utils::read.csv(path_aniso, stringsAsFactors = FALSE)
    an <- an[an$scenario == "strong_aniso", ]
    boxplot(log_ratio ~ method, data = an, las = 2, cex.axis = 0.7,
            main = "strong_aniso: log10(l1/l2)",
            ylab = "log10(lambda1/lambda2)")
    abline(h = 0, lty = 2)
  } else {
    plot.new(); title("no anisotropy file")
  }
})
fig("phase1b_cgcv_starts.png", {
  if (file.exists(path_stab)) {
    st <- utils::read.csv(path_stab, stringsAsFactors = FALSE)
    st <- st[isTRUE(st$ok) | st$ok == TRUE, ]
    boxplot(test_mse ~ start + scenario, data = st, las = 2, cex.axis = 0.65,
            main = "cGCV test-MSE by start", ylab = "test-MSE")
  } else plot.new()
})
fig("phase1b_cost_vs_performance.png", {
  # approximate: use elapsed if present
  if ("elapsed_sec" %in% names(rep_tab)) {
    agg <- aggregate(cbind(test_mse, elapsed_sec) ~ method + scenario,
                     data = rep_tab, mean, na.rm = TRUE)
    plot(agg$elapsed_sec, agg$test_mse, pch = 16,
         col = as.factor(agg$scenario),
         xlab = "mean elapsed (sec)", ylab = "mean test-MSE",
         main = "Cost vs performance")
    text(agg$elapsed_sec, agg$test_mse, labels = agg$method, cex = 0.6, pos = 3)
  } else plot.new()
})
fig("phase1b_flat_region.png", {
  if (file.exists(path_surf)) {
    sf <- utils::read.csv(path_surf, stringsAsFactors = FALSE)
    op <- par(mfrow = c(1, length(unique(sf$scenario))), mar = c(4, 4, 3, 1))
    for (sc in unique(sf$scenario)) {
      s <- sf[sf$scenario == sc, ]
      plot(s$log10_l1, s$log10_l2, pch = ifelse(s$flat_1pct, 15, 1),
           col = ifelse(s$flat_1pct, "orange", "gray50"),
           xlab = "log10(l1)", ylab = "log10(l2)",
           main = paste(sc, "1% flat region"))
      i <- which.min(s$global_gcv)
      points(s$log10_l1[i], s$log10_l2[i], pch = 8, col = "red", cex = 1.4)
    }
    par(op)
  } else plot.new()
})

# -------------------- Decision helper printed to console --------------------
dec_lines <- character(0)
for (sc in cfg$scenarios) {
  sub <- paired[paired$pair == "cGCV_default vs tt_grid_oracle" &
                  paired$scenario == sc, ]
  st_mse <- .phase1b_paired_stats(sub$delta_test_mse)
  st_ise <- .phase1b_paired_stats(sub$delta_ise)
  st_gcv <- .phase1b_paired_stats(sub$delta_gcv)
  # Multistart gain
  sub2 <- paired[paired$pair == "cGCV_best_global vs cGCV_default" &
                   paired$scenario == sc, ]
  st_ms <- .phase1b_paired_stats(sub2$delta_test_mse)
  dec_lines <- c(dec_lines, sprintf(
    paste0("[%s] mean dMSE(cGCV-oracle)=%.4g (CI %.4g,%.4g); ",
           "prop cGCV better=%.2f; dISE=%.4g; multistart dMSE=%.4g"),
    sc, st_mse$mean, st_mse$mean_lo, st_mse$mean_hi,
    st_mse$prop_neg, st_ise$mean, st_ms$mean
  ))
}
writeLines(dec_lines, file.path(cfg$out_dir, "phase1b_decision_snippets.txt"))

wall <- proc.time()[["elapsed"]] - t_wall0
meta <- data.frame(
  mode = cfg$mode,
  wall_sec = wall,
  n_rep_requested = cfg$n_rep,
  n_completed_keys = length(unique(.phase1b_done_keys(path_rep))),
  n_done_this_run = n_done_now,
  n_skipped = n_skip,
  n_failed = n_fail,
  scenarios = paste(cfg$scenarios, collapse = ","),
  M_search = cfg$M_search,
  M_final = cfg$M_final,
  n = cfg$n,
  k = cfg$k,
  stringsAsFactors = FALSE
)
utils::write.csv(meta, file.path(cfg$out_dir, "phase1b_run_meta.csv"), row.names = FALSE)

message("\n=== Phase 1B ", cfg$mode, " finished ===")
print(meta)
message(paste(dec_lines, collapse = "\n"))
message("Artifacts in ", cfg$out_dir)
