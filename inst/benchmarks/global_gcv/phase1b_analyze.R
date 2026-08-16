# One-off analysis for PHASE1B_FULL_REPORT.md
dir <- "inst/benchmarks/global_gcv/results/phase1b_full"
fig <- "inst/benchmarks/global_gcv/figures/phase1b_full"
rep <- utils::read.csv(file.path(dir, "phase1b_replication_results.csv"),
                       stringsAsFactors = FALSE)
paired <- utils::read.csv(file.path(dir, "phase1b_paired_comparisons.csv"),
                          stringsAsFactors = FALSE)
summ <- utils::read.csv(file.path(dir, "phase1b_paired_summary.csv"),
                        stringsAsFactors = FALSE)
stab <- utils::read.csv(file.path(dir, "phase1b_initialization_stability.csv"),
                        stringsAsFactors = FALSE)
an <- utils::read.csv(file.path(dir, "phase1b_anisotropy.csv"),
                      stringsAsFactors = FALSE)
meta <- utils::read.csv(file.path(dir, "phase1b_run_meta.csv"),
                        stringsAsFactors = FALSE)

cat("META\n")
print(meta)
cat("\nN methods x scenarios:\n")
print(table(rep$method, rep$scenario))

cat("\nPAIRED SUMMARY (cGCV_default vs oracle):\n")
print(summ[grepl("cGCV_default vs tt_grid", summ$pair),
           c("scenario", "metric", "n", "mean", "median", "sd",
             "mean_lo", "mean_hi", "prop_neg", "prop_pos")])

cat("\nMULTISTART vs default (delta test_mse = best - default):\n")
print(summ[grepl("best_global vs cGCV_default", summ$pair) &
             summ$metric == "delta_test_mse",
           c("scenario", "mean", "mean_lo", "mean_hi", "prop_neg")])

for (sc in c("smooth_smooth", "strong_aniso")) {
  o <- rep[rep$scenario == sc & rep$method == "tt_grid_oracle", ]
  c0 <- rep[rep$scenario == sc & rep$method == "cGCV_default", ]
  m <- merge(
    o[, c("rep", "test_mse", "ise", "global_gcv", "log_ratio",
          "log10_l1", "log10_l2", "flat_frac", "oracle_on_edge")],
    c0[, c("rep", "test_mse", "ise", "global_gcv", "log_ratio",
           "log10_l1", "log10_l2")],
    by = "rep", suffixes = c("_o", "_c")
  )
  dm <- m$test_mse_c - m$test_mse_o
  cat(sprintf(
    "\n[%s] mean MSE oracle=%.5g cGCV=%.5g rel_dMSE=%.2f%%\n",
    sc, mean(m$test_mse_o), mean(m$test_mse_c),
    100 * mean(dm) / mean(m$test_mse_o)
  ))
  cat(sprintf(
    "  mean |dGCV|=%.4g mean dist_log10=%.3g flat_frac=%.2f edge_rate=%.2f\n",
    mean(abs(m$global_gcv_c - m$global_gcv_o)),
    mean(sqrt((m$log10_l1_c - m$log10_l1_o)^2 +
                (m$log10_l2_c - m$log10_l2_o)^2)),
    mean(m$flat_frac), mean(as.logical(m$oracle_on_edge))
  ))
  cat(sprintf("  corr(dGCV,dMSE)=%.3f\n",
              stats::cor(m$global_gcv_c - m$global_gcv_o, dm)))
}

aa <- an[an$scenario == "strong_aniso", ]
cat("\nANISOTROPY sign_pos rates:\n")
print(stats::aggregate(sign_pos ~ method, aa, function(z) mean(z == TRUE, na.rm = TRUE)))
cat("mean log_ratio by method:\n")
print(stats::aggregate(log_ratio ~ method, aa, mean, na.rm = TRUE))

st <- stab[stab$ok %in% c(TRUE, "TRUE"), ]
cat("\nINIT STABILITY:\n")
for (sc in unique(st$scenario)) {
  sds <- tapply(st$test_mse[st$scenario == sc], st$rep[st$scenario == sc],
                stats::sd, na.rm = TRUE)
  rng <- tapply(st$test_mse[st$scenario == sc], st$rep[st$scenario == sc],
                function(z) diff(range(z, na.rm = TRUE)))
  cat(sc, " median_sd=", stats::median(sds, na.rm = TRUE),
      " mean_range=", mean(rng, na.rm = TRUE), "\n")
}
cat("\nStart mean test_mse:\n")
print(stats::aggregate(test_mse ~ start + scenario, st, mean))

cat("\nFigures:\n")
print(list.files(fig))
cat("\nFailures exist:",
    file.exists(file.path(dir, "phase1b_failures.csv")), "\n")
