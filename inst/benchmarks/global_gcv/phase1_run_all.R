# Phase 1 — orchestrator.
# Usage:
#   TT_GGCV_QUICK=true  Rscript inst/benchmarks/global_gcv/phase1_run_all.R
#   TT_GGCV_QUICK=false Rscript inst/benchmarks/global_gcv/phase1_run_all.R

t_wall0 <- proc.time()[["elapsed"]]
source(file.path("inst", "benchmarks", "global_gcv", "phase1_config.R"))
.phase1_load_pkg()
cfg <- .phase1_config()
.phase1_ensure_dirs(cfg)

message("=== Phase 1 config ===")
print(cfg[c("quick", "n", "k", "M_search", "M_final", "n_rep_mc", "coarse_by")])

run_step <- function(label, script) {
  message("\n######## ", label, " ########")
  t0 <- proc.time()[["elapsed"]]
  status <- tryCatch({
    source(script, local = new.env(parent = globalenv()))
    "ok"
  }, error = function(e) {
    message("ERROR in ", label, ": ", conditionMessage(e))
    paste0("error: ", conditionMessage(e))
  })
  elapsed <- proc.time()[["elapsed"]] - t0
  data.frame(step = label, status = status, elapsed_sec = elapsed,
             stringsAsFactors = FALSE)
}

timing <- rbind(
  run_step("designs", file.path("inst", "benchmarks", "global_gcv", "phase1_designs.R")),
  run_step("surfaces", file.path("inst", "benchmarks", "global_gcv", "phase1_surfaces.R")),
  run_step("basins", file.path("inst", "benchmarks", "global_gcv", "phase1_basins.R")),
  run_step("optimize", file.path("inst", "benchmarks", "global_gcv", "phase1_optimize.R")),
  run_step("mc", file.path("inst", "benchmarks", "global_gcv", "phase1_mc.R"))
)

timing$total_wall_sec <- proc.time()[["elapsed"]] - t_wall0
utils::write.csv(
  timing,
  file.path(cfg$out_dir, "phase1_timing.csv"),
  row.names = FALSE
)

# Location comparison figure if artifacts exist
mins <- file.path(cfg$out_dir, "phase1_surface_minima.csv")
bas <- file.path(cfg$out_dir, "phase1_basins.csv")
opt <- file.path(cfg$out_dir, "phase1_optimize_best.csv")
if (file.exists(mins) && file.exists(bas)) {
  m <- utils::read.csv(mins, stringsAsFactors = FALSE)
  b <- utils::read.csv(bas, stringsAsFactors = FALSE)
  o <- if (file.exists(opt)) utils::read.csv(opt, stringsAsFactors = FALSE) else NULL
  png(file.path(cfg$fig_dir, "phase1_method_locations.png"), width = 900, height = 700)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (sc in unique(m$scenario)) {
    plot(NA, xlim = c(cfg$coarse_lo, cfg$coarse_hi),
         ylim = c(cfg$coarse_lo, cfg$coarse_hi),
         xlab = "log10(l1)", ylab = "log10(l2)", main = sc)
    ms <- m[m$scenario == sc & m$rank_label == "sufficient", ]
    points(ms$full_min_log10_l1, ms$full_min_log10_l2, pch = 1, col = "blue", cex = 1.4)
    points(ms$log10_l1, ms$log10_l2, pch = 8, col = "red", cex = 1.3)
    bs <- b[b$scenario == sc, ]
    points(bs$log10_l1, bs$log10_l2, pch = 4, col = "darkgreen")
    if (!is.null(o)) {
      os <- o[o$scenario == sc, ]
      if (nrow(os)) points(os$log10_l1, os$log10_l2, pch = 17, col = "purple", cex = 1.2)
    }
    legend("topright",
           c("full-TP", "TT-grid", "cGCV", "Nelder"),
           pch = c(1, 8, 4, 17),
           col = c("blue", "red", "darkgreen", "purple"),
           bty = "n", cex = 0.7)
  }
  par(op)
  dev.off()
}

message("\n=== Phase 1 finished ===")
print(timing)
message(sprintf("Total wall time: %.1fs (quick=%s)", timing$total_wall_sec[1], cfg$quick))
