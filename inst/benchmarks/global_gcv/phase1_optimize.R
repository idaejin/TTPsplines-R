# Phase 1 — joint Nelder-Mead optimization of TT-gGCV (benchmark only).
# Usage: TT_GGCV_QUICK=true Rscript inst/benchmarks/global_gcv/phase1_optimize.R

source(file.path("inst", "benchmarks", "global_gcv", "phase1_config.R"))
.phase1_load_pkg()
cfg <- .phase1_config()
.phase1_ensure_dirs(cfg)

ctrl <- TTPsplines:::.tt_lab_phase1_control(seed = cfg$seed0, max_sweeps = cfg$max_sweeps)
cache <- TTPsplines:::.tt_lab_new_cache()
box <- c(cfg$coarse_lo, cfg$coarse_hi)

# Pull grid minima / cGCV if available for starts
surf_min <- file.path(cfg$out_dir, "phase1_surface_minima.csv")
basin_csv <- file.path(cfg$out_dir, "phase1_basins.csv")

all_opt <- list()
oi <- 0L

for (sc in cfg$scenarios) {
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    scenario = sc, n = cfg$n, n_test = cfg$n_test, k = cfg$k, seed = cfg$seed0
  )
  probes <- .phase1_probe_bank(des$n, cfg$M_bank, cfg$seed0 + 99L)
  rk <- .phase1_ranks(des$k)$sufficient
  common_init <- TTPsplines:::.tt_with_preserved_seed({
    tt_initialize(d = 2L, rank = rk, k = des$k, seed = des$seed)
  })

  starts <- list(
    center = c(0, 0),
    low = c(-3, -3),
    high = c(3, 3),
    aniso = c(-2, 2)
  )
  if (file.exists(surf_min)) {
    sm <- utils::read.csv(surf_min, stringsAsFactors = FALSE)
    sm <- sm[sm$scenario == sc & sm$rank_label == "sufficient", ]
    if (nrow(sm)) {
      starts$grid_best <- c(sm$log10_l1[1], sm$log10_l2[1])
    }
  }
  if (file.exists(basin_csv)) {
    bb <- utils::read.csv(basin_csv, stringsAsFactors = FALSE)
    bb <- bb[bb$scenario == sc, ]
    if (nrow(bb)) {
      starts$cgcv <- c(bb$log10_l1[1], bb$log10_l2[1])
    }
  }

  message(sprintf("=== optimize %s rank=%d starts=%d ===", sc, rk, length(starts)))
  opt <- TTPsplines:::.tt_lab_optimize_ggcv(
    design = des, rank = rk, probes = probes,
    M_search = cfg$M_search, M_final = cfg$M_final,
    starts = starts, box = box, common_init = common_init,
    mc_bank_id = "bank0", cache = cache, control = ctrl,
    epsilon_rel = cfg$epsilon_rel,
    optim_control = list(maxit = if (cfg$quick) 25L else 60L, reltol = 1e-3)
  )
  opt$scenario <- sc
  opt$rank <- rk
  oi <- oi + 1L
  all_opt[[oi]] <- opt
}

tab <- do.call(rbind, all_opt)
.phase1_checkpoint_write(tab, file.path(cfg$out_dir, "phase1_optimize.csv"))

# Best per scenario
best <- do.call(rbind, lapply(split(tab, tab$scenario), function(z) {
  z <- z[isTRUE(z$ok) | z$ok == TRUE, , drop = FALSE]
  if (!nrow(z)) return(NULL)
  z[which.min(z$global_gcv_final), , drop = FALSE]
}))
if (!is.null(best) && nrow(best)) {
  .phase1_checkpoint_write(best, file.path(cfg$out_dir, "phase1_optimize_best.csv"))
}
print(tab)
message("Optimize done.")
