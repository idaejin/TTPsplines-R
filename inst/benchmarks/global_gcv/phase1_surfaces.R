# Phase 1 — TT-gGCV / full-TP surfaces (cold + optional continuation).
# Usage: TT_GGCV_QUICK=true Rscript inst/benchmarks/global_gcv/phase1_surfaces.R

source(file.path("inst", "benchmarks", "global_gcv", "phase1_config.R"))
.phase1_load_pkg()
cfg <- .phase1_config()
.phase1_ensure_dirs(cfg)

coarse <- seq(cfg$coarse_lo, cfg$coarse_hi, by = cfg$coarse_by)
ctrl <- TTPsplines:::.tt_lab_phase1_control(seed = cfg$seed0, max_sweeps = cfg$max_sweeps)
cache <- TTPsplines:::.tt_lab_new_cache()

# Optional: restrict to scenarios with interior candidate if designs csv exists
des_csv <- file.path(cfg$out_dir, "phase1_designs.csv")
focus <- cfg$scenarios
if (file.exists(des_csv)) {
  dtab <- utils::read.csv(des_csv, stringsAsFactors = FALSE)
  message("Loaded design scout; interior candidates: ",
          paste(dtab$scenario[dtab$interior_min_candidate], collapse = ", "))
}

all_rows <- list()
minima_rows <- list()
ai <- 0L
mi <- 0L

for (sc in focus) {
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    scenario = sc, n = cfg$n, n_test = cfg$n_test, k = cfg$k, seed = cfg$seed0
  )
  probes <- .phase1_probe_bank(des$n, cfg$M_bank, cfg$seed0 + 99L)
  ranks <- .phase1_ranks(des$k)
  common_init_suf <- TTPsplines:::.tt_with_preserved_seed({
    tt_initialize(d = 2L, rank = ranks$sufficient, k = des$k, seed = des$seed)
  })

  # Coarse cold surfaces: sufficient + restricted
  for (rk_name in c("sufficient", "restricted")) {
    rk <- ranks[[rk_name]]
    common_init <- TTPsplines:::.tt_with_preserved_seed({
      tt_initialize(d = 2L, rank = rk, k = des$k, seed = des$seed)
    })
    message(sprintf("=== %s rank=%s (%d) cold coarse ===", sc, rk_name, rk))
    tab <- TTPsplines:::.tt_lab_eval_lambda_grid(
      design = des, rank = rk, log10_grid = coarse, probes = probes,
      M = cfg$M_search, init_policy = "cold_common", traverse = "row_asc",
      common_init = common_init, mc_bank_id = "bank0", cache = cache,
      control = ctrl, epsilon_rel = cfg$epsilon_rel
    )
    tab$rank_label <- rk_name
    tab$grid <- "coarse"
    ai <- ai + 1L
    all_rows[[ai]] <- tab

    # Continuation serpentine on sufficient only (cost control)
    if (identical(rk_name, "sufficient") && !cfg$quick) {
      for (tr in c("row_asc", "row_desc", "serpentine")) {
        message(sprintf("=== %s continuation %s ===", sc, tr))
        tabc <- TTPsplines:::.tt_lab_eval_lambda_grid(
          design = des, rank = rk, log10_grid = coarse, probes = probes,
          M = cfg$M_search, init_policy = "continuation", traverse = tr,
          common_init = NULL, mc_bank_id = "bank0", cache = cache,
          control = ctrl, epsilon_rel = cfg$epsilon_rel
        )
        tabc$rank_label <- rk_name
        tabc$grid <- paste0("coarse_", tr)
        ai <- ai + 1L
        all_rows[[ai]] <- tabc
      }
    } else if (identical(rk_name, "sufficient") && cfg$quick) {
      # QUICK: one continuation serpentine only
      tabc <- TTPsplines:::.tt_lab_eval_lambda_grid(
        design = des, rank = rk, log10_grid = coarse, probes = probes,
        M = cfg$M_search, init_policy = "continuation", traverse = "serpentine",
        common_init = NULL, mc_bank_id = "bank0", cache = cache,
        control = ctrl, epsilon_rel = cfg$epsilon_rel
      )
      tabc$rank_label <- rk_name
      tabc$grid <- "coarse_serpentine"
      ai <- ai + 1L
      all_rows[[ai]] <- tabc
    }

    # Minima (cold)
    ok <- is.finite(tab$global_gcv_tt)
    if (any(ok)) {
      i <- which(ok)[which.min(tab$global_gcv_tt[ok])]
      # Final eval M_final
      evf <- TTPsplines:::.tt_lab_eval_ggcv(
        lambda = c(tab$lambda1[i], tab$lambda2[i]), design = des, rank = rk,
        probes = probes, M = cfg$M_final, init_cores = common_init,
        init_policy = "cold_common", mc_bank_id = "bank0", cache = cache,
        control = ctrl, epsilon_rel = cfg$epsilon_rel
      )
      ifull <- which.min(tab$gcv_full)
      mi <- mi + 1L
      minima_rows[[mi]] <- data.frame(
        scenario = sc, rank_label = rk_name, rank = rk,
        method = "grid_TT_gGCV",
        log10_l1 = tab$log10_l1[i], log10_l2 = tab$log10_l2[i],
        global_gcv = evf$global_gcv, gdf = evf$gdf, gdf_se = evf$gdf_mc_se,
        rss = evf$rss, rmse_test = evf$rmse_test, ise_true = evf$ise_true,
        flat_1pct_n = sum(tab$flat_1pct),
        full_min_log10_l1 = tab$log10_l1[ifull],
        full_min_log10_l2 = tab$log10_l2[ifull],
        full_min_gcv = tab$gcv_full[ifull],
        stringsAsFactors = FALSE
      )
    }
  }

  # Heatmap for sufficient cold coarse
  tab_suf <- all_rows[[length(all_rows) - (if (cfg$quick) 1L else 0L)]]
  # find last sufficient cold coarse for this scenario
  for (j in rev(seq_along(all_rows))) {
    if (identical(all_rows[[j]]$scenario[1], sc) &&
        identical(all_rows[[j]]$rank_label[1], "sufficient") &&
        identical(all_rows[[j]]$init_policy[1], "cold_common")) {
      tab_suf <- all_rows[[j]]
      break
    }
  }
  png(file.path(cfg$fig_dir, sprintf("phase1_heatmap_%s.png", sc)),
      width = 900, height = 420)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (colnm in c("global_gcv_tt", "gcv_full")) {
    m <- matrix(NA_real_, length(coarse), length(coarse),
                dimnames = list(coarse, coarse))
    for (r in seq_len(nrow(tab_suf))) {
      m[as.character(tab_suf$log10_l1[r]), as.character(tab_suf$log10_l2[r])] <-
        tab_suf[[colnm]][r]
    }
    image(coarse, coarse, m, xlab = "log10(l1)", ylab = "log10(l2)",
          main = paste(sc, colnm), col = hcl.colors(20, "Blues"))
    # flat region
    flat <- tab_suf[tab_suf$flat_1pct %in% TRUE, ]
    if (nrow(flat)) {
      points(flat$log10_l1, flat$log10_l2, pch = 15, col = adjustcolor("orange", 0.5))
    }
    i <- which.min(tab_suf[[colnm]])
    points(tab_suf$log10_l1[i], tab_suf$log10_l2[i], pch = 8, col = "red", cex = 1.3)
  }
  par(op)
  dev.off()
}

surf <- do.call(rbind, all_rows)
mins <- do.call(rbind, minima_rows)
.phase1_checkpoint_write(surf, file.path(cfg$out_dir, "phase1_surfaces.csv"))
.phase1_checkpoint_write(mins, file.path(cfg$out_dir, "phase1_surface_minima.csv"))

# Cold vs continuation delta (sufficient)
if (any(surf$init_policy == "continuation")) {
  cold <- surf[surf$init_policy == "cold_common" & surf$rank_label == "sufficient", ]
  cont <- surf[surf$init_policy == "continuation" & surf$rank_label == "sufficient", ]
  # merge on scenario + log lambda
  mer <- merge(
    cold[, c("scenario", "log10_l1", "log10_l2", "global_gcv_tt")],
    cont[, c("scenario", "log10_l1", "log10_l2", "global_gcv_tt", "traverse")],
    by = c("scenario", "log10_l1", "log10_l2"),
    suffixes = c("_cold", "_cont")
  )
  mer$delta_gcv <- mer$global_gcv_tt_cont - mer$global_gcv_tt_cold
  .phase1_checkpoint_write(mer, file.path(cfg$out_dir, "phase1_cold_vs_continuation.csv"))
  png(file.path(cfg$fig_dir, "phase1_cold_vs_continuation.png"), width = 640, height = 420)
  boxplot(delta_gcv ~ scenario, data = mer, main = "TT-gGCV continuation - cold",
          ylab = "delta global_gcv")
  abline(h = 0, lty = 2)
  dev.off()
}

message("Surfaces done. Rows=", nrow(surf))
