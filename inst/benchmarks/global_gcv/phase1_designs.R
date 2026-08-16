# Phase 1 — designs + coarse full-TP scout for interior minima.
# Usage: TT_GGCV_QUICK=true Rscript inst/benchmarks/global_gcv/phase1_designs.R

source(file.path("inst", "benchmarks", "global_gcv", "phase1_config.R"))
.phase1_load_pkg()
cfg <- .phase1_config()
.phase1_ensure_dirs(cfg)

log_grid <- seq(cfg$coarse_lo, cfg$coarse_hi, by = cfg$coarse_by)
rows <- list()
idx <- 0L

for (sc in cfg$scenarios) {
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    scenario = sc, n = cfg$n, n_test = cfg$n_test, k = cfg$k, seed = cfg$seed0
  )
  # Scout full-TP GCV on coarse grid
  best_gcv <- Inf
  best_th <- c(NA_real_, NA_real_)
  boundary_hit <- TRUE
  gcv_vals <- numeric(0)
  for (a in log_grid) {
    for (b in log_grid) {
      full <- TTPsplines:::.tt_lab_eval_full_gcv(des, 10^c(a, b))
      gcv_vals <- c(gcv_vals, full$gcv)
      if (is.finite(full$gcv) && full$gcv < best_gcv) {
        best_gcv <- full$gcv
        best_th <- c(a, b)
      }
    }
  }
  # Expand box check: also evaluate one step outside if min on edge
  on_edge <- any(best_th %in% c(cfg$coarse_lo, cfg$coarse_hi))
  if (on_edge) {
    widen <- c(cfg$coarse_lo - cfg$coarse_by, cfg$coarse_hi + cfg$coarse_by)
    for (a in widen) {
      for (b in widen) {
        if (a < cfg$coarse_lo || a > cfg$coarse_hi ||
            b < cfg$coarse_lo || b > cfg$coarse_hi) {
          full <- tryCatch(
            TTPsplines:::.tt_lab_eval_full_gcv(des, 10^c(a, b)),
            error = function(e) NULL
          )
          if (!is.null(full) && is.finite(full$gcv) && full$gcv < best_gcv) {
            best_gcv <- full$gcv
            best_th <- c(a, b)
          }
        }
      }
    }
    on_edge <- any(abs(best_th - cfg$coarse_lo) < 1e-12 |
                     abs(best_th - cfg$coarse_hi) < 1e-12) ||
      any(best_th < cfg$coarse_lo - 1e-12 | best_th > cfg$coarse_hi + 1e-12)
  }
  flat_n <- sum(TTPsplines:::.tt_lab_flatness_mask(gcv_vals, 0.01))
  interior <- !on_edge &&
    best_th[1] > cfg$coarse_lo && best_th[1] < cfg$coarse_hi &&
    best_th[2] > cfg$coarse_lo && best_th[2] < cfg$coarse_hi

  idx <- idx + 1L
  rows[[idx]] <- data.frame(
    scenario = sc,
    dataset_id = des$dataset_id,
    n = des$n, k = des$k, sigma = des$sigma, seed = des$seed,
    best_log10_l1 = best_th[1],
    best_log10_l2 = best_th[2],
    best_gcv_full = best_gcv,
    min_on_boundary = on_edge,
    interior_min_candidate = interior,
    flat_1pct_count = flat_n,
    flat_1pct_frac = flat_n / max(length(gcv_vals), 1),
    coarse_by = cfg$coarse_by,
    truth_note = des$truth_note,
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "[%s] full-TP min at (%.1f,%.1f) gcv=%.4g interior=%s flat1%%=%d/%d",
    sc, best_th[1], best_th[2], best_gcv, interior, flat_n, length(gcv_vals)
  ))
}

tab <- do.call(rbind, rows)
path <- file.path(cfg$out_dir, "phase1_designs.csv")
.phase1_checkpoint_write(tab, path)
print(tab)
message("Wrote ", path)
