# GATE driver for d=3 joint global λ optimizer (experimental).
#
# Reference = coarse Cartesian TT-gGCV grid (Sobol's advantage is fewer
# explore points than |grid| = L^d). Valley / anisotropy metrics as in d=2.
#
# Usage (package root):
#   TT_GGCV_OPT_MODE=SMOKE Rscript inst/benchmarks/global_gcv/run_global_lambda_optimize_d3.R
#   TT_GGCV_OPT_MODE=GATE  Rscript inst/benchmarks/global_gcv/run_global_lambda_optimize_d3.R

root <- Sys.getenv("TT_PKG_ROOT", unset = "")
if (!nzchar(root)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
    root <- normalizePath(file.path(dirname(script), "..", "..", ".."))
  } else {
    root <- getwd()
  }
}
setwd(root)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, export_all = FALSE, quiet = TRUE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, export_all = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
mode <- toupper(Sys.getenv("TT_GGCV_OPT_MODE", "SMOKE"))
out_dir <- file.path(root, "inst", "benchmarks", "global_gcv",
                     "results", paste0("global_opt_d3_", tolower(mode)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- if (identical(mode, "GATE")) {
  list(
    scenarios = c("smooth_smooth", "strong_aniso"),
    n = 250L, n_test = 200L, k = 6L, rank = 6L,
    n_global = 128L, n_refine = 5L, n_diverse = 10L,
    M_search = 12L, M_final = 30L,
    core_starts_final = 3L,
    # 7 levels on [-5,5] → 343 cells (Cartesian reference; Sobol << L^d)
    grid_vals = seq(-5, 5, length.out = 7L),
    box = c(-5, 5),
    seed_data = 1L,
    seed_opt = 1L,
    valley_tol = 0.01
  )
} else {
  list(
    scenarios = c("smooth_smooth"),
    n = 120L, n_test = 80L, k = 5L, rank = 5L,
    n_global = 32L, n_refine = 2L, n_diverse = 4L,
    M_search = 4L, M_final = 8L,
    core_starts_final = 2L,
    # 4^3 = 64 > 32 Sobol + anchors for smoke fewer-points check
    grid_vals = seq(-4, 4, length.out = 4L),
    box = c(-4, 4),
    seed_data = 1L,
    seed_opt = 1L,
    valley_tol = 0.01
  )
}

.fair_reeval_d3 <- function(theta, des, rank, probes, M, common_init, ctrl) {
  lam <- 10^as.numeric(theta)
  TTPsplines:::.tt_lab_eval_ggcv(
    lambda = lam, design = des, rank = rank, probes = probes,
    M = M, init_cores = common_init, init_policy = "cold_common",
    mc_bank_id = "gate_final", cache = NULL, control = ctrl
  )
}

.aniso_ok_d3 <- function(scenario, opt_th, grid_th) {
  # strong_aniso: margin 2 should be rougher → smaller θ than 1 and 3
  if (identical(scenario, "strong_aniso")) {
    opt_ok <- (opt_th[2] < opt_th[1] - 0.25) && (opt_th[2] < opt_th[3] - 0.25)
    # also agree with grid ordering if grid itself is anisotropic
    grid_an <- (grid_th[2] < grid_th[1] - 0.25) && (grid_th[2] < grid_th[3] - 0.25)
    if (grid_an) {
      opt_ok || (sign(opt_th[2] - opt_th[1]) == sign(grid_th[2] - grid_th[1]) &&
                   sign(opt_th[2] - opt_th[3]) == sign(grid_th[2] - grid_th[3]))
    } else {
      opt_ok
    }
  } else {
    # smooth_smooth: near-isotropic (range of θ small)
    diff(range(opt_th)) < 1.5
  }
}

message("MODE=", mode, " d=3")
summary_rows <- list()

for (sc in cfg$scenarios) {
  message("=== d3 scenario ", sc, " ===")
  des <- TTPsplines:::.tt_lab_phase1_make_design_d3(
    scenario = sc, n = cfg$n, n_test = cfg$n_test, k = cfg$k,
    seed = cfg$seed_data
  )
  ctrl <- tt_control(
    max_sweeps = 40L, tol = 1e-8, compute_edf = FALSE, seed = cfg$seed_data
  )
  common_init <- TTPsplines:::.tt_with_preserved_seed({
    tt_initialize(d = 3L, rank = cfg$rank, k = cfg$k, seed = cfg$seed_data)
  })
  probes_search <- TTPsplines:::.tt_lab_rademacher_probes(
    des$n, cfg$M_search, probe_seed = cfg$seed_data + 11L
  )
  probes_final <- TTPsplines:::.tt_lab_rademacher_probes(
    des$n, cfg$M_final, probe_seed = cfg$seed_data + 11L
  )

  # Coarse Cartesian reference
  g <- as.numeric(cfg$grid_vals)
  coords <- expand.grid(t1 = g, t2 = g, t3 = g, KEEP.OUT.ATTRS = FALSE)
  cache <- new.env(parent = emptyenv())
  grid_rows <- vector("list", nrow(coords))
  t0g <- proc.time()[["elapsed"]]
  for (i in seq_len(nrow(coords))) {
    th <- c(coords$t1[i], coords$t2[i], coords$t3[i])
    ev <- TTPsplines:::.tt_lab_eval_ggcv(
      lambda = 10^th, design = des, rank = cfg$rank, probes = probes_search,
      M = cfg$M_search, init_cores = common_init, init_policy = "cold_common",
      mc_bank_id = "bank0", cache = cache, control = ctrl
    )
    grid_rows[[i]] <- data.frame(
      t1 = th[1], t2 = th[2], t3 = th[3],
      global_gcv = ev$global_gcv, valid = ev$valid,
      stringsAsFactors = FALSE
    )
  }
  grid_elapsed <- proc.time()[["elapsed"]] - t0g
  grid <- do.call(rbind, grid_rows)
  i_grid <- which.min(grid$global_gcv)
  grid_theta <- c(grid$t1[i_grid], grid$t2[i_grid], grid$t3[i_grid])

  opt <- TTPsplines:::tt_global_lambda_optimize(
    y = des$y, X = des$X, rank = cfg$rank, k = cfg$k,
    degree = des$degree, penalty_order = des$penalty_order,
    theta_lower = cfg$box[1], theta_upper = cfg$box[2],
    n_global = cfg$n_global, n_refine = cfg$n_refine,
    n_diverse = cfg$n_diverse,
    M_search = cfg$M_search, M_final = cfg$M_final,
    core_starts_final = cfg$core_starts_final,
    seed = cfg$seed_opt, control = ctrl,
    include_cgcv_anchor = TRUE, verbose = TRUE
  )

  ev_grid_f <- .fair_reeval_d3(
    grid_theta, des, cfg$rank, probes_final, cfg$M_final, common_init, ctrl
  )
  ev_opt_f <- .fair_reeval_d3(
    opt$theta, des, cfg$rank, probes_final, cfg$M_final, common_init, ctrl
  )
  q_grid <- ev_grid_f$global_gcv
  q_opt <- ev_opt_f$global_gcv
  in_valley <- is.finite(q_opt) && is.finite(q_grid) &&
    q_opt <= (1 + cfg$valley_tol) * q_grid
  delta_rel <- (q_opt - q_grid) / q_grid
  d_theta <- sqrt(sum((opt$theta - grid_theta)^2))
  aniso_ok <- .aniso_ok_d3(sc, opt$theta, grid_theta)
  grid_on_edge <- any(abs(grid_theta - cfg$box[1]) < 1e-9 |
                        abs(grid_theta - cfg$box[2]) < 1e-9)
  n_explore <- as.integer(opt$cost$n_explore_points %||% nrow(opt$sobol_results))
  n_grid <- nrow(grid)
  n_als_grid <- as.integer(n_grid * (1L + cfg$M_search))
  n_als_opt <- as.integer(opt$cost$n_als_fits_total_approx %||% NA_integer_)

  sc_pass <- isTRUE(in_valley) && isTRUE(aniso_ok) &&
    isTRUE(opt$diagnostics$ranking_stable_alt_bank %||%
             opt$diagnostics$ranking_soft_ok %||% TRUE) &&
    n_explore < n_grid &&
    isTRUE(opt$convergence$final_ok)

  # Soft ranking: if unstable, still pass when valley+aniso hold
  if (!isTRUE(opt$diagnostics$ranking_stable_alt_bank)) {
    sc_pass <- isTRUE(in_valley) && isTRUE(aniso_ok) &&
      n_explore < n_grid && isTRUE(opt$convergence$final_ok)
  }

  row <- data.frame(
    scenario = sc, mode = mode, d = 3L,
    grid_t1 = grid_theta[1], grid_t2 = grid_theta[2], grid_t3 = grid_theta[3],
    grid_gcv_Mfinal = q_grid,
    grid_on_box_edge = grid_on_edge,
    opt_t1 = opt$theta[1], opt_t2 = opt$theta[2], opt_t3 = opt$theta[3],
    opt_gcv_Mfinal = q_opt,
    d_theta = d_theta, delta_rel = delta_rel,
    in_valley_1pct = in_valley, aniso_ok = aniso_ok,
    ranking_stable = opt$diagnostics$ranking_stable_alt_bank,
    winner_source = opt$winner_source,
    n_explore_points = n_explore, n_grid = n_grid,
    n_als_opt = n_als_opt, n_als_grid = n_als_grid,
    fewer_explore_points = n_explore < n_grid,
    fewer_als_total = is.finite(n_als_opt) && n_als_opt < n_als_grid,
    scenario_pass = sc_pass,
    elapsed_grid = grid_elapsed, elapsed_opt = opt$elapsed,
    stringsAsFactors = FALSE
  )
  summary_rows[[sc]] <- row
  utils::write.csv(grid, file.path(out_dir, paste0(sc, "_grid.csv")),
                   row.names = FALSE)
  utils::write.csv(opt$sobol_results,
                   file.path(out_dir, paste0(sc, "_explore.csv")),
                   row.names = FALSE)
  utils::write.csv(opt$final_candidates,
                   file.path(out_dir, paste0(sc, "_final.csv")),
                   row.names = FALSE)
  message(sprintf(
    "[d3|%s] grid=(%.1f,%.1f,%.1f) opt=(%.2f,%.2f,%.2f) d=%.3f drel=%.4f valley=%s aniso=%s edge=%s explore=%d/%d PASS=%s",
    sc, grid_theta[1], grid_theta[2], grid_theta[3],
    opt$theta[1], opt$theta[2], opt$theta[3],
    d_theta, delta_rel, in_valley, aniso_ok, grid_on_edge,
    n_explore, n_grid, sc_pass
  ))
}

sum_tab <- do.call(rbind, summary_rows)
rownames(sum_tab) <- NULL
utils::write.csv(sum_tab, file.path(out_dir, "summary.csv"), row.names = FALSE)

if (identical(mode, "GATE")) {
  verdict <- if (all(sum_tab$scenario_pass)) {
    "PASS"
  } else if (any(sum_tab$scenario == "smooth_smooth" & sum_tab$scenario_pass) &&
             !all(sum_tab$scenario_pass)) {
    "CONDITIONAL_PASS"
  } else {
    "FAIL"
  }
  cat("\n==== D3 GATE VERDICT:", verdict, "====\n")
  writeLines(verdict, file.path(out_dir, "VERDICT.txt"))
}

print(sum_tab)
message("Wrote ", out_dir)
