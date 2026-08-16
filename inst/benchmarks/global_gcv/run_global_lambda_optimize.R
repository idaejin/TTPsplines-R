# GATE driver for tt_global_lambda_optimize (experimental).
#
# Primary question: does Sobol+refine reach the same 1% GCV valley as a dense
# TT-gGCV grid, with fewer evaluations, correct anisotropy, without relying
# only on distance in theta or on the cGCV anchor?
#
# Usage (package root):
#   TT_GGCV_OPT_MODE=GATE   Rscript inst/benchmarks/global_gcv/run_global_lambda_optimize.R
#   TT_GGCV_OPT_MODE=SEEDS  Rscript inst/benchmarks/global_gcv/run_global_lambda_optimize.R
#   TT_GGCV_OPT_MODE=ABLATION Rscript ...
#   TT_GGCV_OPT_MODE=SMOKE  Rscript ...

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

mode <- toupper(Sys.getenv("TT_GGCV_OPT_MODE", "SMOKE"))
out_dir <- file.path(root, "inst", "benchmarks", "global_gcv",
                     "results", paste0("global_opt_", tolower(mode)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

.cfg_for_mode <- function(mode) {
  if (identical(mode, "GATE") || identical(mode, "SEEDS") ||
      identical(mode, "ABLATION")) {
    list(
      scenarios = c("smooth_smooth", "strong_aniso"),
      n = 200L, n_test = 200L, k = 8L, rank = 8L,
      n_global = 64L, n_refine = 5L, n_diverse = 10L,
      M_search = 15L, M_final = 40L,
      core_starts_final = 3L,
      grid_by = 1.0, box = c(-5, 5),
      seed_data = 1L,
      seed_opt = 1L,
      valley_tol = 0.01,
      seeds_opt = c(1L, 11L, 21L, 31L, 41L)
    )
  } else {
    list(
      scenarios = c("smooth_smooth"),
      n = 80L, n_test = 60L, k = 5L, rank = 5L,
      n_global = 24L, n_refine = 3L, n_diverse = 4L,
      M_search = 5L, M_final = 10L,
      core_starts_final = 2L,
      grid_by = 1.5, box = c(-4, 4),
      seed_data = 1L,
      seed_opt = 1L,
      valley_tol = 0.01,
      seeds_opt = 1L
    )
  }
}

.fair_reeval <- function(theta, des, rank, probes, M, common_init, ctrl,
                         epsilon_rel = 1e-3) {
  lam <- 10^as.numeric(theta)
  TTPsplines:::.tt_lab_eval_ggcv(
    lambda = lam, design = des, rank = rank, probes = probes,
    M = M, init_cores = common_init, init_policy = "cold_common",
    mc_bank_id = "gate_final", cache = NULL, control = ctrl,
    epsilon_rel = epsilon_rel
  )
}

.gate_metrics <- function(opt_theta, grid_theta, q_opt, q_grid,
                          q_alt_opt = NA_real_,
                          boundary, box, valley_tol = 0.01,
                          scenario, ranking_stable, winner_source,
                          n_explore_points, n_grid,
                          n_als_opt, n_als_grid,
                          final_ok, refine_nlminb,
                          grid_on_box_edge = FALSE) {
  d_theta <- sqrt(sum((opt_theta - grid_theta)^2))
  delta_rel <- if (is.finite(q_grid) && abs(q_grid) > 1e-15) {
    (q_opt - q_grid) / q_grid
  } else {
    NA_real_
  }
  in_valley <- is.finite(q_opt) && is.finite(q_grid) &&
    q_opt <= (1 + valley_tol) * q_grid
  alt_in_valley <- if (is.finite(q_alt_opt) && is.finite(q_grid)) {
    q_alt_opt <= (1 + valley_tol) * q_grid
  } else {
    NA
  }

  opt_diff <- opt_theta[1] - opt_theta[2]
  grid_diff <- grid_theta[1] - grid_theta[2]
  sign_opt <- sign(opt_diff)
  sign_grid <- sign(grid_diff)
  iso_opt <- abs(opt_diff) < 0.25
  iso_grid <- abs(grid_diff) < 0.25
  aniso_ok <- if (identical(scenario, "smooth_smooth")) {
    iso_opt || (sign_opt == sign_grid) || abs(opt_diff - grid_diff) < 0.5
  } else if (identical(scenario, "strong_aniso")) {
    if (iso_grid) TRUE else (sign_opt == sign_grid && abs(opt_diff) >= 0.25)
  } else {
    sign_opt == sign_grid || (iso_opt && iso_grid)
  }

  false_boundary <- isTRUE(boundary) &&
    all(abs(opt_theta - box[1]) > 0.05) &&
    all(abs(opt_theta - box[2]) > 0.05)

  # Locked soft ranking rule: flip is not FAIL if both candidates stay in
  # the 1% valley (and aniso already OK for the reported winner).
  ranking_strict <- isTRUE(ranking_stable)
  ranking_soft_ok <- ranking_strict ||
    (isTRUE(in_valley) && isTRUE(alt_in_valley) && isTRUE(aniso_ok))

  fewer_explore_points <- n_explore_points < n_grid
  fewer_als_total <- if (is.finite(n_als_opt) && is.finite(n_als_grid)) {
    n_als_opt < n_als_grid
  } else {
    NA
  }

  nlminb_ok <- if (!length(refine_nlminb)) {
    TRUE
  } else {
    mean(refine_nlminb == 0L, na.rm = TRUE) >= 0.5
  }

  list(
    d_theta = d_theta,
    delta_rel = delta_rel,
    in_valley_1pct = in_valley,
    alt_in_valley_1pct = alt_in_valley,
    opt_log_ratio = opt_diff,
    grid_log_ratio = grid_diff,
    aniso_ok = aniso_ok,
    iso_opt = iso_opt,
    ranking_stable = ranking_strict,
    ranking_soft_ok = ranking_soft_ok,
    false_boundary = false_boundary,
    fewer_explore_points = fewer_explore_points,
    fewer_als_total = fewer_als_total,
    n_explore_points = n_explore_points,
    n_grid = n_grid,
    n_als_opt = n_als_opt,
    n_als_grid = n_als_grid,
    grid_on_box_edge = isTRUE(grid_on_box_edge),
    final_ok = isTRUE(final_ok),
    nlminb_ok = nlminb_ok,
    winner_source = winner_source
  )
}

# Strict scenario pass uses soft ranking rule (locked).
# "Fewer ALS total" is reported but not required for PASS (may go either way).
.scenario_pass <- function(m, scenario) {
  isTRUE(m$in_valley_1pct) &&
    isTRUE(m$aniso_ok) &&
    isTRUE(m$ranking_soft_ok) &&
    !isTRUE(m$false_boundary) &&
    isTRUE(m$fewer_explore_points) &&
    isTRUE(m$final_ok) &&
    isTRUE(m$nlminb_ok)
}

.run_one <- function(sc, cfg, seed_opt, ablation = "full", verbose = TRUE) {
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    scenario = sc, n = cfg$n, n_test = cfg$n_test, k = cfg$k,
    seed = cfg$seed_data
  )
  ctrl <- tt_control(
    max_sweeps = 40L, tol = 1e-8, compute_edf = FALSE, seed = cfg$seed_data
  )
  common_init <- TTPsplines:::.tt_with_preserved_seed({
    tt_initialize(d = 2L, rank = cfg$rank, k = cfg$k, seed = cfg$seed_data)
  })
  # Fixed MC banks for fair comparison (dataset-fixed; independent of Sobol seed)
  probes_search <- TTPsplines:::.tt_lab_rademacher_probes(
    des$n, cfg$M_search, probe_seed = cfg$seed_data + 11L
  )
  probes_final <- TTPsplines:::.tt_lab_rademacher_probes(
    des$n, cfg$M_final, probe_seed = cfg$seed_data + 11L
  )

  # Reference grid at M_search (scout), then re-score best + opt at M_final
  g <- seq(cfg$box[1], cfg$box[2], by = cfg$grid_by)
  coords <- expand.grid(t1 = g, t2 = g, KEEP.OUT.ATTRS = FALSE)
  cache <- new.env(parent = emptyenv())
  grid_rows <- vector("list", nrow(coords))
  for (i in seq_len(nrow(coords))) {
    lam <- 10^c(coords$t1[i], coords$t2[i])
    ev <- TTPsplines:::.tt_lab_eval_ggcv(
      lambda = lam, design = des, rank = cfg$rank, probes = probes_search,
      M = cfg$M_search, init_cores = common_init, init_policy = "cold_common",
      mc_bank_id = "bank0", cache = cache, control = ctrl
    )
    grid_rows[[i]] <- data.frame(
      t1 = coords$t1[i], t2 = coords$t2[i],
      global_gcv = ev$global_gcv, valid = ev$valid,
      stringsAsFactors = FALSE
    )
  }
  grid <- do.call(rbind, grid_rows)
  i_grid <- which.min(grid$global_gcv)
  grid_theta <- c(grid$t1[i_grid], grid$t2[i_grid])

  # Ablation knobs
  n_global <- cfg$n_global
  n_refine <- cfg$n_refine
  include_cgcv <- TRUE
  if (identical(ablation, "no_cgcv")) {
    include_cgcv <- FALSE
  } else if (identical(ablation, "anchors_only")) {
    n_global <- 0L
  } else if (identical(ablation, "no_refine")) {
    n_refine <- 0L
  }

  opt <- TTPsplines:::tt_global_lambda_optimize(
    y = des$y, X = des$X, rank = cfg$rank, k = cfg$k,
    degree = des$degree, penalty_order = des$penalty_order,
    theta_lower = cfg$box[1], theta_upper = cfg$box[2],
    n_global = n_global, n_refine = n_refine, n_diverse = cfg$n_diverse,
    M_search = cfg$M_search, M_final = cfg$M_final,
    core_starts_final = cfg$core_starts_final,
    seed = as.integer(seed_opt), control = ctrl,
    include_cgcv_anchor = include_cgcv,
    verbose = verbose
  )

  # Fair M_final re-evaluation of grid best, optimizer, and alt-bank winner
  ev_grid_f <- .fair_reeval(
    grid_theta, des, cfg$rank, probes_final, cfg$M_final, common_init, ctrl
  )
  ev_opt_f <- .fair_reeval(
    opt$theta, des, cfg$rank, probes_final, cfg$M_final, common_init, ctrl
  )
  theta_alt <- opt$diagnostics$alt_winner_theta %||% opt$theta
  ev_alt_f <- .fair_reeval(
    theta_alt, des, cfg$rank, probes_final, cfg$M_final, common_init, ctrl
  )

  n_explore_points <- as.integer(
    opt$cost$n_explore_points %||% nrow(opt$sobol_results)
  )
  n_grid <- nrow(grid)
  # Grid cost: each cell = 1 TT-gGCV at M_search ≈ (1+M_search) ALS fits
  n_als_grid <- as.integer(n_grid * (1L + cfg$M_search))
  n_als_opt <- as.integer(
    opt$cost$n_als_fits_total_approx %||% NA_integer_
  )
  grid_on_edge <- any(abs(grid_theta - cfg$box[1]) < 1e-9 |
                        abs(grid_theta - cfg$box[2]) < 1e-9)

  m <- .gate_metrics(
    opt_theta = opt$theta,
    grid_theta = grid_theta,
    q_opt = ev_opt_f$global_gcv,
    q_grid = ev_grid_f$global_gcv,
    q_alt_opt = ev_alt_f$global_gcv,
    boundary = opt$boundary,
    box = cfg$box,
    valley_tol = cfg$valley_tol,
    scenario = sc,
    ranking_stable = opt$diagnostics$ranking_stable_alt_bank,
    winner_source = opt$winner_source %||% NA_character_,
    n_explore_points = n_explore_points,
    n_grid = n_grid,
    n_als_opt = n_als_opt,
    n_als_grid = n_als_grid,
    final_ok = opt$convergence$final_ok,
    refine_nlminb = opt$convergence$refine_nlminb,
    grid_on_box_edge = grid_on_edge
  )
  sc_pass <- .scenario_pass(m, sc)

  list(
    row = data.frame(
      scenario = sc,
      mode = mode,
      ablation = ablation,
      seed_opt = as.integer(seed_opt),
      grid_t1 = grid_theta[1],
      grid_t2 = grid_theta[2],
      grid_gcv_Msearch = grid$global_gcv[i_grid],
      grid_gcv_Mfinal = ev_grid_f$global_gcv,
      grid_on_box_edge = m$grid_on_box_edge,
      opt_t1 = opt$theta[1],
      opt_t2 = opt$theta[2],
      opt_gcv_reported = opt$gcv,
      opt_gcv_Mfinal = ev_opt_f$global_gcv,
      alt_gcv_Mfinal = ev_alt_f$global_gcv,
      d_theta = m$d_theta,
      delta_rel = m$delta_rel,
      in_valley_1pct = m$in_valley_1pct,
      alt_in_valley_1pct = m$alt_in_valley_1pct,
      opt_log_ratio = m$opt_log_ratio,
      grid_log_ratio = m$grid_log_ratio,
      aniso_ok = m$aniso_ok,
      iso_opt = m$iso_opt,
      ranking_stable = m$ranking_stable,
      ranking_soft_ok = m$ranking_soft_ok,
      boundary = opt$boundary,
      false_boundary = m$false_boundary,
      fewer_explore_points = m$fewer_explore_points,
      fewer_als_total = m$fewer_als_total,
      n_explore_points = n_explore_points,
      n_theta_miss_search = opt$cost$n_theta_miss_search %||% NA_integer_,
      n_als_opt = n_als_opt,
      n_als_grid = n_als_grid,
      n_grid = n_grid,
      winner_source = m$winner_source,
      refine_improved = opt$diagnostics$refine_improved %||% FALSE,
      final_ok = m$final_ok,
      nlminb_ok = m$nlminb_ok,
      scenario_pass = sc_pass,
      elapsed_opt = opt$elapsed,
      stringsAsFactors = FALSE
    ),
    grid = grid,
    opt = opt,
    metrics = m
  )
}

cfg <- .cfg_for_mode(mode)
message("MODE=", mode)

summary_rows <- list()

if (identical(mode, "ABLATION")) {
  ablations <- c("full", "no_cgcv", "anchors_only", "no_refine")
  for (ab in ablations) {
    for (sc in cfg$scenarios) {
      message("=== ablation=", ab, " scenario=", sc, " ===")
      res <- .run_one(sc, cfg, seed_opt = cfg$seed_opt, ablation = ab,
                      verbose = TRUE)
      key <- paste(ab, sc, sep = "|")
      summary_rows[[key]] <- res$row
      utils::write.csv(
        res$opt$sobol_results,
        file.path(out_dir, paste0(ab, "_", sc, "_explore.csv")),
        row.names = FALSE
      )
      utils::write.csv(
        res$opt$final_candidates,
        file.path(out_dir, paste0(ab, "_", sc, "_final.csv")),
        row.names = FALSE
      )
      message(sprintf(
        "[%s|%s] valley=%s aniso=%s d=%.3f drel=%.4f winner=%s PASS=%s",
        ab, sc, res$row$in_valley_1pct, res$row$aniso_ok,
        res$row$d_theta, res$row$delta_rel, res$row$winner_source,
        res$row$scenario_pass
      ))
    }
  }
} else if (identical(mode, "SEEDS")) {
  for (sc in cfg$scenarios) {
    for (s in cfg$seeds_opt) {
      message("=== scenario=", sc, " seed_opt=", s, " ===")
      res <- .run_one(sc, cfg, seed_opt = s, ablation = "full", verbose = TRUE)
      key <- paste(sc, s, sep = "|")
      summary_rows[[key]] <- res$row
      message(sprintf(
        "[%s|seed=%d] valley=%s aniso=%s d=%.3f drel=%.4f winner=%s PASS=%s",
        sc, s, res$row$in_valley_1pct, res$row$aniso_ok,
        res$row$d_theta, res$row$delta_rel, res$row$winner_source,
        res$row$scenario_pass
      ))
    }
  }
} else {
  # SMOKE or GATE
  for (sc in cfg$scenarios) {
    message("=== scenario ", sc, " ===")
    res <- .run_one(sc, cfg, seed_opt = cfg$seed_opt, ablation = "full",
                    verbose = TRUE)
    summary_rows[[sc]] <- res$row
    utils::write.csv(res$grid, file.path(out_dir, paste0(sc, "_grid.csv")),
                     row.names = FALSE)
    utils::write.csv(res$opt$sobol_results,
                     file.path(out_dir, paste0(sc, "_explore.csv")),
                     row.names = FALSE)
    utils::write.csv(res$opt$refined_results,
                     file.path(out_dir, paste0(sc, "_refined.csv")),
                     row.names = FALSE)
    utils::write.csv(res$opt$final_candidates,
                     file.path(out_dir, paste0(sc, "_final.csv")),
                     row.names = FALSE)
    message(sprintf(
      "[%s] grid=(%.2f,%.2f) opt=(%.2f,%.2f) d_theta=%.3f delta_rel=%.4f valley=%s aniso=%s winner=%s PASS=%s",
      sc, res$row$grid_t1, res$row$grid_t2, res$row$opt_t1, res$row$opt_t2,
      res$row$d_theta, res$row$delta_rel, res$row$in_valley_1pct,
      res$row$aniso_ok, res$row$winner_source, res$row$scenario_pass
    ))
  }
}

sum_tab <- do.call(rbind, summary_rows)
rownames(sum_tab) <- NULL
utils::write.csv(sum_tab, file.path(out_dir, "summary.csv"), row.names = FALSE)

# Overall verdict for GATE (seed 1 only — region recovery)
if (identical(mode, "GATE")) {
  ss <- sum_tab[sum_tab$scenario == "smooth_smooth", , drop = FALSE]
  sa <- sum_tab[sum_tab$scenario == "strong_aniso", , drop = FALSE]
  both_valley <- all(sum_tab$in_valley_1pct)
  aniso_ss <- isTRUE(ss$aniso_ok[1])
  aniso_sa <- isTRUE(sa$aniso_ok[1])
  rank_soft <- all(sum_tab$ranking_soft_ok)
  no_false_b <- !any(sum_tab$false_boundary)
  fewer_pts <- all(sum_tab$fewer_explore_points)
  conv_ok <- all(sum_tab$final_ok) && all(sum_tab$nlminb_ok)
  edge_sa <- isTRUE(sa$grid_on_box_edge[1])

  if (both_valley && aniso_ss && aniso_sa && rank_soft && no_false_b &&
      fewer_pts && conv_ok) {
    verdict <- "PASS"
  } else if (isTRUE(ss$in_valley_1pct[1]) && isTRUE(ss$aniso_ok[1]) &&
             (!isTRUE(sa$in_valley_1pct[1]) || !isTRUE(sa$aniso_ok[1]))) {
    verdict <- "CONDITIONAL_PASS"
  } else {
    verdict <- "FAIL"
  }
  note <- sprintf(
    paste0(
      "both_valley=%s aniso_ss=%s aniso_sa=%s ranking_soft=%s false_b=%s ",
      "fewer_explore_pts=%s conv=%s strong_aniso_grid_on_edge=%s ",
      "als_opt=%s als_grid=%s claim=box_near_opt_region_not_interior_global"
    ),
    both_valley, aniso_ss, aniso_sa, rank_soft, any(sum_tab$false_boundary),
    fewer_pts, conv_ok, edge_sa,
    paste(ss$n_als_opt, sa$n_als_opt, sep = ","),
    paste(ss$n_als_grid, sa$n_als_grid, sep = ",")
  )
  cat("\n==== GATE VERDICT:", verdict, "====\n", note, "\n", sep = "")
  writeLines(c(verdict, note), file.path(out_dir, "VERDICT.txt"))
}

if (identical(mode, "SEEDS")) {
  for (sc in unique(sum_tab$scenario)) {
    sub <- sum_tab[sum_tab$scenario == sc, , drop = FALSE]
    n_valley <- sum(sub$in_valley_1pct)
    n_strict <- sum(sub$ranking_stable & sub$scenario_pass)
    n_soft <- sum(sub$ranking_soft_ok & sub$in_valley_1pct & sub$aniso_ok)
    aniso_ok <- all(sub$aniso_ok)
    message(sprintf(
      "SEEDS[%s]: valley %d/%d; soft_pass %d/%d; strict_rank %d/%d; aniso_all=%s",
      sc, n_valley, nrow(sub), n_soft, nrow(sub),
      sum(sub$ranking_stable), nrow(sub), aniso_ok
    ))
  }
  # Global MC-stability conclusion
  n_strict_all <- sum(sum_tab$ranking_stable)
  n_all <- nrow(sum_tab)
  mc_verdict <- if (n_strict_all == n_all) {
    "PASS"
  } else if (n_strict_all >= ceiling(0.8 * n_all)) {
    "CONDITIONAL_PASS"
  } else {
    "FAIL"
  }
  cat(sprintf(
    "\n==== SEEDS MC-STABILITY VERDICT: %s (%d/%d strict ranking) ====\n",
    mc_verdict, n_strict_all, n_all
  ))
  writeLines(
    c(mc_verdict, sprintf("strict_ranking=%d/%d", n_strict_all, n_all)),
    file.path(out_dir, "MC_STABILITY_VERDICT.txt")
  )
}

print(sum_tab)
message("Wrote ", out_dir)
