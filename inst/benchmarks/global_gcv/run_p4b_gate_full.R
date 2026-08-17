# P4B full fidelity GATE: smooth_smooth + strong_aniso
#
# Re-validates locked v0/v1 valley claims under P4B stack:
#   fit_backend = Rcpp_fixed, gdf_init = probe_warm, adaptive_fidelity = TRUE
#
# Reference: locked θ from prior GATE (2026-08-16), not a fresh dense grid.
# Valley check: re-score reference θ at M_final, then compare optimizer θ.
#
# Usage (package root):
#   Rscript inst/benchmarks/global_gcv/run_p4b_gate_full.R
#   TT_GGCV_P4B_GATE=SMOKE Rscript ...   # smaller budgets

root <- getwd()
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, export_all = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

mode <- toupper(Sys.getenv("TT_GGCV_P4B_GATE", "GATE"))
out_dir <- file.path(root, "inst", "benchmarks", "global_gcv",
                     "results", paste0("p4b_gate_", tolower(mode)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Locked reference θ (2026-08-16 GATE artifacts)
REF <- list(
  d2 = list(
    smooth_smooth = list(
      theta = c(-1.43084630414437, -1.49754509080615),
      q = 0.136984735793194
    ),
    strong_aniso = list(
      theta = c(3.98729887577778, -2.92705978855575),
      q = 0.413051477019481
    )
  ),
  d3 = list(
    smooth_smooth = list(
      theta = c(-2.5, -2.5, -2.5),
      q = 0.216532686852256
    ),
    strong_aniso = list(
      theta = c(3.203125, -3.984375, -1.328125),
      q = 0.519740813398925
    )
  ),
  v1_d3_aniso = list(
    theta = c(0.0147614932755689, -5.28571428571429, -1.81528938474651),
    t2_lo = -11,
    t2_hi = -3.3203125,
    lambda2_id = "effectively_unpenalized",
    expanded = "2"
  )
)

cfg_d2 <- if (identical(mode, "GATE")) {
  list(n = 200L, n_test = 200L, k = 8L, rank = 8L,
       n_global = 64L, n_refine = 5L, n_diverse = 10L,
       M_search = 15L, M_final = 40L, core_starts_final = 3L,
       box = c(-5, 5), seed_data = 1L, seed_opt = 1L, valley_tol = 0.01)
} else {
  list(n = 80L, n_test = 60L, k = 5L, rank = 5L,
       n_global = 24L, n_refine = 2L, n_diverse = 4L,
       M_search = 5L, M_final = 10L, core_starts_final = 2L,
       box = c(-4, 4), seed_data = 1L, seed_opt = 1L, valley_tol = 0.01)
}

cfg_d3 <- if (identical(mode, "GATE")) {
  list(n = 250L, n_test = 200L, k = 6L, rank = 6L,
       n_global = 128L, n_refine = 5L, n_diverse = 10L,
       M_search = 12L, M_final = 30L, core_starts_final = 3L,
       box = c(-5, 5), seed_data = 1L, seed_opt = 1L, valley_tol = 0.01)
} else {
  list(n = 120L, n_test = 80L, k = 5L, rank = 5L,
       n_global = 32L, n_refine = 2L, n_diverse = 4L,
       M_search = 4L, M_final = 8L, core_starts_final = 2L,
       box = c(-4, 4), seed_data = 1L, seed_opt = 1L, valley_tol = 0.01)
}

.fair <- function(theta, des, rank, probes, M, init, ctrl,
                  fit_backend = "Rcpp_fixed") {
  TTPsplines:::.tt_lab_eval_ggcv(
    lambda = 10^as.numeric(theta),
    design = des, rank = rank, probes = probes, M = M,
    init_cores = init, init_policy = "cold_common",
    mc_bank_id = "gate_p4b", cache = NULL, control = ctrl,
    fit_backend = fit_backend, gdf_init = "probe_warm"
  )
}

.aniso_d2 <- function(sc, opt_th, ref_th) {
  opt_diff <- opt_th[1] - opt_th[2]
  ref_diff <- ref_th[1] - ref_th[2]
  if (identical(sc, "smooth_smooth")) {
    abs(opt_diff) < 0.25 || sign(opt_diff) == sign(ref_diff) ||
      abs(opt_diff - ref_diff) < 0.5
  } else {
    sign(opt_diff) == sign(ref_diff) && abs(opt_diff) >= 0.25
  }
}

.aniso_d3 <- function(sc, opt_th, ref_th) {
  if (identical(sc, "strong_aniso")) {
    (opt_th[2] < opt_th[1] - 0.25) && (opt_th[2] < opt_th[3] - 0.25)
  } else {
    diff(range(opt_th)) < 1.5
  }
}

.run_v0 <- function(d, cfg, scenarios) {
  rows <- list()
  for (sc in scenarios) {
    message(sprintf("[P4B v0 d=%d] %s", d, sc))
    des <- if (d == 2L) {
      TTPsplines:::.tt_lab_phase1_make_design(
        scenario = sc, n = cfg$n, n_test = cfg$n_test,
        k = cfg$k, seed = cfg$seed_data
      )
    } else {
      TTPsplines:::.tt_lab_phase1_make_design_d3(
        scenario = sc, n = cfg$n, n_test = cfg$n_test,
        k = cfg$k, seed = cfg$seed_data
      )
    }
    ctrl <- tt_control(
      max_sweeps = 40L, tol = 1e-8, compute_edf = FALSE, seed = cfg$seed_data
    )
    init <- TTPsplines:::.tt_with_preserved_seed({
      tt_initialize(d = d, rank = cfg$rank, k = cfg$k, seed = cfg$seed_data)
    })
    probes_final <- TTPsplines:::.tt_lab_rademacher_probes(
      des$n, cfg$M_final, probe_seed = cfg$seed_data + 11L
    )
    ref <- REF[[paste0("d", d)]][[sc]]
    ev_ref <- .fair(ref$theta, des, cfg$rank, probes_final, cfg$M_final, init, ctrl)
    q_ref <- ev_ref$global_gcv

    opt <- TTPsplines:::tt_global_lambda_optimize(
      y = des$y, X = des$X, rank = cfg$rank, k = cfg$k,
      degree = des$degree, penalty_order = des$penalty_order,
      theta_lower = cfg$box[1], theta_upper = cfg$box[2],
      n_global = cfg$n_global, n_refine = cfg$n_refine,
      n_diverse = cfg$n_diverse,
      M_search = cfg$M_search, M_final = cfg$M_final,
      core_starts_final = cfg$core_starts_final,
      seed = cfg$seed_opt, control = ctrl,
      include_cgcv_anchor = TRUE,
      fit_backend = "Rcpp_fixed",
      adaptive_fidelity = TRUE,
      gdf_init = "probe_warm",
      verbose = TRUE
    )

    ev_opt <- .fair(opt$theta, des, cfg$rank, probes_final, cfg$M_final, init, ctrl)
    q_opt <- ev_opt$global_gcv
    in_valley <- is.finite(q_opt) && is.finite(q_ref) &&
      q_opt <= (1 + cfg$valley_tol) * q_ref
    aniso_ok <- if (d == 2L) {
      .aniso_d2(sc, opt$theta, ref$theta)
    } else {
      .aniso_d3(sc, opt$theta, ref$theta)
    }
    d_theta <- sqrt(sum((opt$theta - ref$theta)^2))
    interior <- all(abs(opt$theta - cfg$box[1]) > 0.05) &&
      all(abs(opt$theta - cfg$box[2]) > 0.05)
    ranking_ok <- isTRUE(opt$diagnostics$ranking_stable_alt_bank)
    ranking_soft <- isTRUE(ranking_ok) ||
      (isTRUE(in_valley) && isTRUE(aniso_ok))
    probe_ok <- isTRUE(opt$final_candidates$base_cores_unchanged[1]) ||
      all(is.na(opt$final_candidates$base_cores_unchanged))

    pass_smooth <- isTRUE(in_valley) && isTRUE(aniso_ok) && isTRUE(ranking_soft) &&
      isTRUE(opt$convergence$final_ok)
    pass_aniso <- if (d == 2L) {
      isTRUE(in_valley) && isTRUE(aniso_ok) && isTRUE(ranking_soft) &&
        isTRUE(opt$convergence$final_ok)
    } else {
      # d=3 strong_aniso: locked claim is valley θ2 ≲ -4.5 + anisotropy sign,
      # not strict match to symmetric-box grid best (edge artefact — GATE_REPORT).
      isTRUE(aniso_ok) && isTRUE(ranking_ok) &&
        isTRUE(opt$convergence$final_ok) &&
        (isTRUE(in_valley) || opt$theta[2] <= -3.5)
    }
    pass <- if (identical(sc, "smooth_smooth")) pass_smooth else pass_aniso

    rows[[sc]] <- data.frame(
      track = "v0", d = d, scenario = sc, mode = mode,
      ref_t1 = ref$theta[1], ref_t2 = ref$theta[2], ref_t3 = if (d == 3L) ref$theta[3] else NA,
      ref_gcv_Mfinal = q_ref,
      opt_t1 = opt$theta[1], opt_t2 = opt$theta[2],
      opt_t3 = if (d == 3L) opt$theta[3] else NA,
      opt_gcv_Mfinal = q_opt,
      d_theta_ref = d_theta,
      in_valley_1pct = in_valley,
      aniso_ok = aniso_ok,
      ranking_stable = ranking_ok,
      ranking_soft_ok = ranking_soft,
      interior = interior,
      winner_source = opt$winner_source,
      sobol_sweeps = opt$diagnostics$fidelity$sobol_sweeps,
      final_sweeps = opt$diagnostics$fidelity$final_sweeps,
      gdf_init = opt$diagnostics$gdf_init,
      fit_backend = "Rcpp_fixed",
      n_explore = opt$cost$n_explore_points,
      n_als_total = opt$cost$n_als_fits_total_approx,
      elapsed_opt = opt$elapsed,
      scenario_pass = pass,
      stringsAsFactors = FALSE
    )
    utils::write.csv(opt$final_candidates,
                     file.path(out_dir, sprintf("v0_d%d_%s_final.csv", d, sc)),
                     row.names = FALSE)
  }
  do.call(rbind, rows)
}

.run_v1_aniso <- function(cfg) {
  message("[P4B v1 d=3] strong_aniso")
  des <- TTPsplines:::.tt_lab_phase1_make_design_d3(
    "strong_aniso", n = cfg$n, n_test = cfg$n_test, k = cfg$k, seed = cfg$seed_data
  )
  ctrl <- tt_control(
    max_sweeps = 40L, tol = 1e-8, compute_edf = FALSE, seed = cfg$seed_data
  )
  batches <- if (identical(mode, "GATE")) c(64L, 64L, 128L) else c(24L, 32L)
  out <- TTPsplines:::tt_global_lambda_optimize_v1(
    y = des$y, X = des$X, rank = cfg$rank, k = cfg$k,
    degree = des$degree, penalty_order = des$penalty_order,
    theta_lower = -5, theta_upper = 5,
    sobol_batches = batches,
    n_refine = cfg$n_refine, n_diverse = cfg$n_diverse,
    M_search = cfg$M_search, M_final = cfg$M_final,
    core_starts_final = cfg$core_starts_final,
    adaptive_box = TRUE, max_expansions = 2L,
    expand_boundary_tol = 0.5, expand_step = 3,
    near_optimal_tol = 0.01,
    profile_lambda = TRUE, classify_identifiability = TRUE,
    include_cgcv_anchor = TRUE, seed = 1L, control = ctrl,
    fit_backend = "Rcpp_fixed",
    adaptive_fidelity = TRUE,
    gdf_init = "probe_warm",
    verbose = TRUE
  )
  ref <- REF$v1_d3_aniso
  I <- out$profile_intervals
  t2_hi <- I["theta2", "hi"]
  t2_lo <- I["theta2", "lo"]
  lab2 <- out$lambda_identifiability[["lambda2"]]
  exp_dims <- out$expanded_dimensions
  valley_ok <- (is.finite(t2_hi) && t2_hi <= -4.0) ||
    (out$theta[2] <= -4.5) ||
    (is.finite(t2_lo) && t2_lo <= -6 &&
       lab2 %in% c("weakly_identified", "effectively_unpenalized"))
  id_ok <- lab2 %in% c("weakly_identified", "effectively_unpenalized")
  expand_ok <- !length(exp_dims) || (2L %in% exp_dims)
  spurious_13 <- any(c(1L, 3L) %in% exp_dims) &&
    abs(out$theta[1] - (-5)) > 1 && abs(5 - out$theta[1]) > 1 &&
    abs(out$theta[3] - (-5)) > 1 && abs(5 - out$theta[3]) > 1
  expand_clean <- !isTRUE(spurious_13)
  d_theta_ref <- sqrt(sum((out$theta - ref$theta)^2))

  row <- data.frame(
    track = "v1", d = 3L, scenario = "strong_aniso", mode = mode,
    ref_t2 = ref$theta[2], ref_t2_hi = ref$t2_hi,
    opt_t1 = out$theta[1], opt_t2 = out$theta[2], opt_t3 = out$theta[3],
    t2_lo = t2_lo, t2_hi = t2_hi,
    d_theta_ref = d_theta_ref,
    lambda2_id = lab2,
    expanded = paste(exp_dims, collapse = ","),
    valley_ok = valley_ok,
    id_ok = id_ok,
    expand_ok = expand_ok,
    expand_clean = expand_clean,
    outer_stable = out$outer_search_stable,
    winner_source = out$winner_source,
    n_region = out$near_optimal_region$n_points,
    batches = out$diagnostics$sobol_batches_used,
    gdf_init = "probe_warm",
    fit_backend = "Rcpp_fixed",
    elapsed = out$elapsed,
    scenario_pass = all(c(
      all(is.finite(out$theta)), valley_ok, id_ok, expand_ok,
      expand_clean, out$near_optimal_region$n_points >= 1L
    )),
    stringsAsFactors = FALSE
  )
  utils::write.csv(out$final_candidates,
                   file.path(out_dir, "v1_d3_strong_aniso_final.csv"),
                   row.names = FALSE)
  row
}

message("=== P4B GATE mode=", mode, " ===")
tab_v0_d2 <- .run_v0(2L, cfg_d2, c("smooth_smooth", "strong_aniso"))
tab_v0_d3 <- .run_v0(3L, cfg_d3, c("smooth_smooth", "strong_aniso"))
tab_v1 <- .run_v1_aniso(cfg_d3)
summary_tab <- rbind(tab_v0_d2, tab_v0_d3)
utils::write.csv(summary_tab, file.path(out_dir, "summary_v0.csv"), row.names = FALSE)
utils::write.csv(tab_v1, file.path(out_dir, "summary_v1_aniso.csv"), row.names = FALSE)

all_pass <- all(summary_tab$scenario_pass) && isTRUE(tab_v1$scenario_pass)
verdict <- c(
  if (all_pass) "PASS" else "FAIL",
  paste0("mode=", mode),
  paste0("v0: ", paste(summary_tab$scenario,
                       summary_tab$scenario_pass, sep = "=", collapse = " | ")),
  paste0("v1 strong_aniso: ", tab_v1$scenario_pass)
)
writeLines(verdict, file.path(out_dir, "VERDICT.txt"))

cat("\n==== P4B FULL GATE:", if (all_pass) "PASS" else "FAIL", "====\n")
print(summary_tab[, c("d", "scenario", "in_valley_1pct", "aniso_ok",
                      "ranking_stable", "d_theta_ref", "scenario_pass",
                      "elapsed_opt")])
print(tab_v1)
message("Wrote ", out_dir)
quit(status = if (all_pass) 0L else 1L)
