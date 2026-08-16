# v1 GATE: d=4 validation (known scenarios — no new scientific cases).
#
# Extends the d=3 adaptive gate:
#   - strong_aniso: auto-expand only rough margin 2; recover low-θ2 valley;
#     label weakly_identified / effectively_unpenalized
#   - smooth_smooth: finishes; no spurious expand of interior margins
#
# Modes (env TT_GGCV_V1_MODE):
#   SMOKE — smaller n / batches / M (default for first pass)
#   GATE  — full budgets aligned with d=3 aniso gate
#
# Usage (package root):
#   TT_GGCV_V1_MODE=SMOKE Rscript inst/benchmarks/global_gcv/run_v1_gate_d4.R
#   TT_GGCV_V1_MODE=GATE  Rscript inst/benchmarks/global_gcv/run_v1_gate_d4.R
#   TT_GGCV_V1_SCENARIOS=strong_aniso Rscript ...

root <- getwd()
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, export_all = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

mode <- toupper(Sys.getenv("TT_GGCV_V1_MODE", "SMOKE"))
scen_env <- Sys.getenv("TT_GGCV_V1_SCENARIOS", "strong_aniso,smooth_smooth")
scenarios <- trimws(strsplit(scen_env, ",", fixed = TRUE)[[1]])
scenarios <- scenarios[nzchar(scenarios)]

cfg <- if (identical(mode, "GATE")) {
  list(
    n = 250L, n_test = 200L, k = 6L, rank = 5L,
    sobol_batches = c(64L, 64L, 128L),
    n_refine = 5L, n_diverse = 10L,
    M_search = 12L, M_final = 30L, core_starts_final = 3L,
    max_sweeps = 40L,
    stop_on_stable_region = FALSE
  )
} else {
  list(
    n = 180L, n_test = 100L, k = 6L, rank = 4L,
    sobol_batches = c(48L, 48L, 96L),
    n_refine = 4L, n_diverse = 8L,
    M_search = 10L, M_final = 20L, core_starts_final = 2L,
    max_sweeps = 30L,
    stop_on_stable_region = FALSE
  )
}

out_dir <- file.path(root, "inst", "benchmarks", "global_gcv",
                     "results", sprintf("global_opt_v1_d4_%s", tolower(mode)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed_data <- 1L
ctrl <- tt_control(
  max_sweeps = cfg$max_sweeps, tol = 1e-8, compute_edf = FALSE, seed = seed_data
)

run_one <- function(scenario) {
  message(sprintf("=== v1 %s d=4 %s ===", mode, scenario))
  des <- TTPsplines:::.tt_lab_phase1_make_design_dn(
    d = 4L, scenario = scenario,
    n = cfg$n, n_test = cfg$n_test, k = cfg$k, seed = seed_data
  )
  out <- TTPsplines:::tt_global_lambda_optimize_v1(
    y = des$y, X = des$X, rank = cfg$rank, k = cfg$k,
    degree = des$degree, penalty_order = des$penalty_order,
    theta_lower = -5, theta_upper = 5,
    sobol_batches = cfg$sobol_batches,
    n_refine = cfg$n_refine, n_diverse = cfg$n_diverse,
    M_search = cfg$M_search, M_final = cfg$M_final,
    core_starts_final = cfg$core_starts_final,
    adaptive_box = TRUE,
    max_expansions = 2L,
    expand_boundary_tol = 0.5,
    expand_step = 3,
    near_optimal_tol = 0.01,
    profile_lambda = TRUE,
    classify_identifiability = TRUE,
    include_cgcv_anchor = TRUE,
    stop_on_stable_region = isTRUE(cfg$stop_on_stable_region),
    seed = 1L,
    control = ctrl,
    verbose = TRUE
  )

  I <- out$profile_intervals
  labs <- out$lambda_identifiability
  exp_dims <- out$expanded_dimensions
  th <- out$theta
  d <- length(th)

  if (identical(scenario, "strong_aniso")) {
    lab2 <- labs[["lambda2"]]
    t2_hi <- I["theta2", "hi"]
    t2_lo <- I["theta2", "lo"]
    valley_ok <- (is.finite(t2_hi) && t2_hi <= -3.5) ||
      (th[2] <= -4.0) ||
      (is.finite(t2_lo) && t2_lo <= -6 &&
         lab2 %in% c("weakly_identified", "effectively_unpenalized"))
    id_ok <- lab2 %in% c("weakly_identified", "effectively_unpenalized")
    expand_ok <- !length(exp_dims) || (2L %in% exp_dims)
    # Spurious: expanded a clearly interior margin other than 2
    other <- setdiff(exp_dims, 2L)
    spurious <- FALSE
    if (length(other)) {
      for (j in other) {
        if (abs(th[j] - (-5)) > 1 && abs(5 - th[j]) > 1) spurious <- TRUE
      }
    }
    expand_clean <- !isTRUE(spurious)
    checks <- c(
      finite = all(is.finite(th)),
      valley_ok = valley_ok,
      id_ok = id_ok,
      expand_ok = expand_ok,
      expand_clean = expand_clean,
      has_region = out$near_optimal_region$n_points >= 1L
    )
  } else {
    # smooth_smooth: finite + region; no expand of clearly interior margins
    spurious <- FALSE
    if (length(exp_dims)) {
      for (j in exp_dims) {
        if (abs(th[j] - (-5)) > 1.25 && abs(5 - th[j]) > 1.25) spurious <- TRUE
      }
    }
    checks <- c(
      finite = all(is.finite(th)),
      has_region = out$near_optimal_region$n_points >= 1L,
      expand_clean = !isTRUE(spurious),
      outer_stable = isTRUE(out$outer_search_stable)
    )
    lab2 <- labs[["lambda2"]]
    t2_hi <- I["theta2", "hi"]
    t2_lo <- I["theta2", "lo"]
  }

  pass <- all(checks)
  row <- data.frame(
    mode = mode,
    scenario = scenario,
    pass = pass,
    d = d,
    theta = paste(sprintf("%.3f", th), collapse = ","),
    t2_lo = t2_lo, t2_hi = t2_hi,
    lambda2_id = lab2 %||% NA_character_,
    expanded = paste(exp_dims, collapse = ","),
    n_expansions = out$diagnostics$n_expansions,
    winner_source = out$winner_source,
    outer_stable = out$outer_search_stable,
    n_region = out$near_optimal_region$n_points,
    batches = out$diagnostics$sobol_batches_used,
    n_explore = nrow(out$explore),
    elapsed = out$elapsed,
    stringsAsFactors = FALSE
  )

  scen_dir <- file.path(out_dir, scenario)
  dir.create(scen_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(row, file.path(scen_dir, "summary.csv"), row.names = FALSE)
  utils::write.csv(out$explore, file.path(scen_dir, "explore.csv"),
                   row.names = FALSE)
  utils::write.csv(out$final_candidates, file.path(scen_dir, "final.csv"),
                   row.names = FALSE)
  writeLines(
    c(if (pass) "PASS" else "FAIL",
      paste(names(checks), checks, sep = "=", collapse = " | "),
      paste(capture.output(print(row)), collapse = "\n")),
    file.path(scen_dir, "VERDICT.txt")
  )
  list(pass = pass, checks = checks, row = row, out = out)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

rows <- list()
all_pass <- TRUE
for (sc in scenarios) {
  res <- run_one(sc)
  rows[[length(rows) + 1L]] <- res$row
  all_pass <- all_pass && isTRUE(res$pass)
  cat(sprintf("\n---- %s: %s ----\n", sc, if (res$pass) "PASS" else "FAIL"))
  print(res$checks)
  print(res$row)
}

sum_tab <- do.call(rbind, rows)
utils::write.csv(sum_tab, file.path(out_dir, "summary_all.csv"), row.names = FALSE)
writeLines(
  c(if (all_pass) "PASS" else "FAIL",
    paste(sum_tab$scenario, sum_tab$pass, sep = "=", collapse = " | ")),
  file.path(out_dir, "VERDICT.txt")
)

cat("\n==== V1 D4", mode, "GATE:", if (all_pass) "PASS" else "FAIL", "====\n")
print(sum_tab)
message("Wrote ", out_dir)
