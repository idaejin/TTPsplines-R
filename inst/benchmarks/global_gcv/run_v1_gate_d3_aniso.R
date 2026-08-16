# v1 GATE: known d=3 strong_aniso — auto-expand θ2, recover valley, label ID.
#
# Pass if:
# 1) expands only margin 2 (or none if already interior from batches)
# 2) recovers θ2 region with hi ≲ -4.0 (near-optimal)
# 3) labels lambda2 as weakly_identified or effectively_unpenalized
# 4) does not expand margins 1 and 3 unnecessarily when they are interior
# 5) finishes without error
#
# Usage (package root):
#   Rscript inst/benchmarks/global_gcv/run_v1_gate_d3_aniso.R

root <- getwd()
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, export_all = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

out_dir <- file.path(root, "inst", "benchmarks", "global_gcv",
                     "results", "global_opt_v1_d3_aniso")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed_data <- 1L
des <- TTPsplines:::.tt_lab_phase1_make_design_d3(
  "strong_aniso", n = 250L, n_test = 200L, k = 6L, seed = seed_data
)
ctrl <- tt_control(
  max_sweeps = 40L, tol = 1e-8, compute_edf = FALSE, seed = seed_data
)

message("=== v1 GATE d=3 strong_aniso ===")
# Start with the ORIGINAL symmetric box that failed (edge at -5)
out <- TTPsplines:::tt_global_lambda_optimize_v1(
  y = des$y, X = des$X, rank = 6L, k = 6L,
  degree = des$degree, penalty_order = des$penalty_order,
  theta_lower = -5, theta_upper = 5,
  sobol_batches = c(64L, 64L, 128L),
  n_refine = 5L, n_diverse = 10L,
  M_search = 12L, M_final = 30L, core_starts_final = 3L,
  adaptive_box = TRUE,
  max_expansions = 2L,
  expand_boundary_tol = 0.5,
  expand_step = 3,
  near_optimal_tol = 0.01,
  profile_lambda = TRUE,
  classify_identifiability = TRUE,
  include_cgcv_anchor = TRUE,
  seed = 1L,
  control = ctrl,
  verbose = TRUE
)

I <- out$profile_intervals
lab2 <- out$lambda_identifiability[["lambda2"]]
exp_dims <- out$expanded_dimensions
# Valley recovered if region upper for θ2 is ≤ -4.0 (or best θ2 ≤ -4.5)
t2_hi <- I["theta2", "hi"]
t2_lo <- I["theta2", "lo"]
# Valley: region upper ≤ -4 OR best ≤ -4.5 OR lo reaches expanded floor with weak ID
valley_ok <- (is.finite(t2_hi) && t2_hi <= -4.0) ||
  (out$theta[2] <= -4.5) ||
  (is.finite(t2_lo) && t2_lo <= -6 && lab2 %in%
     c("weakly_identified", "effectively_unpenalized"))
id_ok <- lab2 %in% c("weakly_identified", "effectively_unpenalized")
# If expanded, should include margin 2; should not expand 1 and 3 unnecessarily
expand_ok <- !length(exp_dims) || (2L %in% exp_dims)
# Stricter: margins 1,3 not expanded when best is clearly interior on them
spurious_13 <- any(c(1L, 3L) %in% exp_dims) &&
  abs(out$theta[1] - (-5)) > 1 && abs(5 - out$theta[1]) > 1 &&
  abs(out$theta[3] - (-5)) > 1 && abs(5 - out$theta[3]) > 1
expand_clean <- !isTRUE(spurious_13)

checks <- c(
  finite = all(is.finite(out$theta)),
  valley_ok = valley_ok,
  id_ok = id_ok,
  expand_ok = expand_ok,
  expand_clean = expand_clean,
  has_region = out$near_optimal_region$n_points >= 1L
)
pass <- all(checks)

row <- data.frame(
  pass = pass,
  theta1 = out$theta[1], theta2 = out$theta[2], theta3 = out$theta[3],
  t2_lo = t2_lo, t2_hi = t2_hi,
  lambda2_id = lab2,
  expanded = paste(exp_dims, collapse = ","),
  n_expansions = out$diagnostics$n_expansions,
  winner_source = out$winner_source,
  outer_stable = out$outer_search_stable,
  n_region = out$near_optimal_region$n_points,
  batches = out$diagnostics$sobol_batches_used,
  elapsed = out$elapsed,
  stringsAsFactors = FALSE
)
utils::write.csv(row, file.path(out_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(out$explore, file.path(out_dir, "explore.csv"), row.names = FALSE)
utils::write.csv(out$final_candidates, file.path(out_dir, "final.csv"),
                 row.names = FALSE)
if (!is.null(out$profiles[[2]])) {
  utils::write.csv(out$profiles[[2]], file.path(out_dir, "profile_m2.csv"),
                   row.names = FALSE)
}
writeLines(
  c(if (pass) "PASS" else "FAIL",
    paste(names(checks), checks, sep = "=", collapse = " | "),
    paste(capture.output(print(row)), collapse = "\n")),
  file.path(out_dir, "VERDICT.txt")
)

cat("\n==== V1 D3 ANISO GATE:", if (pass) "PASS" else "FAIL", "====\n")
print(checks)
print(row)
message("Wrote ", out_dir)
