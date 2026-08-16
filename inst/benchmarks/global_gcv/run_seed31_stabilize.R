# Stabilize strong_aniso seed_opt=31 under soft ranking rule.
# Dataset + MC bank fixed; only Sobol seed changes (31).
#
# Usage (package root):
#   Rscript inst/benchmarks/global_gcv/run_seed31_stabilize.R

root <- getwd()
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, export_all = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Reuse gate helpers by sourcing function definitions only up through .run_one
# Avoid executing the driver body: parse and eval selected functions.
src <- readLines(file.path(root, "inst", "benchmarks", "global_gcv",
                           "run_global_lambda_optimize.R"))
# Cut before cfg <- .cfg_for_mode
cut_at <- grep("^cfg <- \\.cfg_for_mode", src)[1]
stopifnot(is.finite(cut_at))
eval(parse(text = src[seq_len(cut_at - 1L)]), envir = environment())

out_dir <- file.path(root, "inst", "benchmarks", "global_gcv",
                     "results", "global_opt_seed31")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- .cfg_for_mode("GATE")
message("=== seed31 stabilize: strong_aniso seed_opt=31 M_final=", cfg$M_final, " ===")
res <- .run_one("strong_aniso", cfg, seed_opt = 31L, ablation = "full",
                verbose = TRUE)

row <- res$row
utils::write.csv(row, file.path(out_dir, "seed31_row.csv"), row.names = FALSE)
utils::write.csv(res$opt$final_candidates,
                 file.path(out_dir, "seed31_final.csv"), row.names = FALSE)

soft <- isTRUE(row$ranking_soft_ok) && isTRUE(row$in_valley_1pct) &&
  isTRUE(row$alt_in_valley_1pct) && isTRUE(row$aniso_ok)
strict <- isTRUE(row$ranking_stable)

verdict <- if (strict) {
  "STRICT_PASS"
} else if (soft) {
  "SOFT_PASS"
} else {
  "FAIL"
}

# If soft fails, retry once with higher M_final
if (identical(verdict, "FAIL")) {
  message("=== retry with M_final=80 ===")
  cfg2 <- cfg
  cfg2$M_final <- 80L
  res2 <- .run_one("strong_aniso", cfg2, seed_opt = 31L, ablation = "full",
                   verbose = TRUE)
  row2 <- res2$row
  utils::write.csv(row2, file.path(out_dir, "seed31_row_M80.csv"),
                   row.names = FALSE)
  soft2 <- isTRUE(row2$ranking_soft_ok) && isTRUE(row2$in_valley_1pct) &&
    isTRUE(row2$alt_in_valley_1pct) && isTRUE(row2$aniso_ok)
  verdict <- if (isTRUE(row2$ranking_stable)) {
    "STRICT_PASS_M80"
  } else if (soft2) {
    "SOFT_PASS_M80"
  } else {
    "FAIL_M80"
  }
  row <- row2
}

cat("\n==== SEED31 VERDICT:", verdict, "====\n")
print(row[, c(
  "seed_opt", "in_valley_1pct", "alt_in_valley_1pct",
  "ranking_stable", "ranking_soft_ok", "aniso_ok",
  "scenario_pass", "grid_on_box_edge", "winner_source",
  "d_theta", "delta_rel"
)])
writeLines(verdict, file.path(out_dir, "VERDICT.txt"))
message("Wrote ", out_dir)
