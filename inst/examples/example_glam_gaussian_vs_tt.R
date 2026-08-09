# Small GLAM-vs-TT validation (d=3, d=5) with honest timings
#   Rscript inst/examples/example_glam_gaussian_vs_tt.R

root <- if (file.exists("DESCRIPTION") &&
             identical(unname(read.dcf("DESCRIPTION")[, "Package"]), "TTPsplines")) {
  normalizePath(".")
} else if (file.exists("../../DESCRIPTION")) {
  normalizePath("../..")
} else NULL

if (!is.null(root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, quiet = TRUE)
} else library(TTPsplines)

show <- function(cmp, label) {
  cat("\n===", label, "===\n")
  cols <- c("method", "rank", "npar", "compression", "rmse_truth",
            "fit_time", "rank_selection_time", "total_procedure_time",
            "oracle_probe_time", "delta_rmse_oracle", "status")
  print(as.data.frame(cmp)[, cols], digits = 4, row.names = FALSE)
  cat("\nSummary:\n")
  print(summarize_glam_tt_compare(cmp), digits = 4, row.names = FALSE)
}

cmp3 <- compare_glam_tt_gaussian(
  d = 3, n_grid = c(10, 10, 10), k = 6, ranks = 1:4,
  lambda = 1, folds = 4, max_sweeps = 5, seed = 1, n_test = 800
)
show(cmp3, "d=3")

cmp5 <- compare_glam_tt_gaussian(
  d = 5, n_grid = 6, k = 5, ranks = 1:4,
  lambda = 1, folds = 4, max_sweeps = 5, seed = 1, n_test = 1000
)
show(cmp5, "d=5")

cat("\nTiming interpretation (d=5, TT-1SE):\n")
se1 <- cmp5[cmp5$method == "TT-1SE", ][1, ]
cat(sprintf("  fit_time              = %.3fs  (final model only)\n", se1$fit_time))
cat(sprintf("  rank_selection_time   = %.3fs  (K-fold CV over ranks)\n", se1$rank_selection_time))
cat(sprintf("  total_procedure_time  = %.3fs  (what a user actually pays)\n", se1$total_procedure_time))
ora <- cmp5[cmp5$method == "TT-oracle", ][1, ]
cat(sprintf("  oracle_probe_time     = %.3fs  (simulation only; NOT in totals)\n",
            ora$oracle_probe_time))
