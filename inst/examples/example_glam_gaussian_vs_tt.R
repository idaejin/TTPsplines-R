# Gaussian GLAM vs TTPsplines on grids (d = 3 and d = 5)
#   Rscript inst/examples/example_glam_gaussian_vs_tt.R

root <- if (file.exists("DESCRIPTION") &&
             identical(unname(read.dcf("DESCRIPTION")[, "Package"]), "TTPsplines")) {
  normalizePath(".")
} else if (file.exists("../../DESCRIPTION")) {
  normalizePath("../..")
} else {
  NULL
}

if (!is.null(root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, quiet = TRUE)
} else {
  library(TTPsplines)
}

cols <- c("method", "npar", "npar_dense", "compression",
          "rmse_truth", "rmse_y", "rmse_vs_glam", "time_s", "time_vs_glam",
          "optimizer")

cat("=== d = 3 (Currie–Durbán–Eilers GLAM vs TT) ===\n")
tab3 <- compare_glam_tt_gaussian(d = 3, ranks = 1:4, max_sweeps = 10)
print(tab3[, cols], digits = 4, row.names = FALSE)

cat("\n=== d = 5 ===\n")
tab5 <- compare_glam_tt_gaussian(d = 5, ranks = 1:3, max_sweeps = 8)
print(tab5[, cols], digits = 4, row.names = FALSE)

cat("\nSummary (best TT by rmse_truth):\n")
for (tab in list(tab3, tab5)) {
  tt <- tab[!is.na(tab$rank), ]
  best <- tt[which.min(tt$rmse_truth), ]
  glam <- tab[tab$method == "GLAM", ]
  cat(sprintf(
    "  d=%d | GLAM npar=%d time=%.3fs RMSE=%.4f | best %s npar=%d CR=%.1fx time=%.3fs RMSE=%.4f\n",
    glam$d, glam$npar, glam$time_s, glam$rmse_truth,
    best$method, best$npar, best$compression, best$time_s, best$rmse_truth
  ))
}
