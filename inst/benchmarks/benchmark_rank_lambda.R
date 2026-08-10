# Rank × λ interaction (fixed λ grid, then cGCV)
#   Rscript inst/benchmarks/benchmark_rank_lambda.R
#
# Separates structural capacity (r) from smoothness (λ).

root <- if (file.exists("DESCRIPTION") &&
             identical(unname(read.dcf("DESCRIPTION")[, "Package"]), "TTPsplines")) {
  normalizePath(".")
} else NULL
if (!is.null(root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, quiet = TRUE)
} else library(TTPsplines)

out_dir <- file.path(if (!is.null(root)) root else ".", "inst/benchmarks/output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data(friedman)
idx <- seq_len(400)
X <- as.matrix(friedman[idx, paste0("x", 1:5)])
y <- friedman$y[idx]
f <- friedman$f[idx]
ranks <- 1:4
k <- 5L
ctrl <- tt_control(max_sweeps = 8L, compute_edf = FALSE,
                   warn_lambda_boundary = FALSE)

row_fit <- function(phase, lambda_lab, rank, fit, time_s) {
  data.frame(
    phase = phase,
    lambda_lab = lambda_lab,
    rank = as.integer(rank),
    npar = as.integer(fit$npar_tt),
    compression = fit$compression_ratio,
    rmse_truth = sqrt(mean((fitted(fit) - f)^2)),
    rmse_y = sqrt(mean((fitted(fit) - y)^2)),
    time_s = time_s,
    stringsAsFactors = FALSE
  )
}

# --- Phase A: fixed λ grid -------------------------------------------------
lambdas <- c(0.1, 1, 10)
rows <- list()
for (lam in lambdas) {
  for (r in ranks) {
    t0 <- proc.time()[["elapsed"]]
    fit <- ttps(y, X, rank = r, k = k, lambda = lam, control = ctrl)
    rows[[length(rows) + 1L]] <- row_fit("fixed", as.character(lam), r, fit,
                                         proc.time()[["elapsed"]] - t0)
  }
  t1 <- proc.time()[["elapsed"]]
  sel <- tt_rank_select(y, X, ranks = ranks, k = k, lambda = lam, folds = 4,
                        rule = "1se", seed = 1, control = ctrl)
  fit_min <- tt_rank_refit(sel, rank = sel$rank_min)
  fit_se1 <- tt_rank_refit(sel)
  elapsed <- proc.time()[["elapsed"]] - t1
  rows[[length(rows) + 1L]] <- row_fit("fixed-minCV", as.character(lam),
                                       sel$rank_min, fit_min, elapsed)
  rows[[length(rows) + 1L]] <- row_fit("fixed-1SE", as.character(lam),
                                       sel$rank_1se, fit_se1, NA_real_)
}
tab_fixed <- do.call(rbind, rows)

# --- Phase B: cGCV ---------------------------------------------------------
rows_c <- list()
for (r in ranks) {
  t0 <- proc.time()[["elapsed"]]
  fit <- ttps(y, X, rank = r, k = k, lambda = "cGCV", control = ctrl)
  rows_c[[length(rows_c) + 1L]] <- row_fit(
    "cGCV", paste(sprintf("%.3g", fit$lambda), collapse = ","), r, fit,
    proc.time()[["elapsed"]] - t0
  )
}
t1 <- proc.time()[["elapsed"]]
sel_c <- tt_rank_select(y, X, ranks = ranks, k = k, lambda = "cGCV", folds = 4,
                        rule = "1se", seed = 1, control = ctrl)
fit_minc <- tt_rank_refit(sel_c, rank = sel_c$rank_min)
fit_se1c <- tt_rank_refit(sel_c)
elapsed <- proc.time()[["elapsed"]] - t1
rows_c[[length(rows_c) + 1L]] <- row_fit("cGCV-minCV", "cGCV-per-fold",
                                         sel_c$rank_min, fit_minc, elapsed)
rows_c[[length(rows_c) + 1L]] <- row_fit("cGCV-1SE", "cGCV-per-fold",
                                         sel_c$rank_1se, fit_se1c, NA_real_)
tab_cgcv <- do.call(rbind, rows_c)

tab <- rbind(tab_fixed, tab_cgcv)
write.csv(tab, file.path(out_dir, "rank_lambda_friedman.csv"), row.names = FALSE)

cat("=== Fixed λ × rank grid (Friedman) ===\n")
print(tab_fixed[tab_fixed$phase == "fixed", ], digits = 4, row.names = FALSE)
cat("\nCV-selected ranks at each fixed λ:\n")
print(tab_fixed[tab_fixed$phase %in% c("fixed-minCV", "fixed-1SE"),
                c("phase", "lambda_lab", "rank", "npar", "rmse_truth")],
      digits = 4, row.names = FALSE)
cat("\n=== cGCV × rank ===\n")
print(tab_cgcv, digits = 4, row.names = FALSE)
