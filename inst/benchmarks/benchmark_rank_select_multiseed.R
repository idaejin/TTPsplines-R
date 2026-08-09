# Multi-seed audit: oracle vs min-CV vs 1-SE
#   Rscript inst/benchmarks/benchmark_rank_select_multiseed.R
#
# Default: 30 seeds × {Ishigami, Sobol-g, Friedman} + one GLAM-vs-TT cell.
# Override: N_SEEDS=50 Rscript ...

root <- if (file.exists("DESCRIPTION") &&
             identical(unname(read.dcf("DESCRIPTION")[, "Package"]), "TTPsplines")) {
  normalizePath(".")
} else NULL
if (!is.null(root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, quiet = TRUE)
} else library(TTPsplines)

n_seeds <- as.integer(Sys.getenv("N_SEEDS", "30"))
out_dir <- file.path(root %||% ".", "inst/benchmarks/output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

audit_surface_seed <- function(name, sim_fun, ranks = 1:5, k = 6,
                               folds = 5, max_sweeps = 8, seed) {
  sim <- sim_fun(n = 400L, seed = seed)
  X <- as.matrix(sim$X)
  y <- as.numeric(sim$y)
  f <- as.numeric(sim$f)
  ctrl <- tt_control(max_sweeps = max_sweeps, compute_edf = FALSE)
  t0 <- proc.time()[["elapsed"]]
  sel <- tt_rank_select(y, X, ranks = ranks, k = k, lambda = 1, folds = folds,
                        rule = "1se", seed = seed + 17L, control = ctrl)
  cv_time <- proc.time()[["elapsed"]] - t0
  oracle_rmse <- vapply(ranks, function(r) {
    fit <- ttpspline(y, X, rank = r, k = k, lambda = 1, control = ctrl)
    sqrt(mean((fitted(fit) - f)^2))
  }, numeric(1))
  r_ora <- ranks[which.min(oracle_rmse)]
  fit_ora <- ttpspline(y, X, rank = r_ora, k = k, lambda = 1, control = ctrl)
  fit_min <- tt_rank_refit(sel, rank = sel$rank_min)
  fit_se1 <- tt_rank_refit(sel)
  rmse <- function(fit) sqrt(mean((fitted(fit) - f)^2))
  data.frame(
    surface = name, seed = seed,
    oracle_rank = r_ora, minCV_rank = sel$rank_min, se1_rank = sel$rank_1se,
    oracle_rmse = rmse(fit_ora), minCV_rmse = rmse(fit_min), se1_rmse = rmse(fit_se1),
    se1_npar = fit_se1$npar_tt, se1_CR = fit_se1$compression_ratio,
    cv_time_s = cv_time, stringsAsFactors = FALSE
  )
}

surfaces <- list(
  list(name = "Ishigami", fun = function(n, seed) simulate_ishigami(n = n, seed = seed), k = 8),
  list(name = "Sobol-g", fun = function(n, seed) simulate_sobol_g(n = n, seed = seed), k = 6),
  list(name = "Friedman", fun = function(n, seed) simulate_friedman(n = n, seed = seed), k = 5)
)

cat(sprintf("Multi-seed surface audit: %d seeds\n", n_seeds))
rows <- list()
for (s in surfaces) {
  for (seed in seq_len(n_seeds)) {
    rows[[length(rows) + 1L]] <- audit_surface_seed(
      s$name, s$fun, k = s$k, seed = seed
    )
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(out_dir, "rank_select_multiseed_raw.csv"), row.names = FALSE)

summarize_surface <- function(df) {
  data.frame(
    surface = df$surface[1],
    n_seeds = nrow(df),
    # selection frequencies
    freq_se1_r1 = mean(df$se1_rank == 1),
    freq_se1_r2 = mean(df$se1_rank == 2),
    freq_se1_r3 = mean(df$se1_rank == 3),
    freq_se1_r4 = mean(df$se1_rank == 4),
    freq_se1_r5 = mean(df$se1_rank == 5),
    mean_oracle_rank = mean(df$oracle_rank),
    mean_minCV_rank = mean(df$minCV_rank),
    mean_se1_rank = mean(df$se1_rank),
    mean_oracle_rmse = mean(df$oracle_rmse),
    mean_minCV_rmse = mean(df$minCV_rmse),
    mean_se1_rmse = mean(df$se1_rmse),
    mean_se1_penalty = mean(df$se1_rmse - df$oracle_rmse),
    mean_se1_npar = mean(df$se1_npar),
    mean_se1_CR = mean(df$se1_CR),
    mean_cv_time = mean(df$cv_time_s),
    stringsAsFactors = FALSE
  )
}
summ <- do.call(rbind, lapply(split(tab, tab$surface), summarize_surface))
rownames(summ) <- NULL
write.csv(summ, file.path(out_dir, "rank_select_multiseed_summary.csv"), row.names = FALSE)
cat("\nSurface summary:\n")
print(summ, digits = 3, row.names = FALSE)

# One GLAM-vs-TT cell across seeds
cat(sprintf("\nGLAM-vs-TT multi-seed (d=3, k=8, %d seeds):\n", n_seeds))
grows <- lapply(seq_len(n_seeds), function(seed) {
  tabg <- compare_glam_tt_gaussian(
    d = 3, n_grid = 12, k = 8, ranks = 1:3, max_sweeps = 5,
    select_folds = 3, seed = 1000L + seed
  )
  s <- summarize_glam_tt_compare(tabg)
  s$seed <- seed
  s
})
gtab <- do.call(rbind, grows)
write.csv(gtab, file.path(out_dir, "glam_tt_multiseed_d3k8.csv"), row.names = FALSE)
cat(sprintf(
  "  mean RMSE  glam=%.4f  oracle=%.4f  minCV=%.4f  1SE=%.4f\n  mean ranks oracle=%.2f minCV=%.2f 1SE=%.2f\n  mean 1SE penalty vs oracle=%.4f | mean select_time=%.2fs\n",
  mean(gtab$glam_rmse), mean(gtab$oracle_rmse), mean(gtab$minCV_rmse), mean(gtab$se1_rmse),
  mean(gtab$oracle_rank), mean(gtab$minCV_rank), mean(gtab$se1_rank),
  mean(gtab$se1_rmse - gtab$oracle_rmse), mean(gtab$select_time)
))
