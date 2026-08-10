#!/usr/bin/env Rscript
# Ishigami rank × initialization audit
#
# Scientific question: Ishigami is approximately TT-rank ≤ 2 (sum of two
# separable terms). If random-init ALS with few sweeps shows r=3 << r=2 in
# RSS, is that structural need or optimizer trapping?
#
# Protocol:
#   - fixed λ, Gaussian ALS
#   - multi-seed random init for r ∈ {2, 3}
#   - warm-start: best r=3 → tt_truncate_rank(r=2) → ALS
#
# Usage (from package root):
#   Rscript inst/benchmarks/benchmark_ishigami_rank_init.R
#   TTPSPLINES_BENCH_QUICK=1 Rscript inst/benchmarks/benchmark_ishigami_rank_init.R
#
# Outputs: inst/benchmarks/results/ishigami_rank_init_*.csv

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
bench_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  normalizePath("inst/benchmarks")
}
source(file.path(bench_dir, "helpers.R"))
.ttps_bench_ensure_pkg()

quick <- identical(Sys.getenv("TTPSPLINES_BENCH_QUICK", ""), "1")
seeds <- if (quick) 1:5 else 1:20
max_sweeps <- if (quick) 20L else 50L
k <- 8L
lambda <- 1
n_train <- if (quick) 400L else 800L
n_test <- 1000L
sigma <- 0.1

out_dir <- Sys.getenv("TTPSPLINES_BENCH_OUT", unset = file.path(bench_dir, "results"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
# Packaged Ishigami is n=800; redraw for train/test with known truth
train <- simulate_ishigami(n = n_train, sigma = sigma, seed = 1)
test <- simulate_ishigami(n = n_test, sigma = 0, seed = 99) # noiseless test truth
Xtr <- train$X
ytr <- train$y
Xte <- test$X
fte <- test$f

rmse <- function(a, b) sqrt(mean((as.numeric(a) - as.numeric(b))^2))

fit_one <- function(rank, seed, init = NULL, tag = "random") {
  ctrl <- tt_control(
    max_sweeps = max_sweeps,
    backend = "R",
    compute_edf = TRUE,
    seed = seed,
    init_sd = 0.1
  )
  if (is.null(init)) {
    init <- tt_initialize(Xtr, rank = rank, k = k, seed = seed, sd = 0.1)
  }
  t0 <- proc.time()[["elapsed"]]
  fit <- ttps(
    ytr, Xtr,
    family = gaussian(),
    rank = rank,
    k = k,
    lambda = lambda,
    init = init,
    control = ctrl
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  obj <- tryCatch(tt_objective(fit, Xtr, ytr), error = function(e) NULL)
  data.frame(
    tag = tag,
    rank = as.integer(rank),
    seed = as.integer(seed),
    max_sweeps = max_sweeps,
    n_sweeps = fit$n_sweeps %||% NA_integer_,
    converged = isTRUE(fit$converged),
    rss = fit$deviance,
    objective = if (is.null(obj)) NA_real_ else obj$value,
    edf = fit$edf %||% NA_real_,
    npar_tt = fit$npar_tt,
    lambda = paste(sprintf("%.4g", fit$lambda), collapse = ","),
    rmse_train_y = rmse(fitted(fit), ytr),
    rmse_test_f = rmse(predict(fit, Xte, type = "response"), fte),
    elapsed_s = elapsed,
    stringsAsFactors = FALSE
  )
}

message(sprintf(
  "Ishigami rank-init audit | seeds=%d sweeps=%d k=%d n_train=%d quick=%s",
  length(seeds), max_sweeps, k, n_train, quick
))

rows <- list()
for (r in c(2L, 3L)) {
  for (s in seeds) {
    message(sprintf("  random r=%d seed=%d ...", r, s))
    rows[[length(rows) + 1L]] <- fit_one(r, s, tag = "random")
  }
}
tab <- do.call(rbind, rows)

# Best r=3 by test RMSE → truncate → r=2 warm start
tab3 <- tab[tab$rank == 3L & tab$tag == "random", , drop = FALSE]
best_i <- which.min(tab3$rmse_test_f)
best_seed <- tab3$seed[best_i]
message(sprintf("Warm-start from best r=3 seed=%d (test RMSE=%.4g)",
                best_seed, tab3$rmse_test_f[best_i]))

ctrl3 <- tt_control(
  max_sweeps = max_sweeps, backend = "R", compute_edf = FALSE,
  seed = best_seed, init_sd = 0.1
)
fit3_best <- ttps(
  ytr, Xtr, family = gaussian(), rank = 3, k = k, lambda = lambda,
  init = tt_initialize(Xtr, rank = 3, k = k, seed = best_seed, sd = 0.1),
  control = ctrl3
)
init2_ws <- tt_truncate_rank(fit3_best$cores, rank = 2)
rows[[length(rows) + 1L]] <- fit_one(2L, seed = best_seed, init = init2_ws, tag = "warm_trunc_r3")
# Also refit same truncated init with a few more seeds of ALS noise? keep single warm path.

tab <- do.call(rbind, rows)
summary_path <- file.path(out_dir, "ishigami_rank_init_runs.csv")
write.csv(tab, summary_path, row.names = FALSE)

agg <- aggregate(
  cbind(rss, rmse_test_f, edf, objective, converged) ~ tag + rank,
  data = tab,
  FUN = function(z) c(mean = mean(z, na.rm = TRUE),
                      median = stats::median(z, na.rm = TRUE),
                      min = min(z, na.rm = TRUE),
                      max = max(z, na.rm = TRUE))
)
# flatten
agg_flat <- do.call(rbind, lapply(seq_len(nrow(agg)), function(i) {
  data.frame(
    tag = agg$tag[i],
    rank = agg$rank[i],
    rss_mean = agg$rss[i, "mean"],
    rss_median = agg$rss[i, "median"],
    rss_min = agg$rss[i, "min"],
    rmse_test_mean = agg$rmse_test_f[i, "mean"],
    rmse_test_min = agg$rmse_test_f[i, "min"],
    edf_mean = agg$edf[i, "mean"],
    obj_mean = agg$objective[i, "mean"],
    conv_rate = agg$converged[i, "mean"],
    stringsAsFactors = FALSE
  )
}))
agg_path <- file.path(out_dir, "ishigami_rank_init_summary.csv")
write.csv(agg_flat, agg_path, row.names = FALSE)

message("Wrote:\n  ", summary_path, "\n  ", agg_path)
print(agg_flat, row.names = FALSE)

# Interpretation hints
r2 <- tab[tab$tag == "random" & tab$rank == 2L, ]
r3 <- tab[tab$tag == "random" & tab$rank == 3L, ]
ws <- tab[tab$tag == "warm_trunc_r3", ]
message("\n--- Audit read ---")
message(sprintf("random r=2: median RSS=%.4g  min RSS=%.4g  min test RMSE=%.4g",
                median(r2$rss), min(r2$rss), min(r2$rmse_test_f)))
message(sprintf("random r=3: median RSS=%.4g  min RSS=%.4g  min test RMSE=%.4g",
                median(r3$rss), min(r3$rss), min(r3$rmse_test_f)))
if (nrow(ws)) {
  message(sprintf("warm trunc r=2: RSS=%.4g  test RMSE=%.4g",
                  ws$rss[1], ws$rmse_test_f[1]))
}
succ <- function(df, thr) mean(df$rss < thr)
message(sprintf(
  "Success rates P(RSS<100): r2=%.0f%%  r3=%.0f%% | P(RSS<150): r2=%.0f%%  r3=%.0f%%",
  100 * succ(r2, 100), 100 * succ(r3, 100),
  100 * succ(r2, 150), 100 * succ(r3, 150)
))
message(
  "If warm/min r=2 RSS approaches r=3: structural capacity OK (init/sweeps).\n",
  "If r=2 stays ~10x worse: investigate parameterization/penalty before paper claim.\n",
  "Paper framing: minimal representational rank vs computationally robust working rank."
)

# Success-rate table for paper figure
succ_tab <- data.frame(
  rank = c(2L, 3L),
  n_starts = c(nrow(r2), nrow(r3)),
  p_rss_lt_100 = c(succ(r2, 100), succ(r3, 100)),
  p_rss_lt_150 = c(succ(r2, 150), succ(r3, 150)),
  rss_median = c(median(r2$rss), median(r3$rss)),
  rss_min = c(min(r2$rss), min(r3$rss)),
  stringsAsFactors = FALSE
)
succ_path <- file.path(out_dir, "ishigami_rank_init_success.csv")
write.csv(succ_tab, succ_path, row.names = FALSE)
message("Also wrote: ", succ_path)
