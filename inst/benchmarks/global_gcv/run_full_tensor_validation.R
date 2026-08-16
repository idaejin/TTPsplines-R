# Experimental lab: validate MC global GDF vs exact full-TP EDF (Gaussian, d=2).
# NOT a public method. Source from package root via helpers.
#
# Usage (from package root):
#   Rscript inst/benchmarks/global_gcv/run_full_tensor_validation.R

source(file.path("inst", "benchmarks", "helpers.R"))
.ttps_bench_ensure_pkg()

out_dir <- file.path("inst", "benchmarks", "global_gcv", "results")
fig_dir <- file.path("inst", "benchmarks", "global_gcv", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260816)
n <- 80L
d <- 2L
k <- 5L
degree <- 3L
penalty_order <- 2L
lambda <- c(0.5, 2)
X <- matrix(runif(n * d), n, d)
colnames(X) <- paste0("x", seq_len(d))
f_true <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.3 * X[, 1] * X[, 2]
y <- f_true + rnorm(n, 0, 0.2)

ctrl <- tt_control(
  max_sweeps = 40L, tol = 1e-10, compute_edf = FALSE, seed = 1L,
  backend = "R"
)
probes <- TTPsplines:::.tt_lab_rademacher_probes(n, M = 20L, probe_seed = 11L)

# Full TP reference (same knots / bases as ttps would build)
bs <- TTPsplines:::build_marginal_bases(
  X, k = k, degree = degree, knots = NULL, cyclic = NULL, period = NULL
)
full <- TTPsplines:::.tt_lab_full_tp_gaussian(
  y, bs$basis, lambda = lambda, penalty_order = penalty_order
)

run_case <- function(rank, label) {
  t0 <- proc.time()[["elapsed"]]
  fit <- ttps(
    y, X, family = gaussian(), rank = rank, k = k, degree = degree,
    penalty_order = penalty_order, lambda = lambda, optimizer = "ALS",
    backend = "R", control = ctrl
  )
  g_fwd <- TTPsplines:::tt_global_gdf_mc(
    fit, y = y, M = 10L, probes = probes[, 1:10, drop = FALSE],
    scheme = "forward", epsilon_rel = 1e-3, warm_start = TRUE,
    on_nonconverged = "na", control = ctrl
  )
  g_cen <- TTPsplines:::tt_global_gdf_mc(
    fit, y = y, M = 10L, probes = probes[, 1:10, drop = FALSE],
    scheme = "central", epsilon_rel = 1e-3, warm_start = TRUE,
    on_nonconverged = "na", control = ctrl
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  yhat <- fitted(fit)
  data.frame(
    label = label,
    rank = rank,
    n = n, k = k, d = d,
    lambda = paste(lambda, collapse = ","),
    edf_full = full$edf,
    rss_full = full$rss,
    gcv_full = full$gcv,
    rss_tt = sum((y - yhat)^2),
    rmse_fit_vs_full = sqrt(mean((yhat - full$fitted)^2)),
    cor_fit_vs_full = suppressWarnings(cor(yhat, full$fitted)),
    gdf_forward = g_fwd$gdf,
    gdf_forward_se = g_fwd$gdf_mc_se,
    gdf_central = g_cen$gdf,
    gdf_central_se = g_cen$gdf_mc_se,
    M_ok_fwd = g_fwd$M_ok,
    M_ok_cen = g_cen$M_ok,
    n_refits_fwd = g_fwd$n_refits,
    n_refits_cen = g_cen$n_refits,
    used_conditional_edf = g_fwd$used_conditional_edf,
    converged = fit$converged,
    n_sweeps = fit$n_sweeps,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

# Rank sufficient for d=2 matrices: r = k; rank-restricted: r = 1
tab <- rbind(
  run_case(rank = k, label = "rank_sufficient"),
  run_case(rank = 1L, label = "rank_restricted")
)

csv_path <- file.path(out_dir, "full_tensor_validation.csv")
utils::write.csv(tab, csv_path, row.names = FALSE)

# Simple comparison figure
png(file.path(fig_dir, "gdf_vs_edf_full.png"), width = 720, height = 420)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (i in seq_len(nrow(tab))) {
  vals <- c(tab$edf_full[i], tab$gdf_forward[i], tab$gdf_central[i])
  names(vals) <- c("EDF_full", "GDF_fwd", "GDF_cen")
  barplot(vals, main = tab$label[i], ylab = "degrees of freedom", col = "gray80")
  abline(h = tab$edf_full[i], lty = 2, col = "blue")
}
par(op)
dev.off()

message("Wrote ", csv_path)
print(tab)
