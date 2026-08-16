# Experimental lab: 2D global TT-gGCV surface vs full-TP GCV and outer cGCV.
# Shared probes across all 49 lambda grid points.
#
# Usage:
#   Rscript inst/benchmarks/global_gcv/run_lambda_surface_d2.R

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
rank <- 5L
degree <- 3L
penalty_order <- 2L
X <- matrix(runif(n * d), n, d)
f_true <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) + 0.25 * X[, 1] * X[, 2]
y <- f_true + rnorm(n, 0, 0.2)

ctrl_fixed <- tt_control(
  max_sweeps = 35L, tol = 1e-10, compute_edf = FALSE, seed = 1L, backend = "R"
)
# Shared probes for the whole surface
probes <- TTPsplines:::.tt_lab_rademacher_probes(n, M = 8L, probe_seed = 11L)
M_use <- 5L
eps_rel <- 1e-3
scheme <- "forward"

log10_grid <- seq(-3, 3, by = 1)
grid <- expand.grid(log10_l1 = log10_grid, log10_l2 = log10_grid, KEEP.OUT.ATTRS = FALSE)

bs <- TTPsplines:::build_marginal_bases(X, k = k, degree = degree)

rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  lam <- 10^c(grid$log10_l1[i], grid$log10_l2[i])
  t0 <- proc.time()[["elapsed"]]
  gcv_tt <- TTPsplines:::tt_global_gcv(
    lambda = lam, y = y, X = X, rank = rank, probes = probes,
    M = M_use, epsilon_rel = eps_rel, scheme = scheme,
    control = ctrl_fixed, backend = "R", k = k, degree = degree,
    penalty_order = penalty_order, on_nonconverged = "na", warm_start = TRUE
  )
  full <- TTPsplines:::.tt_lab_full_tp_gaussian(
    y, bs$basis, lambda = lam, penalty_order = penalty_order
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  yhat <- if (!is.null(gcv_tt$fit)) fitted(gcv_tt$fit) else rep(NA_real_, n)
  rows[[i]] <- data.frame(
    log10_l1 = grid$log10_l1[i],
    log10_l2 = grid$log10_l2[i],
    lambda1 = lam[1],
    lambda2 = lam[2],
    rss_tt = gcv_tt$rss,
    gdf_tt = gcv_tt$gdf,
    gdf_se = gcv_tt$gdf_mc_se,
    global_gcv_tt = gcv_tt$global_gcv,
    valid_tt = gcv_tt$valid,
    invalid = paste(gcv_tt$invalid_reasons, collapse = "|"),
    n_sweeps = gcv_tt$fit$n_sweeps %||% NA_integer_,
    converged = gcv_tt$fit_ok$ok %||% FALSE,
    rss_full = full$rss,
    edf_full = full$edf,
    gcv_full = full$gcv,
    rmse_vs_true = sqrt(mean((yhat - f_true)^2)),
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "(%g,%g) gGCV_TT=%.4g GCV_full=%.4g GDF=%.2f EDF=%.2f valid=%s t=%.1fs",
    grid$log10_l1[i], grid$log10_l2[i],
    gcv_tt$global_gcv, full$gcv, gcv_tt$gdf, full$edf,
    gcv_tt$valid, elapsed
  ))
}
surf <- do.call(rbind, rows)
utils::write.csv(surf, file.path(out_dir, "lambda_surface_d2.csv"), row.names = FALSE)

# Outer simultaneous cGCV reference (public method; for location comparison only)
ctrl_cgcv <- tt_control(
  max_sweeps = 20L, outer_maxit = 8L, compute_edf = FALSE, seed = 1L,
  backend = "R", cgcv_update = "outer_simultaneous",
  lambda_start = 1, warn_lambda_boundary = FALSE
)
fit_cgcv <- ttps(
  y, X, family = gaussian(), rank = rank, k = k, lambda = "cGCV",
  optimizer = "ALS", backend = "R", control = ctrl_cgcv
)
cgcv_row <- data.frame(
  lambda1 = fit_cgcv$lambda[1],
  lambda2 = fit_cgcv$lambda[2],
  log10_l1 = log10(fit_cgcv$lambda[1]),
  log10_l2 = log10(fit_cgcv$lambda[2]),
  lambda_boundary = paste(fit_cgcv$lambda_boundary, collapse = ","),
  stringsAsFactors = FALSE
)
utils::write.csv(cgcv_row, file.path(out_dir, "outer_cgcv_point.csv"), row.names = FALSE)

# Locate minima
ok <- is.finite(surf$global_gcv_tt) & isTRUE(surf$valid_tt %in% TRUE)
if (!any(ok)) ok <- is.finite(surf$global_gcv_tt)
i_tt <- which.min(surf$global_gcv_tt[ok])
i_tt <- which(ok)[i_tt]
i_full <- which.min(surf$gcv_full)
i_true <- which.min(surf$rmse_vs_true)
minima <- data.frame(
  criterion = c("global_TT_gGCV", "full_TP_GCV", "rmse_vs_true", "outer_cGCV"),
  log10_l1 = c(surf$log10_l1[i_tt], surf$log10_l1[i_full],
               surf$log10_l1[i_true], cgcv_row$log10_l1),
  log10_l2 = c(surf$log10_l2[i_tt], surf$log10_l2[i_full],
               surf$log10_l2[i_true], cgcv_row$log10_l2),
  value = c(surf$global_gcv_tt[i_tt], surf$gcv_full[i_full],
            surf$rmse_vs_true[i_true], NA_real_),
  stringsAsFactors = FALSE
)
utils::write.csv(minima, file.path(out_dir, "lambda_surface_minima.csv"), row.names = FALSE)

# Heatmap helper
mat_from <- function(col) {
  m <- matrix(NA_real_, length(log10_grid), length(log10_grid),
              dimnames = list(log10_grid, log10_grid))
  for (r in seq_len(nrow(surf))) {
    m[as.character(surf$log10_l1[r]), as.character(surf$log10_l2[r])] <- surf[[col]][r]
  }
  m
}

png(file.path(fig_dir, "heatmap_global_gcv_tt.png"), width = 560, height = 480)
m <- mat_from("global_gcv_tt")
m[!is.finite(m)] <- NA
image(log10_grid, log10_grid, m, xlab = "log10(lambda1)", ylab = "log10(lambda2)",
      main = "global TT-gGCV (experimental)", col = hcl.colors(20, "Blues"))
points(minima$log10_l1[1], minima$log10_l2[1], pch = 8, col = "red", cex = 1.4)
points(cgcv_row$log10_l1, cgcv_row$log10_l2, pch = 4, col = "darkgreen", cex = 1.4, lwd = 2)
legend("topright", c("min gGCV_TT", "outer cGCV"), pch = c(8, 4),
       col = c("red", "darkgreen"), bty = "n")
dev.off()

png(file.path(fig_dir, "heatmap_compare_tt_vs_full.png"), width = 900, height = 420)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
image(log10_grid, log10_grid, mat_from("global_gcv_tt"),
      xlab = "log10(l1)", ylab = "log10(l2)", main = "TT global-gGCV",
      col = hcl.colors(20, "Blues"))
points(minima$log10_l1[1], minima$log10_l2[1], pch = 8, col = "red")
image(log10_grid, log10_grid, mat_from("gcv_full"),
      xlab = "log10(l1)", ylab = "log10(l2)", main = "full-TP GCV",
      col = hcl.colors(20, "Purples"))
points(minima$log10_l1[2], minima$log10_l2[2], pch = 8, col = "red")
points(cgcv_row$log10_l1, cgcv_row$log10_l2, pch = 4, col = "darkgreen", lwd = 2)
par(op)
dev.off()

message("Minima:")
print(minima)
