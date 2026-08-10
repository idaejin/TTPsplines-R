# Chicago TT Poisson — diagnostic and slice plots
#   Rscript inst/examples/example_chicago_plots.R
#
# Notes on scales:
# - We do NOT center covariates in this script. Axes use the columns as stored
#   in gamair::chicago. Pollution (*median) and time already take values that
#   straddle zero in that dataset (residualized / day-index style coding).
# - Partial-slice x-grids are trimmed to the central 2.5%–97.5% of each
#   covariate so sparse extremes (esp. so2) do not dominate Poisson response
#   bands via exp(η ± z SE).
#
# Outputs: inst/examples/chicago_plots/ (+ lab manuscript/figures/ copy)

root <- (function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), "../..")))
  }
  normalizePath(".")
})()

stopifnot(requireNamespace("devtools", quietly = TRUE))
stopifnot(requireNamespace("gamair", quietly = TRUE))
devtools::load_all(root, quiet = TRUE)
library(gamair)
data(chicago)

d <- stats::na.omit(chicago[, c(
  "death", "tmpd", "o3median", "pm10median", "so2median", "time"
)])
y <- d$death
# Keep gamair column names / scales exactly (no centering / scaling)
X <- as.matrix(d[, c("tmpd", "o3median", "pm10median", "so2median", "time")])
xlab_j <- c(
  "tmpd (°F)",
  "o3median (gamair scale)",
  "pm10median (gamair scale)",
  "so2median (gamair scale)",
  "time (days, gamair)"
)

fit <- ttps(
  y, X,
  family = poisson(),
  rank = 3L,
  k = 5L,
  lambda = "cGCV",
  control = tt_control(
    pirls_maxit = 25L,
    als_sweeps_per_pirls = 3L,
    compute_edf = TRUE,
    seed = 1L,
    trace = FALSE
  )
)
cat("cGCV lambda:", paste(sprintf("%.4g", fit$lambda), collapse = ", "), "\n")
cat(sprintf("EDF=%.2f | deviance=%.1f | boundary=%s\n",
            fit$edf, fit$deviance, paste(fit$lambda_boundary, collapse = ",")))

out_dir <- file.path(root, "inst", "examples", "chicago_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
lab_fig <- file.path(dirname(root), "TTPsplines", "manuscript", "figures")
dir.create(lab_fig, recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, file.path(out_dir, "fit_chicago_cgcv.rds"))

mu <- fitted(fit)
pear <- residuals(fit, type = "pearson")
med <- apply(X, 2, median)

# Data-supported grid on the *observed* covariate scale (no recentering)
.grid_j <- function(j, n = 120L, probs = c(0.025, 0.975)) {
  xj <- X[, j]
  lohi <- stats::quantile(xj, probs = probs, names = FALSE, type = 7)
  seq(lohi[1], lohi[2], length.out = n)
}

heatmap2 <- function(j1, j2, n = 50L) {
  g1 <- .grid_j(j1, n)
  g2 <- .grid_j(j2, n)
  grid <- as.matrix(expand.grid(g1, g2))
  Xm <- matrix(rep(med, each = nrow(grid)), nrow(grid), 5L, byrow = FALSE)
  colnames(Xm) <- colnames(X)
  Xm[, j1] <- grid[, 1]
  Xm[, j2] <- grid[, 2]
  z <- matrix(predict(fit, Xm, type = "link"), n, n)
  list(g1 = g1, g2 = g2, z = z)
}

# ---- 1. time series with pointwise CI (x = gamair time as stored) ----
ord <- order(X[, "time"])
pr_all <- predict(
  fit, X, type = "response",
  se.fit = TRUE, interval = "confidence", level = 0.95
)

png(file.path(out_dir, "chicago_time_fit.png"), width = 1100, height = 450, res = 120)
op <- par(mar = c(4, 4, 2.5, 1))
plot(X[ord, "time"], y[ord], type = "l", col = adjustcolor("grey50", 0.7), lwd = 0.6,
     xlab = "time (days, gamair::chicago as stored)", ylab = "daily deaths",
     main = sprintf(
       "Chicago: TT Poisson + cGCV (rank 3, k = 5; EDF = %.0f)",
       fit$edf
     ))
polygon(
  c(X[ord, "time"], rev(X[ord, "time"])),
  c(pr_all$lower[ord], rev(pr_all$upper[ord])),
  col = adjustcolor("#C0392B", 0.22), border = NA
)
lines(X[ord, "time"], pr_all$fit[ord], col = "#C0392B", lwd = 1.4)
legend("topright", c("observed", "TT fitted", "pointwise 95% CI"),
       col = c("grey50", "#C0392B", adjustcolor("#C0392B", 0.4)),
       lwd = c(1, 1.5, 8), bty = "n")
par(op)
dev.off()

# ---- 2. 1D partial slices with bands (trimmed to central 95% of each x) ----
slice1d_ci <- function(j, n = 120L) {
  g <- .grid_j(j, n)
  Xm <- matrix(rep(med, each = length(g)), length(g), ncol(X), byrow = FALSE)
  colnames(Xm) <- colnames(X)
  Xm[, j] <- g
  pr <- predict(fit, Xm, type = "response", interval = "confidence", level = 0.95)
  list(x = g, fit = pr$fit, lower = pr$lower, upper = pr$upper)
}

slices <- lapply(seq_len(5L), slice1d_ci)
# Common y-scale from data-supported fits (avoid one sparse-tail CI warping panels)
y_cap <- max(stats::quantile(y, 0.999), max(vapply(slices, function(s) max(s$fit), 1)) * 1.25)
ylim_common <- c(
  max(0, min(vapply(slices, function(s) min(s$lower), 1))),
  y_cap
)

png(file.path(out_dir, "chicago_partial_slices.png"), width = 1100, height = 700, res = 120)
op <- par(mfrow = c(2, 3), mar = c(4.2, 4, 2.5, 1))
for (j in 1:5) {
  sl <- slices[[j]]
  # Per-panel: clip displayed band to common death scale; keep curve
  lo <- pmax(sl$lower, 0)
  up <- pmin(sl$upper, y_cap)
  plot(sl$x, sl$fit, type = "n", ylim = ylim_common, xlim = range(sl$x),
       xlab = xlab_j[j], ylab = expression(hat(mu)),
       main = sprintf("slice vs %s (others at median)", colnames(X)[j]))
  polygon(c(sl$x, rev(sl$x)), c(lo, rev(up)),
          col = adjustcolor("#1F4E79", 0.25), border = NA)
  lines(sl$x, sl$fit, lwd = 2, col = "#1F4E79")
  rug(X[, j], col = adjustcolor("grey40", 0.25), quiet = TRUE)
  # Mark that x-grid is central 95% of observed values
  abline(v = range(sl$x), col = adjustcolor("grey50", 0.4), lty = 3)
}
plot.new()
legend(
  "center",
  c(
    sprintf("n = %d | deviance = %.0f | TT params = %d (full %d)",
            fit$n, fit$deviance, fit$npar_tt, fit$npar_dense),
    "x: gamair scales as stored (not re-centered here)",
    "slice grid: central 95% of each covariate",
    "bands: pointwise 95% CI (link SE → inverse link); y capped"
  ),
  bty = "n", cex = 0.9
)
par(op)
dev.off()

# ---- 3. 2D heatmaps on log-mean ----
png(file.path(out_dir, "chicago_heatmaps.png"), width = 1100, height = 900, res = 120)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
pairs <- list(c(1L, 2L), c(1L, 5L), c(2L, 3L), c(3L, 4L))
titles <- c("tmpd x o3 (log mean)", "tmpd x time", "o3 x pm10", "pm10 x so2")
for (i in seq_along(pairs)) {
  h <- heatmap2(pairs[[i]][1], pairs[[i]][2])
  image(h$g1, h$g2, h$z, col = hcl.colors(40, "Blues"),
        xlab = xlab_j[pairs[[i]][1]], ylab = xlab_j[pairs[[i]][2]],
        main = titles[i])
  contour(h$g1, h$g2, h$z, add = TRUE, col = adjustcolor("black", 0.45), labcex = 0.7)
}
par(op)
dev.off()

# ---- 4. diagnostics ----
png(file.path(out_dir, "chicago_diagnostics.png"), width = 1100, height = 700, res = 120)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
plot(mu, y, pch = 16, cex = 0.35, col = adjustcolor("black", 0.35),
     xlab = "fitted mean", ylab = "observed deaths",
     main = "Observed vs fitted")
abline(0, 1, col = "#C0392B", lwd = 2)
smoothScatter(mu, pear, xlab = "fitted mean", ylab = "Pearson residual",
              main = "Pearson residuals vs fitted")
abline(h = 0, col = "#C0392B")
qqnorm(pear, pch = 16, cex = 0.4, main = "Pearson residual QQ")
qqline(pear, col = "#C0392B")
hist(pear, breaks = 40, col = "grey80", border = "white",
     main = "Pearson residual histogram", xlab = "Pearson residual")
par(op)
dev.off()

pngs <- list.files(out_dir, pattern = "[.]png$", full.names = TRUE)
if (dir.exists(lab_fig)) {
  for (f in pngs) file.copy(f, file.path(lab_fig, basename(f)), overwrite = TRUE)
}

cat("Wrote:\n")
cat(paste0("  ", pngs), sep = "\n")
if (dir.exists(lab_fig)) cat("Copied to ", lab_fig, "\n", sep = "")
cat("\nScale note: covariates plotted as in gamair::chicago (no extra centering).\n")
cat("so2 99% quantile ≈", round(stats::quantile(X[, "so2median"], 0.99), 2),
    "; max =", round(max(X[, "so2median"]), 2),
    "(old full-range slice hit the sparse tail).\n")
