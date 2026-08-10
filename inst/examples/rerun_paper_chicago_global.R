# Re-run Paper-1 preliminary illustration + Chicago analysis
# under the classical global discrete P-spline penalty (P_k^full).
#
#   cd 01_PROJECTS/ttpsplines-pkg
#   Rscript inst/examples/rerun_paper_chicago_global.R
#
# Outputs:
#   - lab manuscript/figures/paper1_*.png, chicago_*.png
#   - lab outputs/paper1_*.csv
#   - inst/examples/chicago_plots/*.rds + copies of figures

root <- (function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), "../..")))
  }
  normalizePath(".")
})()

stopifnot(requireNamespace("devtools", quietly = TRUE))
devtools::load_all(root, quiet = TRUE)

lab <- normalizePath(file.path(dirname(root), "TTPsplines"), mustWork = TRUE)
fig_dir <- file.path(lab, "manuscript", "figures")
out_dir <- file.path(lab, "outputs")
chi_dir <- file.path(root, "inst", "examples", "chicago_plots")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chi_dir, recursive = TRUE, showWarnings = FALSE)

ctrl_base <- function(max_sweeps = 25L, ...) {
  tt_control(
    compute_edf = TRUE,
    edf_max_npar = 5000L,
    seed = 1L,
    trace = FALSE,
    max_sweeps = max_sweeps,
    ...
  )
}

cat("=== Penalty mode (default) ===\n")
cat("tt_control()$penalty_mode =", tt_control()$penalty_mode, "\n\n")

########################################################################
# A. Paper-1 preliminary fixed-lambda illustration (global J_λ)
########################################################################
cat("=== A. Paper-1 fixed-λ illustration (global penalty) ===\n")

true_surface <- function(x1, x2, x3) {
  1.5 * sin(2 * pi * x1) * (x2 - 0.5) +
    2 * (x2 - 0.5) * (x3 - 0.5) +
    8 * sin(2 * pi * x1) * (x2 - 0.5) * (x3 - 0.5)
}

make_example_data <- function(n = 1000L, sigma = 0.3, seed = 4L) {
  set.seed(seed)
  X <- cbind(runif(n), runif(n), runif(n))
  truth <- true_surface(X[, 1], X[, 2], X[, 3])
  list(X = X, truth = truth, y = truth + rnorm(n, 0, sigma))
}

# Dense scattered TP P-spline (same discrete difference penalties)
fit_dense_tp3 <- function(y, X, k = 8L, lambda = 1, degree = 3L, penalty_order = 2L) {
  d <- ncol(X)
  bases <- vector("list", d)
  pens <- vector("list", d)
  for (j in seq_len(d)) {
    knots <- TTPsplines:::make_knots(X[, j], k = k, degree = degree)
    bases[[j]] <- TTPsplines:::bspline_basis(X[, j], knots = knots, degree = degree)
    D <- diff(diag(k), differences = penalty_order)
    pens[[j]] <- crossprod(D)
  }
  n <- length(y)
  p1 <- ncol(bases[[1]]); p2 <- ncol(bases[[2]]); p3 <- ncol(bases[[3]])
  B <- matrix(0, n, p1 * p2 * p3)
  for (i in seq_len(n)) {
    B[i, ] <- as.numeric(kronecker(bases[[3]][i, ],
                          kronecker(bases[[2]][i, ], bases[[1]][i, ])))
  }
  S <- lambda * (
    kronecker(diag(p3), kronecker(diag(p2), pens[[1]])) +
      kronecker(diag(p3), kronecker(pens[[2]], diag(p1))) +
      kronecker(pens[[3]], kronecker(diag(p2), diag(p1)))
  )
  XtX <- crossprod(B)
  ridge <- 1e-8 * mean(diag(XtX))
  coef <- solve(XtX + S + diag(ridge, nrow(XtX)), crossprod(B, y))
  mu <- as.numeric(B %*% coef)
  list(mu = mu, npar = length(coef), elapsed = NA_real_)
}

dat <- make_example_data()
k_paper <- 8L
lambda_paper <- 1
ranks <- 1:5

tt_rows <- lapply(ranks, function(r) {
  t0 <- proc.time()[[3L]]
  fit <- ttps(
    dat$y, dat$X,
    family = gaussian(),
    rank = r,
    k = k_paper,
    lambda = lambda_paper,
    control = ctrl_base(max_sweeps = 30L)
  )
  elapsed <- proc.time()[[3L]] - t0
  stopifnot(identical(fit$penalty_mode, "global"))
  data.frame(
    model = "TT",
    rank = r,
    npar = fit$npar_tt,
    edf = fit$edf,
    rmse_y = sqrt(mean((dat$y - fitted(fit))^2)),
    rmse_truth = sqrt(mean((dat$truth - fitted(fit))^2)),
    elapsed = elapsed,
    penalty_mode = fit$penalty_mode,
    fit = I(list(fit))
  )
})
tt_tab <- do.call(rbind, lapply(tt_rows, function(z) z[, setdiff(names(z), "fit")]))
print(tt_tab, digits = 4, row.names = FALSE)

t0 <- proc.time()[[3L]]
dense <- fit_dense_tp3(dat$y, dat$X, k = k_paper, lambda = lambda_paper)
dense$elapsed <- proc.time()[[3L]] - t0
dense_row <- data.frame(
  model = "dense Array",
  rank = NA_integer_,
  npar = dense$npar,
  edf = NA_real_,
  rmse_y = sqrt(mean((dat$y - dense$mu)^2)),
  rmse_truth = sqrt(mean((dat$truth - dense$mu)^2)),
  elapsed = dense$elapsed,
  penalty_mode = "global-dense"
)
print(dense_row, digits = 4, row.names = FALSE)

paper_tab <- rbind(tt_tab, dense_row)
write.csv(paper_tab, file.path(out_dir, "paper1_tt_rank_path_global.csv"),
          row.names = FALSE)
write.csv(dense_row, file.path(out_dir, "paper1_dense_baseline_global.csv"),
          row.names = FALSE)
# also overwrite legacy filenames used by Brain snapshots
write.csv(tt_tab[, c("rank", "npar", "rmse_y", "rmse_truth")],
          file.path(out_dir, "paper1_tt_rank_path.csv"), row.names = FALSE)
write.csv(dense_row[, c("model", "rank", "npar", "edf", "rmse_y", "rmse_truth",
                        "elapsed")],
          file.path(out_dir, "paper1_dense_baseline.csv"), row.names = FALSE)

demo_fit <- tt_rows[[3]]$fit[[1]]
mu_demo <- fitted(demo_fit)

png(file.path(fig_dir, "paper1_diagnostics.png"), width = 1100, height = 420, res = 120)
par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1))
plot(tt_tab$rank, tt_tab$rmse_truth, type = "b", pch = 16, col = "#1F4E79",
     xlab = "TT rank r", ylab = "RMSE (truth)",
     main = expression(paste("Global ", J[lambda], ": RMSE vs rank")),
     ylim = range(c(tt_tab$rmse_truth, dense_row$rmse_truth)))
abline(h = dense_row$rmse_truth, col = "grey40", lty = 2)
legend("topright", c("TT", "dense"), col = c("#1F4E79", "grey40"),
       lty = c(1, 2), pch = c(16, NA), bty = "n", cex = 0.85)
plot(tt_tab$rank, tt_tab$npar, type = "b", pch = 16, col = "#C0392B",
     xlab = "TT rank r", ylab = expression(N[TT]),
     main = "Stored TT parameters")
abline(h = dense_row$npar, col = "grey40", lty = 2)
plot(mu_demo, dat$y, pch = 16, cex = 0.4, col = adjustcolor("#1F4E79", 0.35),
     xlab = expression(hat(mu)), ylab = "y",
     main = sprintf("r=3 fitted vs y (cor=%.2f)", cor(mu_demo, dat$y)))
abline(0, 1, col = "grey30")
dev.off()

# surface slices at median of third coord
png(file.path(fig_dir, "paper1_surface.png"), width = 1100, height = 420, res = 120)
par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1))
g <- seq(0.05, 0.95, length.out = 40)
med3 <- median(dat$X[, 3])
grid12 <- as.matrix(expand.grid(g, g, med3))
truth12 <- matrix(true_surface(grid12[, 1], grid12[, 2], grid12[, 3]), 40, 40)
fit12 <- matrix(predict(demo_fit, grid12, type = "response"), 40, 40)
zlim <- range(c(truth12, fit12))
image(g, g, truth12, zlim = zlim, col = hcl.colors(40, "Blues"),
      xlab = "x1", ylab = "x2", main = "Truth (x3=median)")
image(g, g, fit12, zlim = zlim, col = hcl.colors(40, "Blues"),
      xlab = "x1", ylab = "x2", main = "TT r=3 fitted")
image(g, g, fit12 - truth12, col = hcl.colors(40, "RdBu", rev = TRUE),
      xlab = "x1", ylab = "x2", main = "Fitted − truth")
dev.off()
cat("Wrote paper1 figures + CSVs.\n\n")

########################################################################
# B. Chicago Poisson: CV rank + r=1 vs r=3 (global penalty)
########################################################################
cat("=== B. Chicago Poisson (global penalty) ===\n")
stopifnot(requireNamespace("gamair", quietly = TRUE))
library(gamair)
data(chicago)
dchi <- stats::na.omit(chicago[, c(
  "death", "tmpd", "o3median", "pm10median", "so2median", "time"
)])
y <- dchi$death
X <- as.matrix(dchi[, c("tmpd", "o3median", "pm10median", "so2median", "time")])

ctrl_cv <- ctrl_base(
  pirls_maxit = 25L,
  als_sweeps_per_pirls = 3L,
  lambda_bounds = c(1e-6, 1e6)
)

cat("Running tt_rank_select (ranks 1:5, folds=5, cGCV) ...\n")
t_sel <- system.time({
  sel <- tt_rank_select(
    y = y, X = X,
    ranks = 1:5,
    family = poisson(),
    k = 10L,
    lambda = "cGCV",
    folds = 5L,
    rule = "1se",
    control = ctrl_cv
  )
})
cat(sprintf("Selection wall time: %.1fs\n", t_sel[["elapsed"]]))
print(sel)

cat("Refitting selected rank on full data ...\n")
t_refit <- system.time({
  fit_sel <- tt_rank_refit(sel, control = ctrl_cv)
})
stopifnot(identical(fit_sel$penalty_mode, "global"))
cat(sprintf("Selected-rank refit: r=%d in %.2fs | penalty_mode=%s\n",
            fit_sel$rank_max, t_refit[["elapsed"]], fit_sel$penalty_mode))

saveRDS(
  list(selection = sel, fit = fit_sel,
       select_time = t_sel[["elapsed"]],
       refit_time = t_refit[["elapsed"]]),
  file.path(chi_dir, "fit_chicago_rank_cgcv.rds")
)

# Full-data refits r=1 and r=3
refit_rank <- function(r) {
  t0 <- proc.time()[[3L]]
  fit <- ttps(
    y, X,
    family = poisson(),
    rank = r,
    k = 10L,
    lambda = "cGCV",
    control = ctrl_cv
  )
  elapsed <- proc.time()[[3L]] - t0
  stopifnot(identical(fit$penalty_mode, "global"))
  list(fit = fit, time = elapsed)
}

cat("Refitting r=1 and r=3 full data ...\n")
r1 <- refit_rank(1L)
r3 <- refit_rank(3L)
f1 <- r1$fit; f3 <- r3$fit

tt_ic <- function(fit, criterion = c("AIC", "BIC")) {
  criterion <- match.arg(criterion)
  df <- fit$edf + 1
  pen <- if (identical(criterion, "AIC")) 2 else log(fit$n)
  fit$deviance + pen * df
}

cv_tab <- sel$cv_results
cv1 <- cv_tab$mean_cv[cv_tab$rank == 1]
cv3 <- cv_tab$mean_cv[cv_tab$rank == 3]
pr1 <- predict(f1, X, type = "response", interval = "confidence")
pr3 <- predict(f3, X, type = "response", interval = "confidence")

summary_tab <- data.frame(
  rank = c(1, 3),
  npar_tt = c(f1$npar_tt, f3$npar_tt),
  compression = c(f1$compression_ratio, f3$compression_ratio),
  deviance = c(f1$deviance, f3$deviance),
  edf = c(f1$edf, f3$edf),
  AIC = c(tt_ic(f1, "AIC"), tt_ic(f3, "AIC")),
  BIC = c(tt_ic(f1, "BIC"), tt_ic(f3, "BIC")),
  CV = c(cv1, cv3),
  cor = c(cor(y, fitted(f1)), cor(y, fitted(f3))),
  RMSE = c(sqrt(mean((y - fitted(f1))^2)), sqrt(mean((y - fitted(f3))^2))),
  mu_min = c(min(fitted(f1)), min(fitted(f3))),
  mu_max = c(max(fitted(f1)), max(fitted(f3))),
  mean_CIw = c(mean(pr1$upper - pr1$lower), mean(pr3$upper - pr3$lower)),
  fit_time_s = c(r1$time, r3$time),
  penalty_mode = c(f1$penalty_mode, f3$penalty_mode)
)
print(summary_tab, digits = 4)

lambda_tab <- data.frame(
  margin = colnames(X),
  lambda_r1 = signif(f1$lambda, 4),
  bound_r1 = f1$lambda_boundary,
  lambda_r3 = signif(f3$lambda, 4),
  bound_r3 = f3$lambda_boundary
)
print(lambda_tab)

write.csv(summary_tab, file.path(out_dir, "chicago_summary_global.csv"),
          row.names = FALSE)
write.csv(lambda_tab, file.path(out_dir, "chicago_lambda_global.csv"),
          row.names = FALSE)
write.csv(cv_tab, file.path(out_dir, "chicago_cv_rank_global.csv"),
          row.names = FALSE)

saveRDS(
  list(f1 = f1, f3 = f3, tab = summary_tab, lambda = lambda_tab, cv = cv_tab),
  file.path(chi_dir, "fit_chicago_r1_vs_r3.rds")
)

# ---- plots ----
copy_fig <- function(name) {
  file.copy(file.path(fig_dir, name), file.path(chi_dir, name), overwrite = TRUE)
}

png(file.path(fig_dir, "chicago_rank_cv.png"), width = 720, height = 480, res = 120)
plot(sel)
title(main = "Chicago CV rank path (global J_λ)", line = 2.2)
dev.off()
copy_fig("chicago_rank_cv.png")

t <- X[, "time"]; ord <- order(t)
png(file.path(fig_dir, "chicago_time_fit.png"), width = 1200, height = 520, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 2.8, 1))
plot_one <- function(pr, col, title) {
  plot(t[ord], y[ord], type = "l", col = adjustcolor("grey50", 0.7), lwd = 0.55,
       xlab = "time (days, gamair::chicago as stored)", ylab = "daily deaths",
       main = title, ylim = c(0, max(y) * 1.02))
  polygon(c(t[ord], rev(t[ord])),
          c(pmax(pr$lower[ord], 0), rev(pr$upper[ord])),
          col = adjustcolor(col, 0.22), border = NA)
  lines(t[ord], pr$fit[ord], col = col, lwd = 1.35)
  legend("topright", c("observed", "TT fitted", "pointwise 95% CI"),
         col = c("grey50", col, adjustcolor(col, 0.4)),
         lwd = c(1, 1.5, 8), bty = "n", cex = 0.85)
}
plot_one(pr1, "#1F4E79",
         sprintf("r=1, k=10; EDF=%.0f; cor=%.2f (global)", f1$edf, cor(y, pr1$fit)))
plot_one(pr3, "#C0392B",
         sprintf("r=3, k=10; EDF=%.0f; cor=%.2f (global)", f3$edf, cor(y, pr3$fit)))
dev.off()
copy_fig("chicago_time_fit.png")

# partial slices r1 vs r3
med <- apply(X, 2, median)
.slice <- function(fit, j, n = 120L) {
  xj <- X[, j]
  lohi <- stats::quantile(xj, c(0.025, 0.975), names = FALSE)
  g <- seq(lohi[1], lohi[2], length.out = n)
  Xm <- matrix(rep(med, each = length(g)), length(g), ncol(X), byrow = FALSE)
  colnames(Xm) <- colnames(X)
  Xm[, j] <- g
  pr <- predict(fit, Xm, type = "response", interval = "confidence")
  list(x = g, fit = pr$fit, lower = pr$lower, upper = pr$upper)
}
sl1 <- lapply(seq_len(5), function(j) .slice(f1, j))
sl3 <- lapply(seq_len(5), function(j) .slice(f3, j))
y_cap <- max(stats::quantile(y, 0.999),
             max(vapply(c(sl1, sl3), function(s) max(s$fit), 1)) * 1.25)
ylim <- c(0, y_cap)

png(file.path(fig_dir, "chicago_compare_r1_r3_slices.png"),
    width = 1100, height = 700, res = 120)
par(mfrow = c(2, 3), mar = c(4.2, 4, 2.5, 1))
for (j in 1:5) {
  a <- sl1[[j]]; b <- sl3[[j]]
  plot(a$x, a$fit, type = "n", ylim = ylim, xlim = range(a$x),
       xlab = colnames(X)[j], ylab = expression(hat(mu)),
       main = sprintf("slice vs %s", colnames(X)[j]))
  polygon(c(a$x, rev(a$x)), c(pmax(a$lower, 0), rev(pmin(a$upper, y_cap))),
          col = adjustcolor("#1F4E79", 0.2), border = NA)
  polygon(c(b$x, rev(b$x)), c(pmax(b$lower, 0), rev(pmin(b$upper, y_cap))),
          col = adjustcolor("#C0392B", 0.18), border = NA)
  lines(a$x, a$fit, col = "#1F4E79", lwd = 2)
  lines(b$x, b$fit, col = "#C0392B", lwd = 2)
  rug(X[, j], col = adjustcolor("grey40", 0.25), quiet = TRUE)
}
plot.new()
legend("center", c("r=1", "r=3"), col = c("#1F4E79", "#C0392B"),
       lwd = 2, bty = "n", cex = 1.2)
dev.off()
copy_fig("chicago_compare_r1_r3_slices.png")

# also refresh k=5 cGCV demo used by example_chicago_plots
cat("Refreshing k=5 cGCV Chicago demo fit ...\n")
fit_k5 <- ttps(
  y, X, family = poisson(), rank = 3L, k = 5L, lambda = "cGCV",
  control = ctrl_cv
)
saveRDS(fit_k5, file.path(chi_dir, "fit_chicago_cgcv.rds"))

cat("\n=== DONE (global penalty) ===\n")
cat("Paper-1 table:\n"); print(paper_tab[, c("model","rank","npar","rmse_y","rmse_truth")], digits = 3)
cat("\nChicago summary:\n"); print(summary_tab, digits = 3)
cat("\nChicago lambdas:\n"); print(lambda_tab)
invisible(list(paper = paper_tab, chicago = summary_tab, lambda = lambda_tab, sel = sel))
