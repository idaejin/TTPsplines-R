# Margin Activity Path — demo (see PROTOCOL_MAP.md)
#
# Prefer package API: tt_margin_activity_path() + vignette("margin-activity-path")
#
# Canonical formalization:
#   inst/benchmarks/margin_path/PROTOCOL_MAP.md
#
# This script is an extended simulation demo (+ optional linear lasso baseline).
# Pipeline:
#   1) Isotropic-lambda path; A_j(lambda) = partial range; S_j = mean_lambda A_j
#   2) Rank margins by S_j; nested top-m models
#   3) Choose m by K-fold CV (default: 1-SE rule)
#   4) Refit TT+cGCV on A_hat; compare full-d cGCV + optional linear lasso
#
# Optional: install.packages("glmnet")
# Usage:
#   Rscript inst/benchmarks/margin_path/run_map_example.R

root <- Sys.getenv("TT_PKG_ROOT", unset = "")
if (!nzchar(root)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
    root <- normalizePath(file.path(dirname(script), "..", "..", ".."))
  } else {
    root <- getwd()
  }
}
setwd(root)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, export_all = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

out_dir <- file.path(root, "inst", "benchmarks", "margin_path", "results")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# Config
# ===========================================================================

set.seed(1)
n <- 500L
d <- 10L
n_signal <- 5L          # true active margins: x1..x_{n_signal}
rank_tt <- 2L
k_basis <- 5L
K_folds <- 5L
use_one_se <- TRUE      # FALSE = plain argmin CV
lasso_s <- "lambda.1se" # or "lambda.min" (glmnet CV pick)
lams <- 10^seq(2, -1, by = -0.5)

ctrl <- tt_control(max_sweeps = 8L, backend = "R", compute_edf = FALSE)
ctrl_cgcv <- tt_control(
  max_sweeps = 8L, backend = "R", compute_edf = FALSE, outer_maxit = 8L
)

# ===========================================================================
# Data: train + holdout
# ===========================================================================

X <- matrix(runif(n * d), n, d)
colnames(X) <- paste0("x", seq_len(d))
y <- rowSums(sin(2 * pi * X[, seq_len(n_signal), drop = FALSE])) +
  stats::rnorm(n, sd = 0.3)

set.seed(2)
Xte <- matrix(runif(n * d), n, d)
colnames(Xte) <- colnames(X)
yte <- rowSums(sin(2 * pi * Xte[, seq_len(n_signal), drop = FALSE])) +
  stats::rnorm(n, sd = 0.3)

# ===========================================================================
# Helpers
# ===========================================================================

partial_1d <- function(fit, X, j, n_grid = 80L) {
  med <- apply(X, 2L, stats::median)
  g <- seq(
    stats::quantile(X[, j], 0.025),
    stats::quantile(X[, j], 0.975),
    length.out = n_grid
  )
  Xm <- matrix(med, n_grid, ncol(X), byrow = TRUE)
  colnames(Xm) <- colnames(X)
  Xm[, j] <- g
  data.frame(x = as.numeric(g), fhat = as.numeric(predict(fit, Xm)))
}

partial_range <- function(fit, X, j) {
  med <- apply(X, 2L, stats::median)
  g <- seq(0.05, 0.95, length.out = 40L)
  Xm <- matrix(med, length(g), ncol(X), byrow = TRUE)
  colnames(Xm) <- colnames(X)
  Xm[, j] <- g
  diff(range(as.numeric(predict(fit, Xm))))
}

fit_ttps <- function(y, X, lambda, control) {
  p <- ncol(X)
  if (p < 1L) {
    mu <- mean(y)
    return(list(
      type = "intercept",
      mu = mu,
      predict = function(Xnew) rep(mu, nrow(Xnew))
    ))
  }
  if (p == 1L) {
    # ttps requires d >= 2; univariate cubic smooth with GCV
    ss <- stats::smooth.spline(X[, 1L], y, cv = NA)
    return(list(
      type = "univariate",
      ss = ss,
      predict = function(Xnew) as.numeric(stats::predict(ss, x = Xnew[, 1L])$y)
    ))
  }
  ttps(
    y, X,
    rank = min(rank_tt, max(1L, p)),
    k = k_basis,
    lambda = lambda,
    control = control
  )
}

predict_any <- function(fit, Xnew) {
  if (is.list(fit) && !is.null(fit$type) &&
      fit$type %in% c("intercept", "univariate")) {
    return(fit$predict(Xnew))
  }
  as.numeric(predict(fit, Xnew))
}

mse_pred <- function(fit, Xnew, ynew) {
  mean((predict_any(fit, Xnew) - as.numeric(ynew))^2)
}

# K-fold CV MSE for a given column subset (indices into X)
cv_mse_subset <- function(y, X, cols, folds, lambda = "cGCV") {
  K <- max(folds)
  pe <- numeric(K)
  for (k in seq_len(K)) {
    tr <- folds != k
    te <- !tr
    Xtr <- X[tr, cols, drop = FALSE]
    Xte_k <- X[te, cols, drop = FALSE]
    if (length(cols) == 0L) {
      mu <- mean(y[tr])
      pe[k] <- mean((mu - y[te])^2)
    } else {
      fi <- fit_ttps(y[tr], Xtr, lambda = lambda, control = ctrl_cgcv)
      pe[k] <- mean((predict_any(fi, Xte_k) - y[te])^2)
    }
  }
  list(mse = mean(pe), se = stats::sd(pe) / sqrt(K), fold_mse = pe)
}

# ===========================================================================
# Step 1–2: MAP scores S_j
# ===========================================================================

message("=== MAP path (isotropic lambda) ===")
score <- matrix(NA_real_, length(lams), d)
colnames(score) <- colnames(X)

for (i in seq_along(lams)) {
  fi <- fit_ttps(y, X, lambda = lams[i], control = ctrl)
  for (j in seq_len(d)) score[i, j] <- partial_range(fi, X, j)
  message(sprintf("  lambda=%.4g  done", lams[i]))
}

S <- colMeans(score)
ord <- order(S, decreasing = TRUE)
S_ord <- S[ord]
message("Scores (sorted):")
print(round(S_ord, 4))

# Diagnostic: largest-gap cut (not used for final selection)
gaps <- c(S_ord[-length(S_ord)] - S_ord[-1], S_ord[length(S_ord)])
m_gap <- which.max(gaps)
A_gap <- ord[seq_len(m_gap)]
message(sprintf(
  "Gap rule would keep m=%d: %s",
  m_gap, paste(colnames(X)[A_gap], collapse = ", ")
))

# ===========================================================================
# Step 3–4: CV over nested top-m models
# ===========================================================================

message("=== CV over ranked nested models ===")
set.seed(3)
folds <- sample(rep(seq_len(K_folds), length.out = n))

# m = 0..d  (0 = intercept only)
ms <- 0:d
cv_tab <- data.frame(
  m = ms,
  cv_mse = NA_real_,
  cv_se = NA_real_,
  margins = NA_character_,
  stringsAsFactors = FALSE
)

for (ii in seq_along(ms)) {
  m <- ms[ii]
  cols <- if (m == 0L) integer(0) else ord[seq_len(m)]
  res <- cv_mse_subset(y, X, cols, folds, lambda = "cGCV")
  cv_tab$cv_mse[ii] <- res$mse
  cv_tab$cv_se[ii] <- res$se
  cv_tab$margins[ii] <- if (m == 0L) "(none)" else
    paste(colnames(X)[cols], collapse = ",")
  message(sprintf(
    "  m=%2d  CV-MSE=%.4f  (SE=%.4f)  {%s}",
    m, res$mse, res$se, cv_tab$margins[ii]
  ))
}

m_min <- cv_tab$m[which.min(cv_tab$cv_mse)]
cv_star <- cv_tab$cv_mse[cv_tab$m == m_min]
se_star <- cv_tab$cv_se[cv_tab$m == m_min]

if (isTRUE(use_one_se)) {
  # smallest m with CV within 1 SE of the minimum
  ok <- which(cv_tab$cv_mse <= cv_star + se_star)
  m_hat <- min(cv_tab$m[ok])
  rule <- "1-SE"
} else {
  m_hat <- m_min
  rule <- "min-CV"
}

A_hat <- if (m_hat == 0L) integer(0) else ord[seq_len(m_hat)]
message(sprintf(
  "Selected m=%d by %s: %s",
  m_hat, rule,
  if (length(A_hat)) paste(colnames(X)[A_hat], collapse = ", ") else "(none)"
))

# ===========================================================================
# Step 5b: Classical linear lasso (glmnet)
# ===========================================================================
# Linear model y ~ X with L1 penalty; CV chooses lambda.
# Note: DGP here is nonlinear (sin), so lasso is a weak but standard baseline.

have_glmnet <- requireNamespace("glmnet", quietly = TRUE)
A_lasso <- integer(0)
coef_lasso <- setNames(rep(0, d), colnames(X))
train_lasso <- test_lasso <- NA_real_
n_lasso <- NA_integer_

if (have_glmnet) {
  message("=== Classical lasso (glmnet) ===")
  set.seed(4)
  cv_lasso <- glmnet::cv.glmnet(
    x = X, y = y, alpha = 1, nfolds = K_folds,
    standardize = TRUE, intercept = TRUE
  )
  lam_lasso <- cv_lasso[[lasso_s]]
  beta <- as.numeric(stats::coef(cv_lasso, s = lasso_s))[-1L]  # drop intercept
  names(beta) <- colnames(X)
  coef_lasso <- beta
  A_lasso <- which(abs(beta) > 0)
  n_lasso <- length(A_lasso)
  pred_tr <- as.numeric(stats::predict(cv_lasso, newx = X, s = lasso_s))
  pred_te <- as.numeric(stats::predict(cv_lasso, newx = Xte, s = lasso_s))
  train_lasso <- mean((pred_tr - y)^2)
  test_lasso <- mean((pred_te - yte)^2)
  message(sprintf(
    "Lasso %s: lambda=%.4g  selected %d: %s",
    lasso_s, lam_lasso, n_lasso,
    if (n_lasso) paste(colnames(X)[A_lasso], collapse = ", ") else "(none)"
  ))
  utils::write.csv(
    data.frame(
      margin = colnames(X),
      coef = as.numeric(coef_lasso),
      selected = seq_len(d) %in% A_lasso
    ),
    file.path(out_dir, "lasso_coefs.csv"),
    row.names = FALSE
  )
} else {
  message("Skipping lasso: install.packages(\"glmnet\") to enable baseline.")
}

# ===========================================================================
# Step 5: Final fits + holdout
# ===========================================================================

message("=== Final fits ===")
fit_full <- fit_ttps(y, X, lambda = "cGCV", control = ctrl_cgcv)
fit_map <- fit_ttps(
  y, X[, A_hat, drop = FALSE],
  lambda = "cGCV", control = ctrl_cgcv
)
fit_screen <- fit_ttps(y, X, lambda = 1, control = ctrl)

tab <- data.frame(
  method = c(
    "TT+cGCV full",
    sprintf("MAP+CV(%s) then cGCV", rule),
    sprintf("Lasso linear (%s)", lasso_s)
  ),
  n_margins = c(d, length(A_hat), n_lasso),
  train_MSE = c(
    mse_pred(fit_full, X, y),
    mse_pred(fit_map, X[, A_hat, drop = FALSE], y),
    train_lasso
  ),
  test_MSE = c(
    mse_pred(fit_full, Xte, yte),
    mse_pred(fit_map, Xte[, A_hat, drop = FALSE], yte),
    test_lasso
  ),
  stringsAsFactors = FALSE
)
print(tab)

# Truth recovery
truth <- seq_len(n_signal)
tpr <- mean(truth %in% A_hat)
fpr <- mean(setdiff(seq_len(d), truth) %in% A_hat)
message(sprintf(
  "MAP selection TPR=%.2f  FPR=%.2f  (truth x1..x%d)",
  tpr, fpr, n_signal
))
if (have_glmnet) {
  tpr_l <- mean(truth %in% A_lasso)
  fpr_l <- mean(setdiff(seq_len(d), truth) %in% A_lasso)
  message(sprintf("Lasso selection TPR=%.2f  FPR=%.2f", tpr_l, fpr_l))
}

# Save tables
utils::write.csv(cv_tab, file.path(out_dir, "cv_path.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(
    margin = names(S), score = as.numeric(S),
    rank = match(seq_len(d), ord),
    active_map = seq_len(d) %in% A_hat,
    active_lasso = if (have_glmnet) seq_len(d) %in% A_lasso else NA
  ),
  file.path(out_dir, "margin_scores.csv"),
  row.names = FALSE
)
utils::write.csv(
  cbind(log10_lambda = log10(lams), score),
  file.path(out_dir, "activity_path.csv"),
  row.names = FALSE
)
utils::write.csv(tab, file.path(out_dir, "holdout_comparison.csv"), row.names = FALSE)

# ===========================================================================
# Plots
# ===========================================================================

cols_active <- ifelse(seq_len(d) %in% A_hat, "#1F4E79", "#BBBBBB")

# 1) Partial plots (screen fit)
png(file.path(fig_dir, "01_partial_plots.png"), width = 1400, height = 700, res = 120)
op <- par(mfrow = c(2, 5), mar = c(4, 4, 2.2, 0.6))
yl <- range(unlist(lapply(seq_len(d), function(j) partial_1d(fit_screen, X, j)$fhat)))
for (j in seq_len(d)) {
  sl <- partial_1d(fit_screen, X, j)
  plot(sl$x, sl$fhat, type = "l", lwd = 2, col = cols_active[j], ylim = yl,
       xlab = colnames(X)[j], ylab = "partial f",
       main = sprintf("%s%s", colnames(X)[j],
                      if (j %in% A_hat) " [active]" else ""))
  rug(X[, j], col = adjustcolor("grey40", 0.3), quiet = TRUE)
}
par(op)
dev.off()

# 2) Activity path
png(file.path(fig_dir, "02_activity_path.png"), width = 900, height = 550, res = 120)
op <- par(mar = c(4.5, 4.5, 3, 1))
matplot(log10(lams), score, type = "l", lwd = 2, lty = 1, col = cols_active,
        xlab = expression(log[10](lambda)), ylab = "partial range A_j(lambda)",
        main = "Margin activity path")
legend("topright", legend = colnames(X), col = cols_active, lty = 1, lwd = 2,
       cex = 0.7, bty = "n", ncol = 2)
par(op)
dev.off()

# 3) Scores + gap diagnostic
png(file.path(fig_dir, "03_margin_scores.png"), width = 900, height = 450, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
barplot(S[ord], names.arg = names(S)[ord], col = cols_active[ord], border = NA,
        ylab = expression(S[j]), main = "B. Scores (sorted)", las = 2, cex.names = 0.8)
abline(h = S_ord[m_gap], lty = 3, col = "darkorange")
plot(ms, cv_tab$cv_mse, type = "b", lwd = 2, pch = 16, col = "#1F4E79",
     xlab = "m (top-m margins)", ylab = "CV-MSE",
     main = sprintf("C. CV path (%s)", rule))
arrows(ms, cv_tab$cv_mse - cv_tab$cv_se, ms, cv_tab$cv_mse + cv_tab$cv_se,
       length = 0.05, angle = 90, code = 3, col = adjustcolor("#1F4E79", 0.6))
abline(v = m_hat, lty = 2, col = "#C0392B")
abline(h = cv_star + se_star, lty = 3, col = "grey40")
par(op)
dev.off()

# 4) Holdout comparison (full / MAP / lasso)
hold_labs <- c("full cGCV", "MAP+CV", "Lasso")
png(file.path(fig_dir, "04_holdout_mse.png"), width = 800, height = 450, res = 120)
op <- par(mar = c(4.5, 4.5, 3, 1))
barplot(
  rbind(tab$train_MSE, tab$test_MSE),
  beside = TRUE,
  names.arg = hold_labs,
  col = c("#7F8C8D", "#C0392B"),
  ylab = "MSE",
  main = "Train vs holdout test MSE",
  legend.text = c("train", "test"),
  args.legend = list(x = "topright", bty = "n")
)
par(op)
dev.off()

# 4b) Lasso coefficient path (if available)
if (have_glmnet) {
  png(file.path(fig_dir, "04b_lasso_coefs.png"), width = 800, height = 450, res = 120)
  op <- par(mar = c(4.5, 4.5, 3, 1))
  cols_l <- ifelse(seq_len(d) %in% A_lasso, "#1F4E79", "#BBBBBB")
  barplot(
    coef_lasso[ord], names.arg = names(coef_lasso)[ord],
    col = cols_l[ord], border = NA, las = 2, cex.names = 0.8,
    ylab = "lasso coefficient",
    main = sprintf("Classical lasso coefs (%s)", lasso_s)
  )
  abline(h = 0, col = "grey50")
  par(op)
  dev.off()
}

# 5) Summary panel
png(file.path(fig_dir, "05_map_summary.png"), width = 1200, height = 900, res = 120)
op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

matplot(log10(lams), score, type = "l", lwd = 2, lty = 1, col = cols_active,
        xlab = expression(log[10](lambda)), ylab = "partial range",
        main = "A. Activity path")

barplot(S[ord], names.arg = names(S)[ord], col = cols_active[ord], border = NA,
        ylab = expression(S[j]), main = "B. Sorted scores", las = 2, cex.names = 0.7)

plot(ms, cv_tab$cv_mse, type = "b", lwd = 2, pch = 16, col = "#1F4E79",
     xlab = "m", ylab = "CV-MSE", main = paste0("C. CV selects m=", m_hat))
abline(v = m_hat, lty = 2, col = "#C0392B")

barplot(
  rbind(tab$train_MSE, tab$test_MSE),
  beside = TRUE,
  names.arg = hold_labs,
  col = c("#7F8C8D", "#C0392B"),
  ylab = "MSE",
  main = "D. Holdout",
  legend.text = c("train", "test"),
  args.legend = list(x = "topright", bty = "n", cex = 0.8)
)
par(op)
dev.off()

message("Figures: ", fig_dir)
message("Tables:  ", out_dir)
message("Done.")
