# Margin Activity Path + margin drop test
#
# Screens which covariates to keep for a TT–P-spline, then (optionally)
# runs leave-one-margin-out drop diagnostics and reports EDF by margin.
#
# Usage (from package root):
#   Rscript inst/examples/example_margin_activity_path.R

root <- (function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), "../..")))
  }
  if (requireNamespace("pkgload", quietly = TRUE) ||
      requireNamespace("devtools", quietly = TRUE)) {
    return(normalizePath("."))
  }
  normalizePath(".")
})()

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, quiet = TRUE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

fig_dir <- file.path(root, "inst", "examples", "margin_activity_path_figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Data: additive sins on x1..x3; x4..x8 null
# ---------------------------------------------------------------------------

set.seed(1)
n <- 250
d <- 8
n_signal <- 3L
X <- matrix(runif(n * d), n, d)
colnames(X) <- paste0("x", seq_len(d))
y <- rowSums(sin(2 * pi * X[, seq_len(n_signal), drop = FALSE])) +
  stats::rnorm(n, sd = 0.3)

set.seed(2)
Xte <- matrix(runif(n * d), n, d)
colnames(Xte) <- colnames(X)
yte <- rowSums(sin(2 * pi * Xte[, seq_len(n_signal), drop = FALSE])) +
  stats::rnorm(n, sd = 0.3)

ctrl_screen <- tt_control(
  max_sweeps = 6L, compute_edf = FALSE, outer_maxit = 6L, backend = "R"
)
ctrl_edf <- tt_control(
  max_sweeps = 6L, compute_edf = TRUE, outer_maxit = 6L, backend = "R"
)

mse <- function(fit, Xnew, ynew) {
  mean((as.numeric(predict(fit, Xnew)) - as.numeric(ynew))^2)
}

cat("\n=== 1. Full TT + cGCV (all margins) ===\n")
fit_full <- ttps(
  y, X, rank = 2L, k = 5L, lambda = "cGCV", control = ctrl_edf
)
print(fit_full)
cat("Per-margin EDF:\n")
print(round(tt_edf(fit_full)$margin, 3))
cat(sprintf("Holdout MSE (full): %.4f\n", mse(fit_full, Xte, yte)))

# ---------------------------------------------------------------------------
# Margin Activity Path (formal CV / 1-SE selection)
# ---------------------------------------------------------------------------

cat("\n=== 2. Margin Activity Path (select = 1se) ===\n")
path <- tt_margin_activity_path(
  y, X,
  rank = 2L,
  k = 5L,
  lambda_path = 10^seq(1.5, -0.5, by = -0.5),
  select = "1se",
  folds = 5L,
  seed = 1,
  select_lambda = "cGCV",
  refit = TRUE,
  control = ctrl_screen,
  verbose = TRUE
)
print(path)
cat("Selected margins:", paste(path$selected_names, collapse = ", "), "\n")
cat(sprintf(
  "Holdout MSE (activity path): %.4f\n",
  mse(path$fit, Xte[, path$selected, drop = FALSE], yte)
))

png(file.path(fig_dir, "01_activity_path_panels.png"),
    width = 1100, height = 380, res = 120)
plot(path)
dev.off()
cat("Wrote ", file.path(fig_dir, "01_activity_path_panels.png"), "\n", sep = "")

# ---------------------------------------------------------------------------
# Leave-one-margin drop test (diagnostic; often more liberal than 1-SE)
# ---------------------------------------------------------------------------

cat("\n=== 3. Margin drop test (nested, approximate F) ===\n")
tst <- tt_margin_drop_test(
  y, X,
  rank = 2L,
  k = 5L,
  lambda = 1,
  method = "nested",
  alpha = 0.05,
  control = ctrl_edf,
  verbose = TRUE
)
print(tst)

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

cat("\n=== Summary ===\n")
truth <- colnames(X)[seq_len(n_signal)]
tab <- data.frame(
  method = c("Full TT+cGCV", "Margin Activity Path (1-SE)", "Drop-test keep"),
  margins = c(
    paste(colnames(X), collapse = ","),
    paste(path$selected_names, collapse = ","),
    paste(tst$keep_names, collapse = ",")
  ),
  n_margins = c(d, length(path$selected), length(tst$keep)),
  test_MSE = c(
    mse(fit_full, Xte, yte),
    mse(path$fit, Xte[, path$selected, drop = FALSE], yte),
    {
      fit_k <- ttps(
        y, X[, tst$keep, drop = FALSE],
        rank = 2L, k = 5L, lambda = "cGCV", control = ctrl_screen
      )
      mse(fit_k, Xte[, tst$keep, drop = FALSE], yte)
    }
  ),
  stringsAsFactors = FALSE
)
print(tab)
cat("Truth active set:", paste(truth, collapse = ", "), "\n")
cat("Figures: ", fig_dir, "\n", sep = "")
cat("Done.\n")
