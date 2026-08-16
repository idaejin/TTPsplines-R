# Reproducible example: Margin Activity Path + drop test
#
# Fixed seeds and knobs so re-runs match (ALS + R backend).
# From package root:
#   Rscript inst/examples/reprex_margin_activity_path.R
#
# Or paste the block below after:
#   pak::pak("idaejin/TTPsplines-R"); library(TTPsplines)

root <- (function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), "../..")))
  }
  normalizePath(".")
})()
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

out_dir <- file.path(root, "inst", "examples", "margin_activity_path_figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- reprex begin ----------------------------------------------------------

set.seed(20260817)
n <- 200L
d <- 6L
X <- matrix(runif(n * d), n, d)
colnames(X) <- paste0("x", seq_len(d))
# Truth: only x1 and x2 are active
y <- sin(2 * pi * X[, 1L]) + sin(2 * pi * X[, 2L]) + stats::rnorm(n, sd = 0.3)

set.seed(20260818)
Xte <- matrix(runif(n * d), n, d)
colnames(Xte) <- colnames(X)
yte <- sin(2 * pi * Xte[, 1L]) + sin(2 * pi * Xte[, 2L]) +
  stats::rnorm(n, sd = 0.3)

ctrl <- tt_control(
  max_sweeps = 5L,
  outer_maxit = 5L,
  compute_edf = TRUE,
  backend = "R",
  seed = 20260817
)

mse <- function(fit, Xnew, ynew) {
  mean((as.numeric(predict(fit, Xnew)) - as.numeric(ynew))^2)
}

# 1) Full model
fit_full <- ttps(y, X, rank = 2L, k = 5L, lambda = "cGCV", control = ctrl)

# 2) Margin Activity Path (formal selector)
path <- tt_margin_activity_path(
  y, X,
  rank = 2L,
  k = 5L,
  lambda_path = 10^c(1, 0.5, 0, -0.5),
  select = "1se",
  folds = 4L,
  seed = 20260817,
  select_lambda = "cGCV",
  refit = TRUE,
  control = tt_control(
    max_sweeps = 5L, outer_maxit = 5L, compute_edf = FALSE,
    backend = "R", seed = 20260817
  ),
  verbose = FALSE
)

# 3) Leave-one-margin drop test (diagnostic)
tst <- tt_margin_drop_test(
  y, X,
  rank = 2L,
  k = 5L,
  lambda = 1,
  method = "nested",
  alpha = 0.05,
  seed = 20260817,
  control = ctrl,
  verbose = FALSE
)

res <- data.frame(
  method = c("full_cGCV", "activity_path_1se", "drop_test_keep"),
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
        rank = 2L, k = 5L, lambda = "cGCV",
        control = tt_control(
          max_sweeps = 5L, outer_maxit = 5L, compute_edf = FALSE,
          backend = "R", seed = 20260817
        )
      )
      mse(fit_k, Xte[, tst$keep, drop = FALSE], yte)
    }
  ),
  stringsAsFactors = FALSE
)

## ---- reprex end ------------------------------------------------------------

cat("Reproducible Margin Activity Path example\n")
cat("Truth active set: x1, x2\n\n")
print(res)
cat("\nActivity path selected: ", paste(path$selected_names, collapse = ", "),
    "\n", sep = "")
cat("Drop-test keep:         ", paste(tst$keep_names, collapse = ", "),
    "\n", sep = "")
cat("Drop-test drop:         ",
    if (length(tst$drop_candidate_names)) {
      paste(tst$drop_candidate_names, collapse = ", ")
    } else {
      "(none)"
    },
    "\n", sep = "")
cat("\nEDF (full model):\n")
print(tt_edf(fit_full))

utils::write.csv(res, file.path(out_dir, "reprex_holdout.csv"), row.names = FALSE)
utils::write.csv(
  path$cv,
  file.path(out_dir, "reprex_cv_path.csv"),
  row.names = FALSE
)
utils::write.csv(
  tst$results,
  file.path(out_dir, "reprex_drop_test.csv"),
  row.names = FALSE
)

png(file.path(out_dir, "reprex_activity_path.png"),
    width = 1000, height = 360, res = 120)
plot(path)
dev.off()

cat("\nWrote:\n  ", file.path(out_dir, "reprex_holdout.csv"), "\n  ",
    file.path(out_dir, "reprex_cv_path.csv"), "\n  ",
    file.path(out_dir, "reprex_drop_test.csv"), "\n  ",
    file.path(out_dir, "reprex_activity_path.png"), "\n", sep = "")
