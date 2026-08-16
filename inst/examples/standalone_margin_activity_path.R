## Standalone example: Margin Activity Path + drop test (p-values)
##
## Run from anywhere after installing the package, e.g.:
##
##   install.packages("pak")
##   pak::pak("idaejin/TTPsplines-R")   # or remotes::install_github(...)
##
##   Rscript path/to/standalone_margin_activity_path.R
##
## Or open this file in RStudio and Source it.

suppressPackageStartupMessages(library(TTPsplines))

if (!exists("tt_margin_activity_path", mode = "function")) {
  stop(
    "Your installed TTPsplines is too old (missing tt_margin_activity_path).\n",
    "Reinstall from GitHub: pak::pak(\"idaejin/TTPsplines-R\")",
    call. = FALSE
  )
}

## ---- data (reproducible) ---------------------------------------------------

set.seed(20260817)
n <- 200L
d <- 6L
X <- matrix(runif(n * d), n, d)
colnames(X) <- paste0("x", seq_len(d))

# Truth: only x1 and x2 matter
y <- sin(2 * pi * X[, 1L]) + sin(2 * pi * X[, 2L]) + rnorm(n, sd = 0.3)

set.seed(20260818)
Xte <- matrix(runif(n * d), n, d)
colnames(Xte) <- colnames(X)
yte <- sin(2 * pi * Xte[, 1L]) + sin(2 * pi * Xte[, 2L]) + rnorm(n, sd = 0.3)

ctrl <- tt_control(
  max_sweeps = 5L,
  outer_maxit = 5L,
  compute_edf = TRUE,
  backend = "R",
  seed = 20260817L
)

mse <- function(fit, Xnew, ynew) {
  mean((as.numeric(predict(fit, Xnew)) - as.numeric(ynew))^2)
}

## ---- 1) full TT + cGCV -----------------------------------------------------

cat("\n=== Full TT + cGCV (all margins) ===\n")
fit_full <- ttps(y, X, rank = 2L, k = 5L, lambda = "cGCV", control = ctrl)
print(fit_full)
cat("EDF by margin:\n")
print(round(tt_edf(fit_full)$margin, 3))
cat(sprintf("Holdout MSE: %.4f\n", mse(fit_full, Xte, yte)))

## ---- 2) Margin Activity Path (CV 1-SE) -------------------------------------

cat("\n=== Margin Activity Path (select = 1se) ===\n")
path <- tt_margin_activity_path(
  y, X,
  rank = 2L,
  k = 5L,
  lambda_path = 10^c(1, 0.5, 0, -0.5),
  select = "1se",
  folds = 4L,
  seed = 20260817L,
  select_lambda = "cGCV",
  refit = TRUE,
  control = tt_control(
    max_sweeps = 5L, outer_maxit = 5L, compute_edf = FALSE,
    backend = "R", seed = 20260817L
  ),
  verbose = TRUE
)
print(path)
cat("Selected:", paste(path$selected_names, collapse = ", "), "\n")
cat(sprintf(
  "Holdout MSE: %.4f\n",
  mse(path$fit, Xte[, path$selected, drop = FALSE], yte)
))

## Optional plot (comment out if running headless without a device)
if (interactive() || capabilities("png")) {
  plot(path)
}

## ---- 3) Statistical drop test with p-values --------------------------------

cat("\n=== Margin drop test (nested F approx. -> p-values) ===\n")
tst <- tt_margin_drop_test(
  y, X,
  rank = 2L,
  k = 5L,
  lambda = 1,
  method = "nested",   # use "permute" + B = 99 for permutation p-values
  alpha = 0.05,
  seed = 20260817L,
  control = ctrl,
  verbose = TRUE
)
print(tst)

cat("\n=== p-values by margin ===\n")
print(tst$results[, c("margin", "delta_deviance", "stat", "p_value",
                      "drop_candidate")])

## ---- summary ---------------------------------------------------------------

cat("\n=== Summary ===\n")
cat("Truth active set:     x1, x2\n")
cat("Activity path keep:   ", paste(path$selected_names, collapse = ", "), "\n",
    sep = "")
cat("Drop-test keep:       ", paste(tst$keep_names, collapse = ", "), "\n",
    sep = "")
cat("Drop-test candidates: ",
    if (length(tst$drop_candidate_names)) {
      paste(tst$drop_candidate_names, collapse = ", ")
    } else {
      "(none)"
    },
    "\n", sep = "")
cat(sprintf(
  "Holdout MSE full / path: %.4f / %.4f\n",
  mse(fit_full, Xte, yte),
  mse(path$fit, Xte[, path$selected, drop = FALSE], yte)
))
cat("Done.\n")
