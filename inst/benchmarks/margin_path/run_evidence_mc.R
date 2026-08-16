# Mini MC evidence gate — Margin Activity Path
#
# B seeds x 2 DGPs (additive2, interaction).
# Usage (package root):
#   Rscript inst/benchmarks/margin_path/run_evidence_mc.R

root <- (function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), "../../..")))
  }
  normalizePath(".")
})()
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

out_dir <- file.path(root, "inst", "benchmarks", "margin_path", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

B <- 20L
n <- 200L
d <- 6L
truth <- c(1L, 2L)

ctrl_path <- tt_control(
  max_sweeps = 4L, outer_maxit = 4L, compute_edf = FALSE, backend = "R"
)
ctrl_fit <- tt_control(
  max_sweeps = 5L, outer_maxit = 5L, compute_edf = FALSE, backend = "R"
)

mse <- function(fit, Xnew, ynew) {
  mean((as.numeric(predict(fit, Xnew)) - as.numeric(ynew))^2)
}

sim_one <- function(seed, dgp = c("additive2", "interaction")) {
  dgp <- match.arg(dgp)
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  colnames(X) <- paste0("x", seq_len(d))
  f <- switch(
    dgp,
    additive2 = sin(2 * pi * X[, 1]) + sin(2 * pi * X[, 2]),
    interaction = sin(2 * pi * X[, 1] * X[, 2])
  )
  y <- f + stats::rnorm(n, sd = 0.3)

  set.seed(seed + 10000L)
  Xte <- matrix(runif(n * d), n, d)
  colnames(Xte) <- colnames(X)
  fte <- switch(
    dgp,
    additive2 = sin(2 * pi * Xte[, 1]) + sin(2 * pi * Xte[, 2]),
    interaction = sin(2 * pi * Xte[, 1] * Xte[, 2])
  )
  yte <- fte + stats::rnorm(n, sd = 0.3)

  fit_full <- ttps(y, X, rank = 2L, k = 5L, lambda = "cGCV", control = ctrl_fit)
  path <- tt_margin_activity_path(
    y, X, rank = 2L, k = 5L,
    lambda_path = 10^c(1, 0.5, 0, -0.5),
    select = "1se", folds = 4L, seed = seed,
    select_lambda = "cGCV", refit = TRUE,
    control = ctrl_path, verbose = FALSE
  )
  tst <- tt_margin_drop_test(
    y, X, rank = 2L, k = 5L, lambda = 1, method = "nested",
    alpha = 0.05, seed = seed,
    control = tt_control(
      max_sweeps = 4L, compute_edf = TRUE, backend = "R", seed = seed
    ),
    verbose = FALSE
  )

  sel <- as.integer(path$selected)
  keep <- as.integer(tst$keep)
  tpr <- function(A) if (!length(A)) 0 else mean(truth %in% A)
  fpr <- function(A) mean(setdiff(seq_len(d), truth) %in% A)

  mse_drop <- NA_real_
  if (length(keep) >= 2L) {
    fk <- ttps(
      y, X[, keep, drop = FALSE],
      rank = 2L, k = 5L, lambda = "cGCV", control = ctrl_fit
    )
    mse_drop <- mse(fk, Xte[, keep, drop = FALSE], yte)
  } else if (length(keep) == 1L) {
    ss <- stats::smooth.spline(X[, keep], y, cv = NA)
    mse_drop <- mean((stats::predict(ss, x = Xte[, keep])$y - yte)^2)
  }

  data.frame(
    seed = seed,
    dgp = dgp,
    m_hat = length(sel),
    exact_set = identical(sort(sel), sort(truth)),
    tpr_path = tpr(sel),
    fpr_path = fpr(sel),
    tpr_drop = tpr(keep),
    fpr_drop = fpr(keep),
    mse_full = mse(fit_full, Xte, yte),
    mse_path = mse(path$fit, Xte[, sel, drop = FALSE], yte),
    mse_drop = mse_drop,
    selected = paste(colnames(X)[sel], collapse = ","),
    drop_keep = paste(colnames(X)[keep], collapse = ","),
    stringsAsFactors = FALSE
  )
}

rows <- list()
ii <- 1L
for (dgp in c("additive2", "interaction")) {
  for (b in seq_len(B)) {
    seed <- 1000L + b
    message(sprintf("[%s] seed=%d (%d/%d)", dgp, seed, b, B))
    rows[[ii]] <- sim_one(seed, dgp)
    ii <- ii + 1L
  }
}

tab <- do.call(rbind, rows)
utils::write.csv(tab, file.path(out_dir, "evidence_mc_raw.csv"), row.names = FALSE)

summary_dgp <- function(df) {
  data.frame(
    dgp = df$dgp[1],
    B = nrow(df),
    exact_set_rate = mean(df$exact_set),
    mean_tpr_path = mean(df$tpr_path),
    mean_fpr_path = mean(df$fpr_path),
    mean_tpr_drop = mean(df$tpr_drop),
    mean_fpr_drop = mean(df$fpr_drop),
    mean_mse_full = mean(df$mse_full),
    mean_mse_path = mean(df$mse_path),
    mean_mse_drop = mean(df$mse_drop, na.rm = TRUE),
    mean_mse_ratio_path = mean(df$mse_path / df$mse_full),
    stringsAsFactors = FALSE
  )
}

summ <- do.call(rbind, lapply(split(tab, tab$dgp), summary_dgp))
rownames(summ) <- NULL
utils::write.csv(summ, file.path(out_dir, "evidence_mc_summary.csv"), row.names = FALSE)

print(summ)
message("Wrote ", file.path(out_dir, "evidence_mc_raw.csv"))
message("Wrote ", file.path(out_dir, "evidence_mc_summary.csv"))
