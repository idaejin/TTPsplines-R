#!/usr/bin/env Rscript
# Minimal gate: joint vs profiled (Gaussian, d=2, fixed lambda).
# See GATE_PROFILED.md

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
.root <- if (length(.file)) {
  normalizePath(file.path(dirname(.file), "..", "..", ".."))
} else {
  getwd()
}
if (!exists("ttps", mode = "function")) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(.root, quiet = TRUE)
  } else {
    stop("Need pkgload::load_all or installed TTPsplines.")
  }
}

rmse <- function(a, b) sqrt(mean((a - b)^2))
pen_obj <- function(fit) {
  0.5 * sum((fit$y - fitted(fit))^2) +
    tt_global_penalty_value(
      fit$cores, fit$lambda,
      penalty_order = fit$penalty_order, cyclic = fit$cyclic
    )
}

fit_modes <- function(y, X, lambda, rank = 3L, k = 8L, seed = 1L) {
  ctrl <- tt_control(
    max_sweeps = 15L, compute_edf = TRUE, backend = "R", seed = seed,
    warn_lambda_boundary = FALSE
  )
  out <- list()
  for (ns in c("joint", "profiled")) {
    out[[ns]] <- ttps(
      y, X, rank = rank, k = k, lambda = lambda, penalty_order = 2L,
      null_space = ns, control = ctrl
    )
  }
  out
}

summarize <- function(fits, f_true, label) {
  rows <- lapply(names(fits), function(ns) {
    fit <- fits[[ns]]
    nsinfo <- fit$null_space_info
    data.frame(
      scenario = label,
      mode = ns,
      rmse_f = rmse(fitted(fit), f_true),
      obj = pen_obj(fit),
      edf = fit$edf %||% NA_real_,
      edf_tt = fit$edf_tt %||% NA_real_,
      p0 = nsinfo$design_rank %||% NA_real_,
      ortho = if (!is.null(nsinfo$ortho_X0_mu_tt_perp))
        max(abs(nsinfo$ortho_X0_mu_tt_perp)) else NA_real_,
      mu_perp_rms = if (!is.null(nsinfo$mu_tt_perp))
        sqrt(mean(nsinfo$mu_tt_perp^2)) else NA_real_,
      identity = nsinfo$identity_max_abs %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

set.seed(20260817)
n <- 200L
X <- cbind(runif(n), runif(n))
colnames(X) <- c("x1", "x2")
aff <- 1.2 + 1.5 * X[, 1] - 0.8 * X[, 2]
wiggle <- sin(4 * pi * X[, 1]) * sin(4 * pi * X[, 2])
trend_wiggle <- 0.5 * X[, 1] + wiggle

scenarios <- list(
  null_only = list(f = aff, lambda = 1e6, seed = 11L),
  null_plus_rough = list(f = aff + 0.5 * wiggle, lambda = 10, seed = 12L),
  rough_only = list(f = wiggle - mean(wiggle), lambda = 5, seed = 13L),
  correlated_rough = list(f = trend_wiggle, lambda = 10, seed = 14L)
)

tab <- list()
for (nm in names(scenarios)) {
  sc <- scenarios[[nm]]
  y <- sc$f + rnorm(n, sd = 0.1)
  fits <- fit_modes(y, X, lambda = sc$lambda, seed = sc$seed)
  tab[[nm]] <- summarize(fits, sc$f, nm)
  cat("\n====", nm, "lambda=", sc$lambda, "====\n")
  print(tab[[nm]][, c("mode", "rmse_f", "obj", "edf", "ortho", "mu_perp_rms")],
        digits = 4, row.names = FALSE)
}

res <- do.call(rbind, tab)
out_dir <- file.path(.root, "inst", "benchmarks", "null_space")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, "gate_profiled_results.csv")
utils::write.csv(res, out_csv, row.names = FALSE)

prof <- res[res$mode == "profiled", ]
checks <- c(
  ortho_ok = all(prof$ortho < 1e-5, na.rm = TRUE),
  identity_ok = all(prof$identity < 1e-8, na.rm = TRUE),
  null_mu_small = {
    r <- prof[prof$scenario == "null_only", ]
    isTRUE(r$mu_perp_rms < 1e-3)
  }
)
cat("\nGATE CHECKS:\n")
print(checks)
cat("Wrote", out_csv, "\n")
if (!all(unlist(checks))) {
  quit(status = 1L)
}
