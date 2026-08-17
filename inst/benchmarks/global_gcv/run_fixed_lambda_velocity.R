#!/usr/bin/env Rscript
# Velocity: fixed-λ ALS R vs Rcpp (P1 core / P2 sweep / P3 fitter)
#
#   Rscript inst/benchmarks/global_gcv/run_fixed_lambda_velocity.R
#   TT_GGCV_VEL_MODE=QUICK Rscript ...   # smaller grid
#   TT_GGCV_VEL_MODE=FULL  Rscript ...   # default grid below
#
# Does not write large CSVs under deep paths (CRAN path-length hygiene).
# Optional: TT_GGCV_VEL_OUT=path/to/dir to save a short CSV.

mode <- toupper(Sys.getenv("TT_GGCV_VEL_MODE", "FULL"))
out_dir <- Sys.getenv("TT_GGCV_VEL_OUT", unset = "")

suppressPackageStartupMessages({
  root <- Sys.getenv("TT_PKG_ROOT", unset = "")
  if (!nzchar(root)) {
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", args, value = TRUE)
    if (length(f)) {
      root <- normalizePath(file.path(dirname(sub("^--file=", "", f)), "../..", ".."))
    } else {
      root <- normalizePath(".")
    }
  }
  if (file.exists(file.path(root, "DESCRIPTION"))) {
    pkgload::load_all(root, quiet = TRUE)
  } else {
    library(TTPsplines)
  }
})

.bench_fit <- function(n, d, rank, k_basis = 8L, max_sweeps = 20L,
                       tol = 1e-8, reps = 5L, seed = 1L) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  knots <- lapply(seq_len(d), function(j) {
    make_knots(X[, j], k = k_basis, degree = 3L)
  })
  basis <- eval_marginal_bases(X, knots, degree = 3L)
  ranks <- tt_rank(rank, d = d)
  cores <- initialize_tt_cores(ncol(basis[[1]]), ranks, seed = seed, sd = 0.2)
  y <- as.numeric(tt_contraction(cores, basis) + rnorm(n, sd = 0.2))
  lam <- rep(1.0, d)

  invisible(tt_als_fit_fixed_global(y, cores, basis, lam, max_sweeps = 2L, backend = "R"))
  invisible(tt_als_fit_fixed_global(y, cores, basis, lam, max_sweeps = 2L, backend = "Rcpp"))

  t_r <- t_c <- obj_r <- obj_c <- numeric(reps)
  sw_r <- sw_c <- integer(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    fr <- tt_als_fit_fixed_global(
      y, cores, basis, lam, max_sweeps = max_sweeps, tol = tol, backend = "R"
    )
    t_r[i] <- proc.time()[["elapsed"]] - t0
    sw_r[i] <- fr$n_sweeps
    obj_r[i] <- fr$objective

    t0 <- proc.time()[["elapsed"]]
    fc <- tt_als_fit_fixed_global(
      y, cores, basis, lam, max_sweeps = max_sweeps, tol = tol, backend = "Rcpp"
    )
    t_c[i] <- proc.time()[["elapsed"]] - t0
    sw_c[i] <- fc$n_sweeps
    obj_c[i] <- fc$objective
  }

  data.frame(
    component = "fit_P3",
    n = n, d = d, r = rank, k = k_basis, max_sweeps = max_sweeps, reps = reps,
    med_R_s = median(t_r), med_Cpp_s = median(t_c),
    speedup = median(t_r) / median(t_c),
    sweeps_R = median(sw_r), sweeps_Cpp = median(sw_c),
    rel_obj = abs(median(obj_c) - median(obj_r)) / (1 + abs(median(obj_r))),
    stringsAsFactors = FALSE
  )
}

.bench_micro <- function(n = 5000L, d = 5L, rank = 2L, k_basis = 8L,
                         reps = 20L, seed = 99L) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  knots <- lapply(seq_len(d), function(j) {
    make_knots(X[, j], k = k_basis, degree = 3L)
  })
  basis <- eval_marginal_bases(X, knots, degree = 3L)
  cores <- initialize_tt_cores(ncol(basis[[1]]), tt_rank(rank, d), seed = seed, sd = 0.2)
  y <- as.numeric(tt_contraction(cores, basis) + rnorm(n, sd = 0.2))
  lam <- rep(1, d)
  intercept <- mean(y)

  tr <- tc <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    invisible(tt_als_core_update_global(
      y, cores, intercept, basis, 3L, lam, backend = "R"
    ))
    tr[i] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    invisible(tt_als_core_update_global(
      y, cores, intercept, basis, 3L, lam, backend = "Rcpp"
    ))
    tc[i] <- proc.time()[["elapsed"]] - t0
  }
  core_row <- data.frame(
    component = "core_P1", n = n, d = d, r = rank, k = k_basis,
    max_sweeps = NA_integer_, reps = reps,
    med_R_s = median(tr), med_Cpp_s = median(tc),
    speedup = median(tr) / median(tc),
    sweeps_R = NA_real_, sweeps_Cpp = NA_real_, rel_obj = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    invisible(tt_als_sweep_global(
      y, cores, intercept, basis, lam, backend = "R"
    ))
    tr[i] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    invisible(tt_als_sweep_global(
      y, cores, intercept, basis, lam, backend = "Rcpp"
    ))
    tc[i] <- proc.time()[["elapsed"]] - t0
  }
  sweep_row <- data.frame(
    component = "sweep_P2", n = n, d = d, r = rank, k = k_basis,
    max_sweeps = NA_integer_, reps = reps,
    med_R_s = median(tr), med_Cpp_s = median(tc),
    speedup = median(tr) / median(tc),
    sweeps_R = NA_real_, sweeps_Cpp = NA_real_, rel_obj = NA_real_,
    stringsAsFactors = FALSE
  )
  rbind(core_row, sweep_row)
}

grid <- if (identical(mode, "QUICK")) {
  list(
    list(n = 500L, d = 3L, rank = 2L, reps = 5L),
    list(n = 5000L, d = 3L, rank = 2L, reps = 3L)
  )
} else {
  list(
    list(n = 500L, d = 3L, rank = 2L, reps = 7L),
    list(n = 500L, d = 5L, rank = 2L, reps = 5L),
    list(n = 5000L, d = 3L, rank = 2L, reps = 5L),
    list(n = 5000L, d = 5L, rank = 2L, reps = 4L),
    list(n = 5000L, d = 3L, rank = 5L, reps = 4L),
    list(n = 50000L, d = 3L, rank = 2L, reps = 3L),
    list(n = 5000L, d = 10L, rank = 2L, reps = 3L),
    list(n = 5000L, d = 5L, rank = 5L, reps = 3L)
  )
}

cat("TTPsplines fixed-lambda ALS velocity (R vs Rcpp)\n")
cat(sprintf("mode=%s | %s | %s %s\n\n",
            mode, R.version.string, Sys.info()[["sysname"]], Sys.info()[["machine"]]))

fit_rows <- lapply(grid, function(g) {
  cat(sprintf("... fit n=%d d=%d r=%d\n", g$n, g$d, g$rank))
  .bench_fit(n = g$n, d = g$d, rank = g$rank, reps = g$reps)
})
tab <- do.call(rbind, fit_rows)

if (!identical(mode, "QUICK")) {
  cat("... microbench P1/P2\n")
  tab <- rbind(tab, .bench_micro())
}

options(width = 140)
print(tab, row.names = FALSE, digits = 3)

if (nzchar(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_csv <- file.path(out_dir, "fixed_lambda_velocity.csv")
  utils::write.csv(tab, out_csv, row.names = FALSE)
  cat(sprintf("\nWrote %s\n", out_csv))
}

invisible(tab)
