# P4B end-to-end lab velocity: adaptive fidelity + probe_warm
#
#   Rscript inst/benchmarks/global_gcv/run_p4b_lab_velocity.R
#
# Compares:
#   A) adaptive_fidelity=FALSE, flat high sweeps (legacy-like)
#   B) adaptive_fidelity=TRUE  (sobol 12 / refine 25 / final 50)
# Both use fit_backend=Rcpp_fixed and gdf_init=probe_warm.

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))

.make <- function(n, d, seed, aniso = FALSE) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  if (isTRUE(aniso) && d >= 3L) {
    f <- sin(2 * pi * X[, 1]) + 0.05 * X[, 2] + sin(pi * X[, 3])
  } else {
    f <- rowSums(sin(2 * pi * X))
  }
  list(y = f + rnorm(n, 0, 0.2), X = X)
}

.run_opt <- function(dat, adaptive, label) {
  t0 <- proc.time()[["elapsed"]]
  out <- TTPsplines:::tt_global_lambda_optimize(
    y = dat$y, X = dat$X, rank = 2L,
    theta_lower = -2, theta_upper = 2,
    n_global = 16L, n_refine = 2L, n_diverse = 4L,
    M_search = 4L, M_final = 8L,
    core_starts_final = 1L,
    include_cgcv_anchor = FALSE,
    k = 5L,
    control = tt_control(max_sweeps = 40L, tol = 1e-8,
                         compute_edf = FALSE, seed = 11L),
    fit_backend = "Rcpp_fixed",
    adaptive_fidelity = adaptive,
    gdf_init = "probe_warm",
    verbose = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  data.frame(
    label = label,
    adaptive = adaptive,
    elapsed_s = elapsed,
    n_als_search = out$cost$n_als_fits_search,
    n_als_final = out$cost$n_als_fits_final_approx,
    n_als_total = out$cost$n_als_fits_total_approx,
    sobol_sweeps = out$diagnostics$fidelity$sobol_sweeps,
    refine_sweeps = out$diagnostics$fidelity$refine_sweeps,
    final_sweeps = out$diagnostics$fidelity$final_sweeps,
    gcv = out$gcv,
    winner_source = out$winner_source,
    ranking_stable = isTRUE(out$diagnostics$ranking_stable_alt_bank),
    stringsAsFactors = FALSE
  )
}

cat("P4B lab velocity (adaptive fidelity)\n")
cat(sprintf("%s | %s %s\n\n", R.version.string, Sys.info()[["sysname"]], Sys.info()[["machine"]]))

dat2 <- .make(400, 2, seed = 21)
dat3 <- .make(400, 3, seed = 22, aniso = TRUE)

rows <- rbind(
  .run_opt(dat2, FALSE, "d2_flat"),
  .run_opt(dat2, TRUE, "d2_adaptive"),
  .run_opt(dat3, FALSE, "d3_flat"),
  .run_opt(dat3, TRUE, "d3_adaptive")
)
options(width = 160)
print(rows, row.names = FALSE, digits = 3)

out_dir <- "inst/benchmarks/global_gcv/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(rows, file.path(out_dir, "p4b_lab_velocity.csv"), row.names = FALSE)
cat("\nWrote", file.path(out_dir, "p4b_lab_velocity.csv"), "\n")
