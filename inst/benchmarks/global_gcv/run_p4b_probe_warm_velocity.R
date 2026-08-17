# P4B velocity: cold vs probe_warm GDF (same base fit, same probes)
#
#   Rscript inst/benchmarks/global_gcv/run_p4b_probe_warm_velocity.R
#
# Compares wall-clock of tt_global_gcv with gdf_init = cold | probe_warm.
# Does NOT enable θ→θ continuation.

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))

.rel <- function(a, b) abs(a - b) / (1 + abs(b))

.make <- function(n, d, seed) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  y <- rowSums(sin(2 * pi * X)) + rnorm(n, 0, 0.2)
  list(y = y, X = X)
}

.bench <- function(n, d, rank, M, max_sweeps, reps, seed, be = "Rcpp_fixed") {
  dat <- .make(n, d, seed)
  probes <- TTPsplines:::.tt_lab_rademacher_probes(n, M, probe_seed = seed + 11L)
  ctrl <- tt_control(max_sweeps = max_sweeps, tol = 1e-8, compute_edf = FALSE, seed = seed)
  init <- tt_initialize(d = d, rank = rank, k = 6L, seed = seed)
  lam <- rep(1.0, d)

  run_one <- function(gi) {
    TTPsplines:::tt_global_gcv(
      lambda = lam, y = dat$y, X = dat$X, rank = rank, probes = probes,
      M = M, control = ctrl, init = TTPsplines:::.tt_clone_cores(init),
      fit_backend = be, k = 6L, on_nonconverged = "na",
      gdf_init = gi, fidelity = "bench"
    )
  }

  invisible(run_one("probe_warm"))
  invisible(run_one("cold"))

  t_w <- t_c <- numeric(reps)
  out_w <- out_c <- NULL
  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    out_w <- run_one("probe_warm")
    t_w[i] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    out_c <- run_one("cold")
    t_c[i] <- proc.time()[["elapsed"]] - t0
  }

  data.frame(
    n = n, d = d, r = rank, M = M, max_sweeps = max_sweeps, backend = be,
    med_probe_warm_s = median(t_w), med_cold_s = median(t_c),
    speedup_warm_vs_cold = median(t_c) / median(t_w),
    gdf_rel = .rel(out_w$gdf, out_c$gdf),
    gcv_rel = .rel(out_w$global_gcv, out_c$global_gcv),
    base_ok_warm = isTRUE(out_w$gdf_result$base_cores_unchanged),
    stringsAsFactors = FALSE
  )
}

cat("P4B probe_warm vs cold velocity\n")
cat(sprintf("%s | %s %s\n\n", R.version.string, Sys.info()[["sysname"]], Sys.info()[["machine"]]))

rows <- list(
  .bench(500, 2, 2, M = 5, max_sweeps = 15, reps = 4, seed = 1),
  .bench(500, 3, 2, M = 5, max_sweeps = 15, reps = 3, seed = 2),
  .bench(1000, 3, 2, M = 8, max_sweeps = 12, reps = 3, seed = 3)
)
tab <- do.call(rbind, rows)
options(width = 140)
print(tab, row.names = FALSE, digits = 3)

out_dir <- "inst/benchmarks/global_gcv/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(out_dir, "p4b_probe_warm_velocity.csv"), row.names = FALSE)
cat("\nWrote", file.path(out_dir, "p4b_probe_warm_velocity.csv"), "\n")
