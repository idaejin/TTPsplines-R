# Experimental lab: epsilon_rel x M x scheme sensitivity for MC-GDF.
# Nested probes: first M columns of a shared probe bank.
#
# Usage:
#   Rscript inst/benchmarks/global_gcv/run_epsilon_M_sensitivity.R

source(file.path("inst", "benchmarks", "helpers.R"))
.ttps_bench_ensure_pkg()

out_dir <- file.path("inst", "benchmarks", "global_gcv", "results")
fig_dir <- file.path("inst", "benchmarks", "global_gcv", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260816)
n <- 80L
d <- 2L
k <- 5L
rank <- 5L # sufficient for d=2, k=5
lambda <- c(0.5, 2)
X <- matrix(runif(n * d), n, d)
f_true <- sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2])
y <- f_true + rnorm(n, 0, 0.2)

ctrl <- tt_control(
  max_sweeps = 40L, tol = 1e-10, compute_edf = FALSE, seed = 1L, backend = "R"
)
fit <- ttps(
  y, X, family = gaussian(), rank = rank, k = k, lambda = lambda,
  optimizer = "ALS", backend = "R", control = ctrl
)
bs <- TTPsplines:::build_marginal_bases(X, k = k, degree = 3L)
full <- TTPsplines:::.tt_lab_full_tp_gaussian(y, bs$basis, lambda = lambda)
probe_bank <- TTPsplines:::.tt_lab_rademacher_probes(n, M = 20L, probe_seed = 11L)

eps_grid <- c(1e-2, 1e-3, 1e-4)
M_grid <- c(3L, 5L, 10L, 20L)
schemes <- c("forward", "central")

rows <- list()
idx <- 0L
for (scheme in schemes) {
  for (eps_rel in eps_grid) {
    for (M in M_grid) {
      idx <- idx + 1L
      t0 <- proc.time()[["elapsed"]]
      g <- TTPsplines:::tt_global_gdf_mc(
        fit, y = y, M = M,
        probes = probe_bank[, seq_len(M), drop = FALSE],
        scheme = scheme, epsilon_rel = eps_rel, warm_start = TRUE,
        on_nonconverged = "na", control = ctrl
      )
      elapsed <- proc.time()[["elapsed"]] - t0
      rows[[idx]] <- data.frame(
        scheme = scheme,
        epsilon_rel = eps_rel,
        epsilon = g$epsilon,
        M = M,
        gdf = g$gdf,
        gdf_mc_se = g$gdf_mc_se,
        M_ok = g$M_ok,
        n_refits = g$n_refits,
        n_fail = sum(!g$perturbation_diagnostics$ok, na.rm = TRUE),
        edf_full = full$edf,
        abs_err_vs_full = abs(g$gdf - full$edf),
        elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      )
      message(sprintf(
        "[%s] eps_rel=%.0e M=%d gdf=%.3f se=%.3f err=%.3f t=%.1fs",
        scheme, eps_rel, M, g$gdf, g$gdf_mc_se, abs(g$gdf - full$edf), elapsed
      ))
    }
  }
}
tab <- do.call(rbind, rows)
utils::write.csv(tab, file.path(out_dir, "epsilon_M_sensitivity.csv"), row.names = FALSE)

# Provisional exploration pick (not a public default): prefer central, eps_rel=1e-3, M=10
# if finite and low failure count among that slice.
slice <- tab[tab$scheme == "central" & tab$epsilon_rel == 1e-3 & tab$M == 10, ]
prov <- if (nrow(slice) && is.finite(slice$gdf[1])) {
  list(scheme = "central", epsilon_rel = 1e-3, M = 10L, note = "provisional lab only")
} else {
  list(scheme = "forward", epsilon_rel = 1e-3, M = 5L, note = "fallback provisional")
}
utils::write.csv(
  as.data.frame(prov, stringsAsFactors = FALSE),
  file.path(out_dir, "provisional_epsilon_M.csv"),
  row.names = FALSE
)

# Figures: GDF vs epsilon_rel (M fixed = 10), GDF vs M (eps fixed = 1e-3)
png(file.path(fig_dir, "gdf_vs_epsilon_rel.png"), width = 640, height = 420)
sub <- tab[tab$M == 10L, ]
plot(NA, xlim = range(log10(sub$epsilon_rel)),
     ylim = range(sub$gdf, full$edf, na.rm = TRUE),
     xlab = "log10(epsilon_rel)", ylab = "GDF estimate",
     main = "MC-GDF vs epsilon_rel (M=10, nested probes)")
abline(h = full$edf, lty = 2, col = "blue")
for (sch in schemes) {
  s <- sub[sub$scheme == sch, ]
  lines(log10(s$epsilon_rel), s$gdf, type = "b",
        pch = if (sch == "forward") 16 else 17,
        col = if (sch == "forward") "black" else "darkred")
}
legend("topright", c("forward", "central", "EDF_full"),
       pch = c(16, 17, NA), lty = c(1, 1, 2),
       col = c("black", "darkred", "blue"), bty = "n")
dev.off()

png(file.path(fig_dir, "gdf_vs_M.png"), width = 640, height = 420)
sub <- tab[tab$epsilon_rel == 1e-3, ]
plot(NA, xlim = range(sub$M), ylim = range(sub$gdf, full$edf, na.rm = TRUE),
     xlab = "M (nested probes)", ylab = "GDF estimate",
     main = "MC-GDF vs M (epsilon_rel=1e-3)")
abline(h = full$edf, lty = 2, col = "blue")
for (sch in schemes) {
  s <- sub[sub$scheme == sch, ]
  lines(s$M, s$gdf, type = "b",
        pch = if (sch == "forward") 16 else 17,
        col = if (sch == "forward") "black" else "darkred")
}
legend("topright", c("forward", "central", "EDF_full"),
       pch = c(16, 17, NA), lty = c(1, 1, 2),
       col = c("black", "darkred", "blue"), bty = "n")
dev.off()

message("Provisional lab settings: ", paste(names(prov), prov, sep = "=", collapse = ", "))
print(tab)
