## validate_gdf.R — accuracy of the array gGCV GDF away from saturating rank
##
## tt_ggcv_array() estimates GDF = tr(d mu_hat / d y) of the fixed-lambda array
## fit by Monte Carlo (Rademacher, non-negative perturbations). The package
## tests check it only at saturating rank in d = 2, where TT = dense GLAM and
## the exact EDF is available. Here: low TT rank, d = 3 and d = 5, Gaussian
## and zero-heavy Poisson. References per case:
##   gdf_unit_B     exact finite-difference trace of the SAME deterministic map
##                  (one non-negative unit perturbation per cell; budget B)
##   gdf_unit_3B    same with budget 3B: is the map converged at B?
##   edf_T          package T-EDF (left-orthogonal linearized EDF, compute_edf)
##                  at the budget-3B fit: an independent analytic comparator
##   gdf_mc_cold    Hutchinson, M probes, cold common init (default)
##   gdf_mc_warm    Hutchinson, M probes, warm probes with budget B / 3
##
## Run from the package root:  Rscript inst/benchmarks/ggcv_array/validate_gdf.R
## Env: GGV_CORES (6), GGV_M (60). Output: validate_gdf_results.csv (this folder).
## ---------------------------------------------------------------------------
pkg <- normalizePath(".")
suppressMessages(pkgload::load_all(pkg, quiet = TRUE))
out_dir <- file.path(pkg, "inst", "benchmarks", "ggcv_array")
n_cores <- as.integer(Sys.getenv("GGV_CORES", "6"))
M <- as.integer(Sys.getenv("GGV_M", "60"))

grid_axes <- function(ng) {
  ax <- lapply(ng, function(m) seq(0, 1, length.out = m))
  names(ax) <- paste0("x", seq_along(ng))
  ax
}
truth <- function(g) {
  d <- ncol(g)
  e <- sin(2 * pi * g[, 1]) + 0.8 * cos(2 * pi * g[, 2]) + 0.6 * g[, 1] * g[, 3]
  if (d >= 5) e <- e + 0.6 * g[, 3] * g[, 4] - 0.5 * g[, 5] + 0.7 * sin(pi * g[, 1]) * g[, 5]
  e
}
cases <- list(
  list(id = "gauss_d3_r2",    family = "gaussian", ng = c(8L, 7L, 6L), K = 5L, rank = 2L,      lambda = c(0.3, 1, 3)),
  list(id = "gauss_d3_rsat",  family = "gaussian", ng = c(8L, 7L, 6L), K = 5L, rank = c(5L, 5L), lambda = c(0.3, 1, 3)),
  list(id = "pois_d3_r2",     family = "poisson",  ng = c(8L, 7L, 6L), K = 5L, rank = 2L,      lambda = c(0.3, 1, 3)),
  list(id = "gauss_d5_r2",    family = "gaussian", ng = c(5L, 5L, 4L, 4L, 4L), K = 4L, rank = 2L, lambda = 1),
  list(id = "pois_d5_r2",     family = "poisson",  ng = c(5L, 5L, 4L, 4L, 4L), K = 4L, rank = 2L, lambda = 1),
  list(id = "pois_d5_r3",     family = "poisson",  ng = c(5L, 5L, 4L, 4L, 4L), K = 4L, rank = 3L, lambda = 1)
)

rows <- list()
for (cs in cases) {
  axes <- grid_axes(cs$ng)
  g <- as.matrix(expand.grid(axes))
  set.seed(20261001)
  Y <- if (identical(cs$family, "poisson")) {
    array(rpois(nrow(g), exp(log(0.3) + truth(g))), cs$ng)
  } else {
    array(truth(g) + rnorm(nrow(g), sd = 0.4), cs$ng)
  }
  fam <- if (identical(cs$family, "poisson")) poisson() else gaussian()
  B <- 10L
  ctl <- function(b) tt_control(max_sweeps = b, pirls_maxit = b, seed = 1L, compute_edf = FALSE)
  gdf <- function(b, probes, init = "cold", m = M) {
    t0 <- proc.time()[["elapsed"]]
    # budget = "fixed": this script studies the fixed budgets B and 3B
    z <- tt_gdf_array(Y, lambda = cs$lambda, axes = axes, family = fam, rank = cs$rank,
                      k = cs$K, M = m, probes = probes, probe_init = init,
                      n_cores = n_cores, budget = "fixed", control = ctl(b))
    z$time_s <- proc.time()[["elapsed"]] - t0
    z
  }
  u1 <- gdf(B, "unit")
  u3 <- gdf(3L * B, "unit")
  hc <- gdf(B, "rademacher", "cold")
  hw <- gdf(B, "rademacher", "warm")
  ce <- ctl(3L * B); ce$compute_edf <- TRUE; ce$tol <- 0
  fe <- ttps(Y, axes = axes, array = TRUE, family = fam, rank = cs$rank, k = cs$K,
             lambda = cs$lambda, optimizer = "ALS", control = ce)
  rows[[length(rows) + 1L]] <- data.frame(
    case = cs$id, family = cs$family, d = length(cs$ng), n_cells = prod(cs$ng),
    rank = paste(fe$rank, collapse = "-"), npar = fe$npar_tt, K = cs$K,
    frac_zero = mean(Y == 0), budget_B = B,
    gdf_unit_B = u1$gdf, gdf_unit_3B = u3$gdf, edf_T = as.numeric(fe$edf)[1],
    gdf_mc_cold = hc$gdf, gdf_mc_cold_se = hc$gdf_se,
    gdf_mc_warm = hw$gdf, gdf_mc_warm_se = hw$gdf_se, M = M,
    z_cold = (hc$gdf - u1$gdf) / hc$gdf_se,
    rel_budget = (u1$gdf - u3$gdf) / u3$gdf,
    rel_warm = (hw$gdf - u1$gdf) / u1$gdf,
    time_unit_B_s = u1$time_s, time_mc_cold_s = hc$time_s
  )
  print(rows[[length(rows)]])
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(out_dir, "validate_gdf_results.csv"), row.names = FALSE)
print(res[, c("case", "gdf_unit_B", "gdf_unit_3B", "edf_T", "gdf_mc_cold",
              "gdf_mc_cold_se", "z_cold", "gdf_mc_warm", "rel_warm")], digits = 4)
