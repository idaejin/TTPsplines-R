# d=3 strong_aniso: directed expansion of rough margin θ2 only.
#
# Question (locked):
#   Does the expanded margin reveal an interior minimum, or confirm an
#   asymptotic boundary solution with λ2 → 0?
#
# Adaptive box rule (locked):
#   1. min >0.5 from bound → interior
#   2. min <0.5 from lower θ2 → expand once to -10
#   3. if still descending at edge → boundary/asymptotic; stop
#   4. points within 1% of best = same valley
#
# Usage (package root):
#   Rscript inst/benchmarks/global_gcv/run_d3_aniso_expand.R

root <- Sys.getenv("TT_PKG_ROOT", unset = "")
if (!nzchar(root)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
    root <- normalizePath(file.path(dirname(script), "..", "..", ".."))
  } else {
    root <- getwd()
  }
}
setwd(root)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, export_all = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

out_dir <- file.path(root, "inst", "benchmarks", "global_gcv",
                     "results", "global_opt_d3_aniso_expand")
fig_dir <- file.path(root, "inst", "benchmarks", "global_gcv",
                     "figures", "global_opt_d3_aniso_expand")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Fixed design / controls (same as d3 GATE) --------------------------------
seed_data <- 1L
seed_opt <- 1L
n <- 250L
n_test <- 200L
k <- 6L
rank <- 6L
M_search <- 12L
M_final <- 30L
n_global <- 256L
n_refine <- 5L
n_diverse <- 10L
core_starts_final <- 3L
valley_tol <- 0.01
edge_tol <- 0.5

# Prior solution from d3 GATE strong_aniso (anchor)
prior_theta <- matrix(c(3.203125, -3.984375, -1.328125), nrow = 1L)

des <- TTPsplines:::.tt_lab_phase1_make_design_d3(
  scenario = "strong_aniso", n = n, n_test = n_test, k = k, seed = seed_data
)
ctrl <- tt_control(
  max_sweeps = 40L, tol = 1e-8, compute_edf = FALSE, seed = seed_data
)
common_init <- TTPsplines:::.tt_with_preserved_seed({
  tt_initialize(d = 3L, rank = rank, k = k, seed = seed_data)
})
probes_search <- TTPsplines:::.tt_lab_rademacher_probes(
  des$n, M_search, probe_seed = seed_data + 11L
)
probes_final <- TTPsplines:::.tt_lab_rademacher_probes(
  des$n, M_final, probe_seed = seed_data + 11L
)

.fair <- function(theta, probes, M) {
  TTPsplines:::.tt_lab_eval_ggcv(
    lambda = 10^as.numeric(theta), design = des, rank = rank,
    probes = probes, M = M, init_cores = common_init,
    init_policy = "cold_common", mc_bank_id = "expand",
    cache = NULL, control = ctrl
  )
}

.nonuniform_grid <- function(t2_lo) {
  # θ1, θ3: coarse on [-5,5]; θ2: dense on [t2_lo,-4], coarser above
  t1 <- seq(-5, 5, length.out = 5L)
  t3 <- t1
  t2_dense <- seq(t2_lo, -4, by = 0.5)
  t2_coarse <- seq(-3, 5, by = 2)
  t2 <- sort(unique(c(t2_dense, t2_coarse)))
  expand.grid(t1 = t1, t2 = t2, t3 = t3, KEEP.OUT.ATTRS = FALSE)
}

.eval_grid <- function(coords, probes, M, cache) {
  rows <- vector("list", nrow(coords))
  for (i in seq_len(nrow(coords))) {
    th <- c(coords$t1[i], coords$t2[i], coords$t3[i])
    ev <- TTPsplines:::.tt_lab_eval_ggcv(
      lambda = 10^th, design = des, rank = rank, probes = probes,
      M = M, init_cores = common_init, init_policy = "cold_common",
      mc_bank_id = "bank0", cache = cache, control = ctrl
    )
    rows[[i]] <- data.frame(
      t1 = th[1], t2 = th[2], t3 = th[3],
      global_gcv = ev$global_gcv, gdf = ev$gdf, rss = ev$rss,
      valid = ev$valid, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.classify_min <- function(theta, lower, upper, edge_tol = 0.5) {
  dist_lo <- abs(theta - lower)
  dist_hi <- abs(theta - upper)
  near <- (dist_lo < edge_tol) | (dist_hi < edge_tol)
  list(
    interior = !any(near),
    near_lower = dist_lo < edge_tol,
    near_upper = dist_hi < edge_tol,
    min_edge_dist = min(c(dist_lo, dist_hi))
  )
}

.run_stage <- function(t2_lo, label) {
  lower <- c(-5, t2_lo, -5)
  upper <- c(5, 5, 5)
  message(sprintf("=== stage %s: θ2 ∈ [%.1f, 5] ===", label, t2_lo))

  # Non-uniform reference grid
  coords <- .nonuniform_grid(t2_lo)
  cache <- new.env(parent = emptyenv())
  t0 <- proc.time()[["elapsed"]]
  grid <- .eval_grid(coords, probes_search, M_search, cache)
  grid_elapsed <- proc.time()[["elapsed"]] - t0
  i_g <- which.min(grid$global_gcv)
  grid_th <- c(grid$t1[i_g], grid$t2[i_g], grid$t3[i_g])

  opt <- TTPsplines:::tt_global_lambda_optimize(
    y = des$y, X = des$X, rank = rank, k = k,
    degree = des$degree, penalty_order = des$penalty_order,
    theta_lower = lower, theta_upper = upper,
    n_global = n_global, n_refine = n_refine, n_diverse = n_diverse,
    M_search = M_search, M_final = M_final,
    core_starts_final = core_starts_final,
    seed = seed_opt, control = ctrl,
    include_cgcv_anchor = TRUE,
    extra_theta = prior_theta,
    verbose = TRUE
  )

  ev_g <- .fair(grid_th, probes_final, M_final)
  ev_o <- .fair(opt$theta, probes_final, M_final)
  q_g <- ev_g$global_gcv
  q_o <- ev_o$global_gcv
  in_valley <- is.finite(q_o) && is.finite(q_g) &&
    q_o <= (1 + valley_tol) * q_g
  # also: grid best within 1% of opt (same valley either way)
  same_valley <- is.finite(q_o) && is.finite(q_g) &&
    (q_o <= (1 + valley_tol) * q_g || q_g <= (1 + valley_tol) * q_o)

  cls <- .classify_min(opt$theta, lower, upper, edge_tol)
  aniso_ok <- (opt$theta[2] < opt$theta[1] - 0.25) &&
    (opt$theta[2] < opt$theta[3] - 0.25)

  list(
    label = label,
    lower = lower, upper = upper,
    grid = grid, grid_theta = grid_th,
    opt = opt,
    q_grid = q_g, q_opt = q_o,
    delta_rel = (q_o - q_g) / q_g,
    in_valley = in_valley,
    same_valley = same_valley,
    class = cls,
    aniso_ok = aniso_ok,
    grid_elapsed = grid_elapsed,
    n_grid = nrow(grid),
    n_explore = opt$cost$n_explore_points %||% NA_integer_,
    n_als_opt = opt$cost$n_als_fits_total_approx %||% NA_integer_,
    n_als_grid = as.integer(nrow(grid) * (1L + M_search))
  )
}

# --- Stage A: θ2 ∈ [-8, 5] ----------------------------------------------------
stageA <- .run_stage(-8, "A_t2lo_-8")
utils::write.csv(stageA$grid, file.path(out_dir, "stageA_grid.csv"),
                 row.names = FALSE)
utils::write.csv(stageA$opt$sobol_results,
                 file.path(out_dir, "stageA_explore.csv"), row.names = FALSE)
utils::write.csv(stageA$opt$final_candidates,
                 file.path(out_dir, "stageA_final.csv"), row.names = FALSE)

need_expand <- isTRUE(stageA$class$near_lower[2]) ||
  (!isTRUE(stageA$class$interior) && stageA$opt$theta[2] < -8 + edge_tol)

stageB <- NULL
if (need_expand) {
  message("θ2 near lower edge → expand once to -10")
  stageB <- .run_stage(-10, "B_t2lo_-10")
  utils::write.csv(stageB$grid, file.path(out_dir, "stageB_grid.csv"),
                   row.names = FALSE)
  utils::write.csv(stageB$opt$final_candidates,
                   file.path(out_dir, "stageB_final.csv"), row.names = FALSE)
}

final_stage <- if (!is.null(stageB)) stageB else stageA

# --- Profile: θ2 ↦ min_{θ1,θ3} Q  (from densest available grid) --------------
prof_src <- final_stage$grid
# Prefer M_final re-score along a θ2 path at fixed (θ1*,θ3*) from opt / grid best
t2_path <- seq(final_stage$lower[2], 5, by = 0.5)
# Use θ1,θ3 from optimizer best for conditional profile
th13 <- final_stage$opt$theta[c(1, 3)]
prof_rows <- vector("list", length(t2_path))
for (i in seq_along(t2_path)) {
  th <- c(th13[1], t2_path[i], th13[2])
  # clamp to box
  th <- pmin(pmax(th, final_stage$lower), final_stage$upper)
  ev <- .fair(th, probes_final, M_final)
  mu_te <- tryCatch(
    as.numeric(predict(ev$fit, newdata = des$X_test, type = "response")),
    error = function(e) rep(NA_real_, des$n_test)
  )
  ise <- if (all(is.finite(mu_te))) mean((mu_te - des$f_test)^2) else NA_real_
  # Also min over θ1,θ3 on the searched grid at this θ2 (nearest)
  near <- which(abs(prof_src$t2 - t2_path[i]) < 0.26)
  q_slice <- if (length(near)) min(prof_src$global_gcv[near], na.rm = TRUE) else NA_real_
  prof_rows[[i]] <- data.frame(
    t2 = t2_path[i],
    q_at_opt13 = ev$global_gcv,
    gdf = ev$gdf,
    rss = ev$rss,
    ise_test = ise,
    q_min_slice_Msearch = q_slice,
    stringsAsFactors = FALSE
  )
}
profile <- do.call(rbind, prof_rows)
utils::write.csv(profile, file.path(out_dir, "theta2_profile.csv"),
                 row.names = FALSE)

# Focused checkpoints θ2 = -5,-6,-7,-8 (and -9,-10 if expanded)
check_t2 <- c(-5, -6, -7, -8)
if (!is.null(stageB)) check_t2 <- c(check_t2, -9, -10)
check_rows <- vector("list", length(check_t2))
for (i in seq_along(check_t2)) {
  th <- c(th13[1], check_t2[i], th13[2])
  th <- pmin(pmax(th, final_stage$lower), final_stage$upper)
  ev <- .fair(th, probes_final, M_final)
  mu_te <- tryCatch(
    as.numeric(predict(ev$fit, newdata = des$X_test, type = "response")),
    error = function(e) rep(NA_real_, des$n_test)
  )
  ise <- if (all(is.finite(mu_te))) mean((mu_te - des$f_test)^2) else NA_real_
  yhat <- if (!is.null(ev$fit)) as.numeric(fitted(ev$fit)) else rep(NA_real_, des$n)
  check_rows[[i]] <- data.frame(
    t2 = check_t2[i],
    global_gcv = ev$global_gcv,
    gdf = ev$gdf,
    rss = ev$rss,
    ise_test = ise,
    yhat_l2 = sqrt(mean(yhat^2)),
    stringsAsFactors = FALSE
  )
}
checks <- do.call(rbind, check_rows)
utils::write.csv(checks, file.path(out_dir, "t2_checkpoints.csv"),
                 row.names = FALSE)

# --- Profile plot (base R, ASCII labels) --------------------------------------
png(file.path(fig_dir, "theta2_profile.png"), width = 900, height = 500, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(profile$t2, profile$q_at_opt13, type = "b", pch = 16,
     xlab = "theta2 = log10(lambda2)", ylab = "TT-gGCV Q",
     main = "Profile Q at fixed (theta1, theta3)")
abline(v = final_stage$opt$theta[2], lty = 2, col = "gray40")
plot(checks$t2, checks$ise_test, type = "b", pch = 16, col = "darkred",
     xlab = "theta2 = log10(lambda2)", ylab = "test ISE",
     main = "ISE vs theta2")
par(op)
dev.off()

# --- Verdict ------------------------------------------------------------------
opt_th <- final_stage$opt$theta
cls <- final_stage$class
# Monotone descent toward lower edge?
prof_lo <- profile[profile$t2 <= -4, , drop = FALSE]
# Does Q keep decreasing as t2 decreases near the edge?
edge_band <- profile[profile$t2 <= final_stage$lower[2] + 1.5, , drop = FALSE]
descending <- if (nrow(edge_band) >= 2L) {
  # correlation of t2 with Q: positive => Q decreases as t2 decreases
  stats::cor(edge_band$t2, edge_band$q_at_opt13, use = "complete.obs") > 0.3
} else {
  FALSE
}

# Stabilization of fit: relative ISE change from -5 to most negative check
ise_vals <- checks$ise_test[is.finite(checks$ise_test)]
ise_stable <- if (length(ise_vals) >= 2L) {
  (max(ise_vals) - min(ise_vals)) / max(mean(ise_vals), 1e-12) < 0.05
} else {
  NA
}

if (isTRUE(cls$interior) && !descending) {
  verdict <- "INTERIOR_MIN"
} else if (isTRUE(cls$near_lower[2]) && (descending || isTRUE(ise_stable))) {
  verdict <- "BOUNDARY_ASYMPTOTIC_LAMBDA2"
} else if (isTRUE(cls$near_lower[2])) {
  verdict <- "BOUNDARY_UNRESOLVED"
} else {
  verdict <- "INTERIOR_CANDIDATE"
}

summary_row <- data.frame(
  verdict = verdict,
  stage = final_stage$label,
  t2_lo = final_stage$lower[2],
  opt_t1 = opt_th[1], opt_t2 = opt_th[2], opt_t3 = opt_th[3],
  grid_t1 = final_stage$grid_theta[1],
  grid_t2 = final_stage$grid_theta[2],
  grid_t3 = final_stage$grid_theta[3],
  q_opt = final_stage$q_opt,
  q_grid = final_stage$q_grid,
  delta_rel = final_stage$delta_rel,
  in_valley = final_stage$in_valley,
  same_valley = final_stage$same_valley,
  aniso_ok = final_stage$aniso_ok,
  interior = cls$interior,
  near_t2_lower = cls$near_lower[2],
  min_edge_dist = cls$min_edge_dist,
  descending_to_edge = descending,
  ise_stable_along_t2 = ise_stable,
  winner_source = final_stage$opt$winner_source,
  n_explore = final_stage$n_explore,
  n_grid = final_stage$n_grid,
  n_als_opt = final_stage$n_als_opt,
  n_als_grid = final_stage$n_als_grid,
  expanded_to_m10 = !is.null(stageB),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_row, file.path(out_dir, "summary.csv"), row.names = FALSE)
writeLines(c(verdict, paste(capture.output(print(summary_row)), collapse = "\n")),
           file.path(out_dir, "VERDICT.txt"))

# Stage A snapshot
stageA_row <- data.frame(
  stage = "A", t2_lo = -8,
  opt_t2 = stageA$opt$theta[2],
  grid_t2 = stageA$grid_theta[2],
  near_t2_lower = stageA$class$near_lower[2],
  interior = stageA$class$interior,
  q_opt = stageA$q_opt, q_grid = stageA$q_grid,
  delta_rel = stageA$delta_rel,
  in_valley = stageA$in_valley,
  stringsAsFactors = FALSE
)
utils::write.csv(stageA_row, file.path(out_dir, "stageA_summary.csv"),
                 row.names = FALSE)

cat("\n==== D3 ANISO EXPAND VERDICT:", verdict, "====\n")
print(summary_row)
message("Wrote ", out_dir)
message("Figure ", file.path(fig_dir, "theta2_profile.png"))
