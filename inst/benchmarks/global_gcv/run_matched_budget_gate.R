# Matched-budget gate: Sobol vs random vs LHS on joint TT-gGCV.
#
# Primary question (PROTOCOL_MATCHED_BUDGET.md):
#   Does Sobol locate the near-optimal Q_r region with fewer evaluations
#   than random search or LHS at the same initial N?
#
# Product lock: KEEP cGCV. This driver does not change defaults.
#
# Usage (package root):
#   TT_GGCV_MB_MODE=SMOKE Rscript inst/benchmarks/global_gcv/run_matched_budget_gate.R
#   TT_GGCV_MB_MODE=GATE  Rscript inst/benchmarks/global_gcv/run_matched_budget_gate.R
#
# Optional:
#   TT_GGCV_MB_METHODS=random,lhs,sobol_pure,sobol_refine
#   TT_GGCV_MB_SCENARIOS=smooth_smooth,strong_aniso
#   TT_GGCV_MB_N_MULT=16,32,64

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
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, export_all = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(TTPsplines))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

mode <- toupper(Sys.getenv("TT_GGCV_MB_MODE", "SMOKE"))
out_dir <- file.path(
  root, "inst", "benchmarks", "global_gcv", "results",
  paste0("matched_budget_", tolower(mode))
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

.cfg_for_mode <- function(mode) {
  if (identical(mode, "GATE")) {
    list(
      scenarios = c("smooth_smooth", "strong_aniso"),
      n = 200L, n_test = 200L, k = 8L, rank = 8L,
      box = c(-5, 5),
      grid_by = 1.0,
      n_mult = c(16L, 32L, 64L),
      seeds_design = c(1L, 11L, 21L, 31L, 41L),
      seed_data = 1L,
      n_refine = 5L, n_diverse = 10L, min_dist = 0.75,
      M_search = 15L, M_final = 40L,
      valley_tol = 0.01,
      methods = c("random", "lhs", "sobol_pure", "sobol_refine"),
      include_cgcv = TRUE,
      include_grid_ref = TRUE
    )
  } else {
    # SMOKE: cheap algorithmic check
    list(
      scenarios = c("smooth_smooth"),
      n = 80L, n_test = 60L, k = 5L, rank = 5L,
      box = c(-5, 5),
      grid_by = 1.0,  # denser ref; 1.5 left cGCV beating the grid
      n_mult = c(16L),
      seeds_design = c(1L, 11L),
      seed_data = 1L,
      n_refine = 3L, n_diverse = 4L, min_dist = 0.75,
      M_search = 5L, M_final = 10L,
      valley_tol = 0.01,
      methods = c("random", "lhs", "sobol_pure", "sobol_refine"),
      include_cgcv = TRUE,
      include_grid_ref = TRUE
    )
  }
}

.cfg_apply_env <- function(cfg) {
  meth <- Sys.getenv("TT_GGCV_MB_METHODS", "")
  if (nzchar(meth)) {
    cfg$methods <- strsplit(meth, ",", fixed = TRUE)[[1]]
    cfg$methods <- trimws(cfg$methods)
  }
  sc <- Sys.getenv("TT_GGCV_MB_SCENARIOS", "")
  if (nzchar(sc)) {
    cfg$scenarios <- strsplit(sc, ",", fixed = TRUE)[[1]]
    cfg$scenarios <- trimws(cfg$scenarios)
  }
  nm <- Sys.getenv("TT_GGCV_MB_N_MULT", "")
  if (nzchar(nm)) {
    cfg$n_mult <- as.integer(strsplit(nm, ",", fixed = TRUE)[[1]])
  }
  cfg
}

# ---------------------------------------------------------------------------
# Design generators (unit cube -> box). Design seed only; not MC/cores.
# ---------------------------------------------------------------------------

.unit_to_box <- function(U, lower, upper) {
  U <- as.matrix(U)
  d <- length(lower)
  stopifnot(ncol(U) == d)
  sweep(U, 2L, upper - lower, `*`) +
    matrix(lower, nrow = nrow(U), ncol = d, byrow = TRUE)
}

.design_random <- function(n, lower, upper, seed) {
  d <- length(lower)
  U <- TTPsplines:::.tt_with_preserved_seed({
    set.seed(as.integer(seed))
    matrix(stats::runif(n * d), n, d)
  })
  .unit_to_box(U, lower, upper)
}

.design_lhs <- function(n, lower, upper, seed) {
  d <- length(lower)
  n <- as.integer(n)
  U <- TTPsplines:::.tt_with_preserved_seed({
    set.seed(as.integer(seed))
    out <- matrix(NA_real_, n, d)
    for (j in seq_len(d)) {
      # centered LHS with random stratum permutation
      out[, j] <- (sample.int(n) - stats::runif(n)) / n
    }
    out
  })
  U <- pmin(pmax(U, .Machine$double.eps), 1 - .Machine$double.eps)
  .unit_to_box(U, lower, upper)
}

.design_sobol <- function(n, lower, upper, seed) {
  d <- length(lower)
  # Builtin Joe–Kuo Sobol + Cranley–Patterson rotation (seeded)
  U0 <- TTPsplines:::.tt_lab_sobol_unit(n, d, skip = 1L)
  shift <- TTPsplines:::.tt_with_preserved_seed({
    set.seed(as.integer(seed))
    stats::runif(d)
  })
  U <- sweep(U0, 2L, shift, `+`) %% 1
  U <- pmin(pmax(U, .Machine$double.eps), 1 - .Machine$double.eps)
  .unit_to_box(U, lower, upper)
}

.make_design <- function(method, n, lower, upper, seed) {
  base <- switch(
    method,
    random = .design_random(n, lower, upper, seed),
    lhs = .design_lhs(n, lower, upper, seed),
    sobol_pure = .design_sobol(n, lower, upper, seed),
    sobol_refine = .design_sobol(n, lower, upper, seed),
    random_refine = .design_random(n, lower, upper, seed),
    lhs_refine = .design_lhs(n, lower, upper, seed),
    stop("Unknown method: ", method, call. = FALSE)
  )
  list(
    theta = base,
    do_refine = grepl("_refine$", method),
    source = method
  )
}

# ---------------------------------------------------------------------------
# Evaluate design + optional refine; fair M_final re-score
# ---------------------------------------------------------------------------

.eval_points <- function(theta_mat, evaluator, M, source_lab) {
  rows <- vector("list", nrow(theta_mat))
  n_failed <- 0L
  for (i in seq_len(nrow(theta_mat))) {
    ev <- evaluator$eval(theta_mat[i, ], M = M, mc_bank_id = "bank0",
                         return_fit = FALSE)
    if (!isTRUE(ev$valid) || !is.finite(ev$global_gcv)) n_failed <- n_failed + 1L
    rows[[i]] <- data.frame(
      i = i,
      source = source_lab,
      matrix(theta_mat[i, ], nrow = 1L,
             dimnames = list(NULL, paste0("theta", seq_len(ncol(theta_mat))))),
      global_gcv = ev$global_gcv,
      gdf = ev$gdf,
      gdf_se = ev$gdf_mc_se %||% NA_real_,
      rss = ev$rss,
      valid = isTRUE(ev$valid),
      stringsAsFactors = FALSE
    )
  }
  list(df = do.call(rbind, rows), n_failed = n_failed)
}

.refine_elites <- function(explore_df, evaluator, lower, upper,
                           n_refine, n_diverse, min_dist, M_search) {
  d <- length(lower)
  th_cols <- paste0("theta", seq_len(d))
  th_mat <- as.matrix(explore_df[, th_cols, drop = FALSE])
  q_vec <- explore_df$global_gcv
  elite_idx <- TTPsplines:::.tt_lab_diverse_top(
    th_mat, q_vec, n_keep = max(as.integer(n_diverse), 1L),
    min_dist = min_dist
  )
  if (!length(elite_idx)) elite_idx <- which.min(q_vec)
  n_ref_use <- max(0L, as.integer(n_refine))
  starts <- elite_idx[seq_len(min(n_ref_use, length(elite_idx)))]
  if (!length(starts)) {
    return(list(
      refined = data.frame(),
      n_refine_starts = 0L,
      n_refine_evals = 0L
    ))
  }
  n_miss0 <- evaluator$cost$n_theta_miss
  rows <- vector("list", length(starts))
  for (j in seq_along(starts)) {
    i0 <- starts[j]
    th0 <- as.numeric(th_mat[i0, ])
    obj <- function(th) {
      ev <- evaluator$eval(th, M = M_search, mc_bank_id = "bank0",
                           return_fit = FALSE, use_cache = FALSE)
      q <- ev$global_gcv
      if (!is.finite(q)) 1e30 else q
    }
    opt <- tryCatch(
      stats::nlminb(
        start = th0, objective = obj, lower = lower, upper = upper,
        control = list(abs.tol = 1e-8, rel.tol = 1e-8, iter.max = 80L)
      ),
      error = function(e) list(
        par = th0, objective = obj(th0), convergence = 99L
      )
    )
    th_ref <- as.numeric(opt$par)
    ev_ref <- evaluator$eval(th_ref, M = M_search, mc_bank_id = "bank0",
                             return_fit = FALSE)
    rows[[j]] <- data.frame(
      start_i = i0,
      start_gcv = q_vec[i0],
      matrix(th_ref, nrow = 1L, dimnames = list(NULL, th_cols)),
      global_gcv = ev_ref$global_gcv,
      nlminb_conv = as.integer(opt$convergence %||% NA_integer_),
      improved = is.finite(ev_ref$global_gcv) && is.finite(q_vec[i0]) &&
        ev_ref$global_gcv < q_vec[i0] - 1e-12,
      stringsAsFactors = FALSE
    )
  }
  list(
    refined = do.call(rbind, rows),
    n_refine_starts = length(starts),
    n_refine_evals = as.integer(evaluator$cost$n_theta_miss - n_miss0)
  )
}

.pick_best_theta <- function(explore_df, refined_df, d) {
  th_cols <- paste0("theta", seq_len(d))
  best_q <- Inf
  best_th <- rep(NA_real_, d)
  best_src <- NA_character_
  i_ex <- which.min(explore_df$global_gcv)
  if (length(i_ex) && is.finite(explore_df$global_gcv[i_ex])) {
    best_q <- explore_df$global_gcv[i_ex]
    best_th <- as.numeric(explore_df[i_ex, th_cols])
    best_src <- "explore_best"
  }
  if (!is.null(refined_df) && nrow(refined_df)) {
    i_r <- which.min(refined_df$global_gcv)
    if (length(i_r) && is.finite(refined_df$global_gcv[i_r]) &&
        refined_df$global_gcv[i_r] < best_q - 1e-12) {
      best_q <- refined_df$global_gcv[i_r]
      best_th <- as.numeric(refined_df[i_r, th_cols])
      best_src <- "refined"
    }
  }
  list(theta = best_th, q_search = best_q, source = best_src)
}

# ---------------------------------------------------------------------------
# Reference grid (d=2) and cGCV anchor
# ---------------------------------------------------------------------------

.grid_reference <- function(des, cfg, common_init, probes_search, probes_final,
                            ctrl) {
  lower <- cfg$box[1]
  upper <- cfg$box[2]
  g <- seq(lower, upper, by = cfg$grid_by)
  coords <- expand.grid(t1 = g, t2 = g, KEEP.OUT.ATTRS = FALSE)
  cache <- new.env(parent = emptyenv())
  q <- numeric(nrow(coords))
  n_failed <- 0L
  for (i in seq_len(nrow(coords))) {
    lam <- 10^c(coords$t1[i], coords$t2[i])
    ev <- TTPsplines:::.tt_lab_eval_ggcv(
      lambda = lam, design = des, rank = cfg$rank, probes = probes_search,
      M = cfg$M_search, init_cores = common_init, init_policy = "cold_common",
      mc_bank_id = "bank0", cache = cache, control = ctrl
    )
    q[i] <- ev$global_gcv
    if (!isTRUE(ev$valid) || !is.finite(q[i])) n_failed <- n_failed + 1L
  }
  i_best <- which.min(q)
  th_grid <- c(coords$t1[i_best], coords$t2[i_best])
  ev_f <- TTPsplines:::.tt_lab_eval_ggcv(
    lambda = 10^th_grid, design = des, rank = cfg$rank, probes = probes_final,
    M = cfg$M_final, init_cores = common_init, init_policy = "cold_common",
    mc_bank_id = "gate_final", cache = NULL, control = ctrl
  )
  list(
    theta = th_grid,
    q_search = q[i_best],
    q_final = ev_f$global_gcv,
    n_grid = nrow(coords),
    n_failed = n_failed,
    n_als_approx = as.integer(nrow(coords) * (1L + cfg$M_search)),
    grid_df = data.frame(t1 = coords$t1, t2 = coords$t2, global_gcv = q)
  )
}

# ---------------------------------------------------------------------------
# One method × seed × N
# ---------------------------------------------------------------------------

.run_method <- function(method, des, cfg, N, seed_design,
                        common_init, probes_search, probes_final, ctrl,
                        q_ref_final, theta_ref, verbose = TRUE) {
  d <- des$d
  lower <- rep(cfg$box[1], d)
  upper <- rep(cfg$box[2], d)
  t0 <- proc.time()[["elapsed"]]

  desn <- .make_design(method, N, lower, upper, seed_design)
  evaluator <- TTPsplines:::.tt_lab_make_theta_evaluator(
    y = des$y, X = des$X, rank = cfg$rank, common_init = common_init,
    probes = probes_search, control = ctrl, k = cfg$k,
    degree = des$degree, penalty_order = des$penalty_order,
    lower = lower, upper = upper
  )

  ex <- .eval_points(desn$theta, evaluator, cfg$M_search, method)
  n_explore <- nrow(ex$df)
  n_failed <- ex$n_failed
  refined <- NULL
  n_refine_starts <- 0L
  n_refine_evals <- 0L
  if (isTRUE(desn$do_refine)) {
    rf <- .refine_elites(
      ex$df, evaluator, lower, upper,
      cfg$n_refine, cfg$n_diverse, cfg$min_dist, cfg$M_search
    )
    refined <- rf$refined
    n_refine_starts <- rf$n_refine_starts
    n_refine_evals <- rf$n_refine_evals
  }
  best <- .pick_best_theta(ex$df, refined, d)

  # Fair final re-eval (common bank, cold_common)
  ev_f <- TTPsplines:::.tt_lab_eval_ggcv(
    lambda = 10^best$theta, design = des, rank = cfg$rank,
    probes = probes_final, M = cfg$M_final, init_cores = common_init,
    init_policy = "cold_common", mc_bank_id = "gate_final",
    cache = NULL, control = ctrl
  )
  q_final <- ev_f$global_gcv
  in_valley <- is.finite(q_final) && is.finite(q_ref_final) &&
    q_final <= (1 + cfg$valley_tol) * q_ref_final
  regret <- if (is.finite(q_final) && is.finite(q_ref_final) &&
                abs(q_ref_final) > 1e-15) {
    (q_final - q_ref_final) / q_ref_final
  } else {
    NA_real_
  }
  d_theta <- if (all(is.finite(best$theta)) && all(is.finite(theta_ref))) {
    sqrt(sum((best$theta - theta_ref)^2))
  } else {
    NA_real_
  }

  # Anisotropy (d=2)
  aniso_ok <- NA
  if (d == 2L && all(is.finite(best$theta)) && all(is.finite(theta_ref))) {
    opt_diff <- best$theta[1] - best$theta[2]
    ref_diff <- theta_ref[1] - theta_ref[2]
    iso_opt <- abs(opt_diff) < 0.25
    iso_ref <- abs(ref_diff) < 0.25
    if (identical(des$scenario, "smooth_smooth")) {
      aniso_ok <- iso_opt || (sign(opt_diff) == sign(ref_diff)) ||
        abs(opt_diff - ref_diff) < 0.5
    } else if (identical(des$scenario, "strong_aniso")) {
      aniso_ok <- if (iso_ref) TRUE else
        (sign(opt_diff) == sign(ref_diff) && abs(opt_diff) >= 0.25)
    } else {
      aniso_ok <- sign(opt_diff) == sign(ref_diff) || (iso_opt && iso_ref)
    }
  }

  n_als <- as.integer(evaluator$cost$n_als_fits_approx %||% NA_integer_)
  # Final re-eval ALS approx
  n_als_final <- as.integer(1L + cfg$M_final)
  n_als_total <- n_als + n_als_final
  elapsed <- proc.time()[["elapsed"]] - t0

  if (verbose) {
    message(sprintf(
      "  [%s N=%d seed=%d] valley=%s regret=%.4f n_als=%d t=%.1fs",
      method, N, seed_design, in_valley, regret, n_als_total, elapsed
    ))
  }

  list(
    row = data.frame(
      scenario = des$scenario,
      method = method,
      d = d,
      N = as.integer(N),
      seed_design = as.integer(seed_design),
      seed_data = as.integer(cfg$seed_data),
      theta1 = best$theta[1],
      theta2 = if (d >= 2L) best$theta[2] else NA_real_,
      q_search = best$q_search,
      q_final = q_final,
      q_ref = q_ref_final,
      in_valley_1pct = in_valley,
      regret_rel = regret,
      d_theta = d_theta,
      aniso_ok = aniso_ok,
      winner_source = best$source,
      n_explore = n_explore,
      n_refine_starts = n_refine_starts,
      n_refine_evals = n_refine_evals,
      n_failed = n_failed,
      n_als_search = n_als,
      n_als_total = n_als_total,
      n_theta_miss = as.integer(evaluator$cost$n_theta_miss),
      elapsed = elapsed,
      stringsAsFactors = FALSE
    ),
    explore = ex$df,
    refined = refined
  )
}

.run_cgcv_row <- function(des, cfg, common_init, probes_final, ctrl,
                          q_ref_final, theta_ref) {
  t0 <- proc.time()[["elapsed"]]
  cg <- TTPsplines:::.tt_lab_cgcv_anchor(
    des$y, des$X, cfg$rank, common_init, ctrl, cfg$k,
    des$degree, des$penalty_order
  )
  if (!isTRUE(cg$ok)) {
    return(data.frame(
      scenario = des$scenario, method = "cgcv", d = des$d,
      N = 0L, seed_design = NA_integer_, seed_data = cfg$seed_data,
      theta1 = NA_real_, theta2 = NA_real_,
      q_search = NA_real_, q_final = NA_real_, q_ref = q_ref_final,
      in_valley_1pct = FALSE, regret_rel = NA_real_, d_theta = NA_real_,
      aniso_ok = NA, winner_source = "cgcv_fail",
      n_explore = 0L, n_refine_starts = 0L, n_refine_evals = 0L,
      n_failed = 1L, n_als_search = NA_integer_, n_als_total = NA_integer_,
      n_theta_miss = NA_integer_, elapsed = proc.time()[["elapsed"]] - t0,
      stringsAsFactors = FALSE
    ))
  }
  th <- as.numeric(cg$theta)
  ev_f <- TTPsplines:::.tt_lab_eval_ggcv(
    lambda = 10^th, design = des, rank = cfg$rank, probes = probes_final,
    M = cfg$M_final, init_cores = common_init, init_policy = "cold_common",
    mc_bank_id = "gate_final", cache = NULL, control = ctrl
  )
  q_final <- ev_f$global_gcv
  in_valley <- is.finite(q_final) && is.finite(q_ref_final) &&
    q_final <= (1 + cfg$valley_tol) * q_ref_final
  regret <- if (is.finite(q_final) && is.finite(q_ref_final) &&
                abs(q_ref_final) > 1e-15) {
    (q_final - q_ref_final) / q_ref_final
  } else NA_real_
  d_theta <- sqrt(sum((th - theta_ref)^2))
  data.frame(
    scenario = des$scenario, method = "cgcv", d = des$d,
    N = 0L, seed_design = NA_integer_, seed_data = cfg$seed_data,
    theta1 = th[1], theta2 = th[2],
    q_search = NA_real_, q_final = q_final, q_ref = q_ref_final,
    in_valley_1pct = in_valley, regret_rel = regret, d_theta = d_theta,
    aniso_ok = NA, winner_source = "cgcv",
    n_explore = 0L, n_refine_starts = 0L, n_refine_evals = 0L,
    n_failed = 0L, n_als_search = NA_integer_,
    n_als_total = NA_integer_, n_theta_miss = NA_integer_,
    elapsed = proc.time()[["elapsed"]] - t0,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Aggregate decision helpers
# ---------------------------------------------------------------------------

.summarize_gate <- function(summary_df) {
  # Exclude cgcv from design comparison
  des_df <- summary_df[summary_df$method != "cgcv", , drop = FALSE]
  if (!nrow(des_df)) return(list(table = data.frame(), verdict = "NO_DATA"))

  agg <- aggregate(
    cbind(in_valley_1pct, regret_rel, n_als_total, elapsed) ~ method + N + scenario,
    data = des_df,
    FUN = function(x) mean(as.numeric(x), na.rm = TRUE)
  )
  # success count
  n_seeds <- aggregate(
    in_valley_1pct ~ method + N + scenario,
    data = des_df,
    FUN = function(x) sum(as.logical(x), na.rm = TRUE)
  )
  names(n_seeds)[4] <- "n_valley"
  n_tot <- aggregate(
    in_valley_1pct ~ method + N + scenario,
    data = des_df,
    FUN = length
  )
  names(n_tot)[4] <- "n_seeds"
  tab <- merge(agg, n_seeds, by = c("method", "N", "scenario"))
  tab <- merge(tab, n_tot, by = c("method", "N", "scenario"))
  tab$valley_rate <- tab$n_valley / tab$n_seeds

  # Simple algorithmic verdict on mid/high N (or only available N)
  N_focus <- max(tab$N)
  focus <- tab[tab$N == N_focus, , drop = FALSE]
  sobol_r <- focus[focus$method == "sobol_refine", , drop = FALSE]
  random <- focus[focus$method == "random", , drop = FALSE]
  lhs <- focus[focus$method == "lhs", , drop = FALSE]

  verdict <- "CONDITIONAL GO"
  notes <- character(0)
  if (nrow(sobol_r)) {
    # Across scenarios: require mean valley rate >= 0.8 (4/5 style)
    rate_s <- mean(sobol_r$valley_rate)
    rate_r <- if (nrow(random)) mean(random$valley_rate) else NA_real_
    rate_l <- if (nrow(lhs)) mean(lhs$valley_rate) else NA_real_
    als_s <- mean(sobol_r$n_als_total)
    als_r <- if (nrow(random)) mean(random$n_als_total) else NA_real_
    if (rate_s >= 0.8 &&
        (is.na(rate_r) || rate_s > rate_r + 1e-9 ||
         (!is.na(als_r) && als_s < als_r)) &&
        (is.na(rate_l) || rate_s > rate_l + 1e-9)) {
      verdict <- "GO"
      notes <- sprintf(
        "sobol_refine valley_rate=%.2f vs random=%.2f lhs=%.2f at N=%d",
        rate_s, rate_r, rate_l, N_focus
      )
    } else if (rate_s >= 0.8) {
      verdict <- "CONDITIONAL GO"
      notes <- sprintf(
        "sobol_refine valley_rate=%.2f but not clearly better than random/LHS",
        rate_s
      )
    } else if (!is.na(rate_r) && abs(rate_s - rate_r) < 0.05 &&
               !is.na(rate_l) && abs(rate_s - rate_l) < 0.05) {
      verdict <- "NO-GO"
      notes <- "random/LHS equivalent valley rates"
    } else {
      verdict <- "CONDITIONAL GO"
      notes <- sprintf("sobol_refine valley_rate=%.2f < 0.8 threshold", rate_s)
    }
  }
  list(table = tab, verdict = verdict, notes = paste(notes, collapse = "; "))
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

cfg <- .cfg_apply_env(.cfg_for_mode(mode))
message("MODE=", mode)
message("methods=", paste(cfg$methods, collapse = ","))
message("scenarios=", paste(cfg$scenarios, collapse = ","))
message("n_mult=", paste(cfg$n_mult, collapse = ","))
message("seeds=", paste(cfg$seeds_design, collapse = ","))

summary_rows <- list()
ref_rows <- list()

for (sc in cfg$scenarios) {
  message("=== scenario=", sc, " ===")
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    scenario = sc, n = cfg$n, n_test = cfg$n_test, k = cfg$k,
    seed = cfg$seed_data
  )
  d <- des$d
  ctrl <- tt_control(
    max_sweeps = 40L, tol = 1e-8, compute_edf = FALSE, seed = cfg$seed_data
  )
  # Fixed cores + MC banks (seed_data only — design seed does not touch these)
  common_init <- TTPsplines:::.tt_with_preserved_seed({
    tt_initialize(d = d, rank = cfg$rank, k = cfg$k, seed = cfg$seed_data)
  })
  M_bank <- max(cfg$M_search, cfg$M_final)
  probes_bank <- TTPsplines:::.tt_lab_rademacher_probes(
    des$n, M_bank, probe_seed = cfg$seed_data + 11L
  )
  probes_search <- probes_bank[, seq_len(cfg$M_search), drop = FALSE]
  probes_final <- probes_bank[, seq_len(cfg$M_final), drop = FALSE]

  # Reference
  if (isTRUE(cfg$include_grid_ref) && d == 2L) {
    message("  building grid reference...")
    ref <- .grid_reference(
      des, cfg, common_init, probes_search, probes_final, ctrl
    )
    q_ref <- ref$q_final
    theta_ref <- ref$theta
    utils::write.csv(
      ref$grid_df,
      file.path(out_dir, paste0(sc, "_grid.csv")),
      row.names = FALSE
    )
    ref_rows[[sc]] <- data.frame(
      scenario = sc,
      theta1 = theta_ref[1], theta2 = theta_ref[2],
      q_ref = q_ref, n_grid = ref$n_grid,
      n_als_grid = ref$n_als_approx,
      stringsAsFactors = FALSE
    )
    message(sprintf(
      "  ref theta=(%.2f,%.2f) Q=%.6g n_grid=%d",
      theta_ref[1], theta_ref[2], q_ref, ref$n_grid
    ))
  } else {
    stop("Matched-budget v0 requires d=2 grid reference.", call. = FALSE)
  }

  if (isTRUE(cfg$include_cgcv)) {
    message("  cGCV baseline...")
    cg_row <- .run_cgcv_row(
      des, cfg, common_init, probes_final, ctrl, q_ref, theta_ref
    )
    summary_rows[[paste(sc, "cgcv", sep = "|")]] <- cg_row
    message(sprintf(
      "  [cgcv] valley=%s regret=%.4f",
      cg_row$in_valley_1pct, cg_row$regret_rel
    ))
  }

  for (mult in cfg$n_mult) {
    N <- as.integer(mult) * d
    for (seed_d in cfg$seeds_design) {
      for (meth in cfg$methods) {
        key <- paste(sc, meth, N, seed_d, sep = "|")
        message(sprintf("--- %s ---", key))
        res <- .run_method(
          method = meth, des = des, cfg = cfg, N = N,
          seed_design = seed_d,
          common_init = common_init,
          probes_search = probes_search,
          probes_final = probes_final,
          ctrl = ctrl,
          q_ref_final = q_ref,
          theta_ref = theta_ref,
          verbose = TRUE
        )
        summary_rows[[key]] <- res$row
        utils::write.csv(
          res$explore,
          file.path(out_dir, sprintf("%s_%s_N%d_s%d_explore.csv",
                                     sc, meth, N, seed_d)),
          row.names = FALSE
        )
        if (!is.null(res$refined) && nrow(res$refined)) {
          utils::write.csv(
            res$refined,
            file.path(out_dir, sprintf("%s_%s_N%d_s%d_refined.csv",
                                       sc, meth, N, seed_d)),
            row.names = FALSE
          )
        }
      }
    }
  }
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

# Empirical reference within the run: min over grid + all finite finals.
# Avoids false FAIL when a coarse grid is beaten by cGCV or a design point.
for (sc in unique(summary_df$scenario)) {
  idx <- which(summary_df$scenario == sc)
  q_emp <- min(summary_df$q_final[idx], na.rm = TRUE)
  if (is.finite(q_emp)) {
    # Prefer documented grid ref if it is better or equal; else lift to empirical
    q_grid <- summary_df$q_ref[idx[1]]
    q_use <- if (is.finite(q_grid)) min(q_grid, q_emp) else q_emp
    summary_df$q_ref[idx] <- q_use
    summary_df$regret_rel[idx] <- (summary_df$q_final[idx] - q_use) / q_use
    summary_df$in_valley_1pct[idx] <- is.finite(summary_df$q_final[idx]) &
      summary_df$q_final[idx] <= (1 + cfg$valley_tol) * q_use
  }
}

utils::write.csv(summary_df, file.path(out_dir, "summary.csv"), row.names = FALSE)

if (length(ref_rows)) {
  ref_df <- do.call(rbind, ref_rows)
  utils::write.csv(ref_df, file.path(out_dir, "reference.csv"), row.names = FALSE)
}

gate <- .summarize_gate(summary_df)
utils::write.csv(gate$table, file.path(out_dir, "aggregate.csv"), row.names = FALSE)

# Console + short report
message("\n======== MATCHED-BUDGET SUMMARY ========")
print(gate$table)
message("VERDICT: ", gate$verdict)
message("NOTES: ", gate$notes)

report_path <- file.path(out_dir, "MATCHED_BUDGET_REPORT.md")
sink(report_path)
cat("# Matched-budget gate report\n\n")
cat("**Mode:** ", mode, "\n\n", sep = "")
cat("**Date:** ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n", sep = "")
cat("**Verdict:** ", gate$verdict, "\n\n", sep = "")
cat(gate$notes, "\n\n", sep = "")
cat("Product lock: KEEP cGCV. Protocol: `PROTOCOL_MATCHED_BUDGET.md`.\n\n")
cat("## Aggregate (mean over design seeds)\n\n")
print(gate$table)
cat("\n## Full summary\n\n")
print(summary_df)
sink()

message("Wrote ", out_dir)
message("Report: ", report_path)
