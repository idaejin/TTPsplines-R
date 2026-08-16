# Phase-1 lab configuration (env overrides).
# Sourced by phase1_*.R scripts.

.tt_ggcv_env_flag <- function(name, default = FALSE) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) return(isTRUE(default))
  tolower(v) %in% c("1", "true", "yes", "y", "t")
}

.tt_ggcv_env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) return(as.integer(default))
  as.integer(v)
}

.tt_ggcv_env_num <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) return(as.numeric(default))
  as.numeric(v)
}

.phase1_config <- function() {
  quick <- .tt_ggcv_env_flag("TT_GGCV_QUICK", TRUE)
  list(
    quick = quick,
    # Grids
    coarse_by = if (quick) 2 else 1,
    coarse_lo = -5,
    coarse_hi = 5,
    refine_halfwidth = if (quick) 1 else 1.5,
    refine_by = if (quick) 0.5 else 0.25,
    # MC
    M_search = .tt_ggcv_env_int("TT_GGCV_M_SEARCH", if (quick) 3L else 5L),
    M_final = .tt_ggcv_env_int("TT_GGCV_M_FINAL", if (quick) 8L else 20L),
    M_bank = .tt_ggcv_env_int("TT_GGCV_M_BANK", if (quick) 20L else 40L),
    epsilon_rel = .tt_ggcv_env_num("TT_GGCV_EPS_REL", 1e-3),
    # Design sizes
    n = .tt_ggcv_env_int("TT_GGCV_N", if (quick) 120L else 200L),
    n_test = .tt_ggcv_env_int("TT_GGCV_NTEST", if (quick) 120L else 200L),
    k = .tt_ggcv_env_int("TT_GGCV_K", if (quick) 6L else 8L),
    # Replicates (diagnostic = 1 seed; compact MC uses more)
    n_rep_diag = .tt_ggcv_env_int("TT_GGCV_R_DIAG", 1L),
    n_rep_mc = .tt_ggcv_env_int("TT_GGCV_R", if (quick) 2L else 5L),
    seed0 = .tt_ggcv_env_int("TT_GGCV_SEED0", 20260816L),
    # Ranks: restricted / mid / sufficient (=k for d=2)
    max_sweeps = if (quick) 25L else 40L,
    scenarios = c("smooth_smooth", "smooth_rough", "strong_aniso", "weak_signal"),
    out_dir = file.path("inst", "benchmarks", "global_gcv", "results", "phase1"),
    fig_dir = file.path("inst", "benchmarks", "global_gcv", "figures", "phase1")
  )
}

.phase1_ensure_dirs <- function(cfg) {
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)
}

.phase1_checkpoint_write <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}

.phase1_load_pkg <- function() {
  if (!exists(".ttps_bench_ensure_pkg", mode = "function")) {
    source(file.path("inst", "benchmarks", "helpers.R"))
  }
  .ttps_bench_ensure_pkg()
  invisible(TRUE)
}

.phase1_ranks <- function(k) {
  k <- as.integer(k)
  list(
    restricted = 1L,
    mid = max(2L, as.integer(ceiling(k / 2))),
    sufficient = k
  )
}

.phase1_probe_bank <- function(n, M_bank, seed) {
  TTPsplines:::.tt_lab_rademacher_probes(n, M = as.integer(M_bank), probe_seed = as.integer(seed))
}
