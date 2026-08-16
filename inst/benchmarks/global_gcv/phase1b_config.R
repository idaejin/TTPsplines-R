# Phase 1B FULL selective configuration (does NOT write into results/phase1/).
# Env:
#   TT_GGCV_MODE=SMOKE|FULL   (default FULL if QUICK=false else SMOKE)
#   TT_GGCV_SCENARIOS=smooth_smooth,strong_aniso
#   TT_GGCV_R=50
#   TT_GGCV_M_SEARCH / TT_GGCV_M_FINAL / TT_GGCV_M_BANK
#   TT_GGCV_DO_NELDER=true|false
#   TT_GGCV_DO_RESTRICTED=true|false

.phase1b_this_dir <- function() {
  # When sourced, prefer the sourcing file path
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    return(dirname(normalizePath(ofile)))
  }
  # Fallback: walk from getwd()
  cand <- c(
    file.path("inst", "benchmarks", "global_gcv"),
    file.path("..", "..", "inst", "benchmarks", "global_gcv"),
    "."
  )
  for (p in cand) {
    if (file.exists(file.path(p, "phase1_config.R"))) return(normalizePath(p))
  }
  normalizePath(".")
}

.source_phase1_config <- function() {
  d <- .phase1b_this_dir()
  cfg_path <- file.path(d, "phase1_config.R")
  if (!file.exists(cfg_path)) {
    # When phase1b_config lives in inst/... and ofile works:
    cfg_path <- file.path("inst", "benchmarks", "global_gcv", "phase1_config.R")
  }
  if (!file.exists(cfg_path)) {
    # From tests/testthat
    cfg_path <- file.path("..", "..", "inst", "benchmarks", "global_gcv", "phase1_config.R")
  }
  source(cfg_path, local = FALSE)
}

.source_phase1_config()

.phase1b_parse_scenarios <- function(default = c("smooth_smooth", "strong_aniso")) {
  raw <- Sys.getenv("TT_GGCV_SCENARIOS", unset = "")
  if (!nzchar(raw)) return(default)
  sc <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  sc <- sc[nzchar(sc)]
  allowed <- c("smooth_smooth", "smooth_rough", "strong_aniso", "weak_signal")
  bad <- setdiff(sc, allowed)
  if (length(bad)) stop("Unknown scenarios: ", paste(bad, collapse = ", "), call. = FALSE)
  sc
}

.phase1b_config <- function() {
  mode_env <- toupper(Sys.getenv("TT_GGCV_MODE", unset = ""))
  if (!nzchar(mode_env)) {
    mode_env <- if (.tt_ggcv_env_flag("TT_GGCV_QUICK", FALSE)) "SMOKE" else "FULL"
  }
  mode_env <- match.arg(mode_env, c("SMOKE", "FULL"))
  smoke <- identical(mode_env, "SMOKE")

  # QUICK phase1 paths must never be used as write targets here
  out_dir <- file.path("inst", "benchmarks", "global_gcv", "results",
                       if (smoke) "phase1b_smoke" else "phase1b_full")
  fig_dir <- file.path("inst", "benchmarks", "global_gcv", "figures",
                       if (smoke) "phase1b_smoke" else "phase1b_full")

  list(
    mode = mode_env,
    smoke = smoke,
    scenarios = .phase1b_parse_scenarios(),
    n_rep = .tt_ggcv_env_int("TT_GGCV_R", if (smoke) 2L else 50L),
    M_search = .tt_ggcv_env_int("TT_GGCV_M_SEARCH", if (smoke) 3L else 5L),
    M_final = .tt_ggcv_env_int("TT_GGCV_M_FINAL", if (smoke) 8L else 20L),
    M_bank = .tt_ggcv_env_int("TT_GGCV_M_BANK", if (smoke) 20L else 40L),
    epsilon_rel = .tt_ggcv_env_num("TT_GGCV_EPS_REL", 1e-3),
    # Keep DGP sizes compatible with Phase 1 FULL defaults (not QUICK n=120)
    n = .tt_ggcv_env_int("TT_GGCV_N", if (smoke) 120L else 200L),
    n_test = .tt_ggcv_env_int("TT_GGCV_NTEST", if (smoke) 150L else 400L),
    k = .tt_ggcv_env_int("TT_GGCV_K", if (smoke) 6L else 8L),
    seed0 = .tt_ggcv_env_int("TT_GGCV_SEED0", 20260816L),
    max_sweeps = if (smoke) 25L else 40L,
    # Refined grid around QUICK interior centers
    refine_by = if (smoke) 1 else 0.5,
    refine_halfwidth = if (smoke) 1 else 1.5,
    widen_once = TRUE,
    do_nelder = .tt_ggcv_env_flag("TT_GGCV_DO_NELDER", !smoke),
    nelder_max_reps = .tt_ggcv_env_int("TT_GGCV_NELDER_MAX_REPS", if (smoke) 1L else 10L),
    do_restricted = .tt_ggcv_env_flag("TT_GGCV_DO_RESTRICTED", TRUE),
    restricted_max_reps = .tt_ggcv_env_int("TT_GGCV_RESTRICTED_MAX_REPS", if (smoke) 1L else 10L),
    # QUICK centers (log10) from phase1_designs / surface minima
    centers = list(
      smooth_smooth = c(-1, -3),
      strong_aniso = c(3, 3)
    ),
    # Expected anisotropy sign for strong_aniso: margin2 needs less smooth
    # (higher roughness) → typically lambda2 < lambda1 → log10(l1/l2) > 0
    aniso_sign_truth = list(
      smooth_smooth = 0,          # roughly isotropic preference
      strong_aniso = 1            # expect log10(l1/l2) > 0
    ),
    out_dir = out_dir,
    fig_dir = fig_dir,
    quick_dir = file.path("inst", "benchmarks", "global_gcv", "results", "phase1")
  )
}

.phase1b_ensure_dirs <- function(cfg) {
  # Refuse to write into QUICK phase1 directory
  norm <- gsub("\\\\", "/", cfg$out_dir)
  if (grepl("/results/phase1$", norm) || grepl("/results/phase1/", norm)) {
    stop("Refusing to write into QUICK results/phase1 — use phase1b_* paths.",
         call. = FALSE)
  }
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

.phase1b_refine_grid <- function(center, halfwidth, by, box = c(-5, 5)) {
  center <- as.numeric(center)
  g1 <- seq(center[1] - halfwidth, center[1] + halfwidth, by = by)
  g2 <- seq(center[2] - halfwidth, center[2] + halfwidth, by = by)
  g1 <- g1[g1 >= box[1] & g1 <= box[2]]
  g2 <- g2[g2 >= box[1] & g2 <= box[2]]
  list(g1 = unique(round(g1, 10)), g2 = unique(round(g2, 10)))
}
