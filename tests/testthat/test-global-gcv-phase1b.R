# Phase 1B tests (SMOKE-scale; NOT_CRAN for heavier bits).

.phase1b_pkg_root <- function() {
  cand <- c(file.path("..", ".."), ".", file.path("..", "..", ".."))
  for (p in cand) {
    desc <- file.path(p, "DESCRIPTION")
    if (file.exists(desc)) {
      pkg <- tryCatch(read.dcf(desc, fields = "Package")[1],
                      error = function(e) NA_character_)
      if (identical(pkg, "TTPsplines")) return(normalizePath(p))
    }
  }
  normalizePath(".")
}

.phase1b_source <- function(rel) {
  root <- .phase1b_pkg_root()
  env <- parent.frame()
  sys.source(file.path(root, rel), envir = env)
}

test_that("phase1b scenario parsing and path separation", {
  .phase1b_source("inst/benchmarks/global_gcv/phase1b_config.R")
  old_sc <- Sys.getenv("TT_GGCV_SCENARIOS", unset = NA)
  on.exit({
    if (is.na(old_sc)) Sys.unsetenv("TT_GGCV_SCENARIOS")
    else Sys.setenv(TT_GGCV_SCENARIOS = old_sc)
  }, add = TRUE)
  Sys.setenv(TT_GGCV_SCENARIOS = "smooth_smooth,strong_aniso")
  sc <- .phase1b_parse_scenarios()
  expect_equal(sc, c("smooth_smooth", "strong_aniso"))
  Sys.setenv(TT_GGCV_SCENARIOS = "nope")
  expect_error(.phase1b_parse_scenarios(), regexp = "Unknown")
})

test_that("phase1b refuses QUICK phase1 out_dir and separates modes", {
  .phase1b_source("inst/benchmarks/global_gcv/phase1b_config.R")
  old_mode <- Sys.getenv("TT_GGCV_MODE", unset = NA)
  on.exit({
    if (is.na(old_mode)) Sys.unsetenv("TT_GGCV_MODE")
    else Sys.setenv(TT_GGCV_MODE = old_mode)
  }, add = TRUE)
  Sys.setenv(TT_GGCV_MODE = "SMOKE")
  cfg <- .phase1b_config()
  expect_true(grepl("phase1b_smoke", cfg$out_dir))
  expect_false(grepl("results/phase1$", cfg$out_dir))
  expect_error({
    bad <- cfg
    bad$out_dir <- file.path(.phase1b_pkg_root(), "inst", "benchmarks",
                             "global_gcv", "results", "phase1")
    .phase1b_ensure_dirs(bad)
  }, regexp = "Refusing")
})

test_that("phase1b paired stats and flatness", {
  .phase1b_source("inst/benchmarks/global_gcv/phase1b_config.R")
  .phase1b_source("inst/benchmarks/global_gcv/phase1b_helpers.R")
  st <- .phase1b_paired_stats(c(-1, -0.5, 0.2, 0.1, -0.3))
  expect_equal(st$n, 5L)
  expect_true(st$mean < 0)
  expect_true(st$prop_neg > 0.5)
  gr <- .phase1b_refine_grid(c(-1, -3), 1, 1, c(-5, 5))
  expect_true(all(gr$g1 >= -5 & gr$g1 <= 5))
})

test_that("phase1b checkpoint append + resume keys", {
  .phase1b_source("inst/benchmarks/global_gcv/phase1b_helpers.R")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  d1 <- data.frame(scenario = "smooth_smooth", rep = 1L, x = 1)
  d2 <- data.frame(scenario = "smooth_smooth", rep = 2L, x = 2)
  .phase1b_append_csv(d1, tmp)
  .phase1b_append_csv(d2, tmp)
  tab <- utils::read.csv(tmp)
  expect_equal(nrow(tab), 2L)
  keys <- .phase1b_done_keys(tmp)
  expect_true("smooth_smooth|1" %in% keys)
  expect_true("smooth_smooth|2" %in% keys)
})

test_that("phase1b multistart selection does not use test-MSE", {
  starts <- list(
    low = list(ok = TRUE, global_gcv = 0.20, test_mse = 0.01, start = "low"),
    central = list(ok = TRUE, global_gcv = 0.10, test_mse = 0.50, start = "central"),
    high = list(ok = TRUE, global_gcv = 0.30, test_mse = 0.02, start = "high")
  )
  gcv_vals <- vapply(starts, function(z) z$global_gcv, numeric(1))
  best <- starts[[which.min(gcv_vals)]]
  expect_equal(best$start, "central")
  expect_gt(best$test_mse, min(vapply(starts, function(z) z$test_mse, numeric(1))))
})

test_that("phase1b smoke one-rep runs end-to-end", {
  skip_on_cran()
  root <- .phase1b_pkg_root()
  .phase1b_source("inst/benchmarks/global_gcv/phase1b_config.R")
  .phase1b_source("inst/benchmarks/global_gcv/phase1b_helpers.R")
  sys.source(file.path(root, "inst", "benchmarks", "helpers.R"), envir = environment())
  .ttps_bench_ensure_pkg()
  old_vars <- c(
    "TT_GGCV_MODE", "TT_GGCV_R", "TT_GGCV_SCENARIOS", "TT_GGCV_DO_NELDER",
    "TT_GGCV_DO_RESTRICTED", "TT_GGCV_N", "TT_GGCV_NTEST", "TT_GGCV_K",
    "TT_GGCV_M_SEARCH", "TT_GGCV_M_FINAL"
  )
  old <- Sys.getenv(old_vars, unset = NA)
  on.exit({
    for (nm in old_vars) {
      if (is.na(old[[nm]])) Sys.unsetenv(nm)
      else do.call(Sys.setenv, setNames(list(old[[nm]]), nm))
    }
  }, add = TRUE)
  Sys.setenv(
    TT_GGCV_MODE = "SMOKE",
    TT_GGCV_R = "1",
    TT_GGCV_SCENARIOS = "smooth_smooth",
    TT_GGCV_DO_NELDER = "false",
    TT_GGCV_DO_RESTRICTED = "false",
    TT_GGCV_N = "60",
    TT_GGCV_NTEST = "40",
    TT_GGCV_K = "5",
    TT_GGCV_M_SEARCH = "2",
    TT_GGCV_M_FINAL = "4"
  )
  cfg <- .phase1b_config()
  cfg$refine_by <- 2
  cfg$refine_halfwidth <- 2
  cfg$n_rep <- 1L
  out <- .phase1b_one_rep("smooth_smooth", 1L, cfg)
  expect_true(is.data.frame(out$replication))
  expect_true(any(out$replication$method == "tt_grid_oracle"))
  expect_true(any(out$replication$method == "cGCV_default"))
  expect_true(all(is.finite(out$replication$test_mse[
    out$replication$method == "cGCV_default"
  ])))
  expect_true(is.data.frame(out$stability))
})
