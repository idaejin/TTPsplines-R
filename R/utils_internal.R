# Small internal utilities shared across files.
#
# .tt_with_preserved_seed() keeps the caller's RNG stream intact around code
# that reseeds (ttps() core initialization, probe draws); it is used by the
# gGCV engine, the gGCV dispatch and tt_geometry.R. .tt_clone_cores() deep-
# copies TT cores (warm-started gGCV probes). .tt_lab_rademacher_probes()
# keeps its old name because lab simulation scripts call it through
# TTPsplines::: (sc_main/run_selector_compare.R, ijoc_exp_B/run_gcv_surface_d2.R,
# ijoc_exp_B/run_ext_pack.R); the package itself does not use it.

#' Save / restore .Random.seed around an expression (package has no withr).
#' @keywords internal
#' @noRd
.tt_with_preserved_seed <- function(expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

#' Deep-copy a list of TT core arrays.
#' @keywords internal
#' @noRd
.tt_clone_cores <- function(cores) {
  lapply(cores, function(a) {
    out <- array(as.numeric(a), dim = dim(a))
    storage.mode(out) <- "double"
    out
  })
}

#' Rademacher probe matrix (n x M) with seed isolation.
#'
#' Kept under its old name for lab scripts that call it through
#' `TTPsplines:::`; the gGCV engine draws its probes with `.ggcv_probe()`.
#' @keywords internal
#' @noRd
.tt_lab_rademacher_probes <- function(n, M, probe_seed = 1L) {
  n <- as.integer(n)
  M <- as.integer(M)
  stopifnot(n >= 1L, M >= 1L)
  .tt_with_preserved_seed({
    set.seed(as.integer(probe_seed))
    mat <- matrix(
      sample(c(-1, 1), size = n * M, replace = TRUE),
      nrow = n, ncol = M
    )
  })
  storage.mode(mat) <- "double"
  mat
}
