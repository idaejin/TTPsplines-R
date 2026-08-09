# Shared helpers for TTPsplines benchmarks.
# Sourced by benchmark_*.R — not run by R CMD check.

.ttps_bench_pkg_root <- function() {
  # When sourced from package install:
  root <- tryCatch(
    system.file(package = "TTPsplines"),
    error = function(e) ""
  )
  if (nzchar(root) && dir.exists(file.path(root, "benchmarks"))) {
    return(root)
  }
  # Development: this file lives in inst/benchmarks/
  normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "../.."),
                mustWork = FALSE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

.ttps_bench_ensure_pkg <- function() {
  if (requireNamespace("TTPsplines", quietly = TRUE) &&
      "ttpspline" %in% getNamespaceExports("TTPsplines")) {
    suppressPackageStartupMessages(library(TTPsplines))
    return(invisible(TRUE))
  }
  # Walk up from CWD for DESCRIPTION (devtools load_all)
  cand <- normalizePath(getwd())
  for (i in 1:8) {
    if (file.exists(file.path(cand, "DESCRIPTION"))) {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        stop("Install/load TTPsplines, or use devtools::load_all()", call. = FALSE)
      }
      devtools::load_all(cand, quiet = TRUE)
      return(invisible(TRUE))
    }
    parent <- dirname(cand)
    if (identical(parent, cand)) break
    cand <- parent
  }
  stop("TTPsplines not loaded. Run: devtools::load_all('ttpsplines-pkg')", call. = FALSE)
}

.ttps_bench_outdir <- function() {
  out <- Sys.getenv("TTPSPLINES_BENCH_OUT", unset = "")
  if (!nzchar(out)) {
    # Prefer package inst/benchmarks/results when developing
    here <- getwd()
    if (basename(here) == "benchmarks") {
      out <- file.path(here, "results")
    } else if (dir.exists(file.path(here, "inst", "benchmarks"))) {
      out <- file.path(here, "inst", "benchmarks", "results")
    } else {
      out <- file.path(tempdir(), "ttpsplines_bench")
    }
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  normalizePath(out)
}

# ----------------------------------------------------------------------
# Non-additive truth (scattered X on [0,1]^d)
# ----------------------------------------------------------------------

#' Genuinely interactive smooth on [0,1]^d (d >= 2).
true_surface_nd <- function(X, amp = 1) {
  X <- as.matrix(X)
  d <- ncol(X)
  stopifnot(d >= 2L)
  f <- amp * (
    sin(2 * pi * X[, 1]) * cos(2 * pi * X[, 2]) +
      0.5 * sin(2 * pi * X[, min(3L, d)]) +
      0.35 * (X[, 1] - 0.5) * (X[, 2] - 0.5)
  )
  if (d >= 4L) {
    f <- f + amp * 0.25 * sin(2 * pi * X[, 4]) * cos(2 * pi * X[, 2])
  }
  f
}

simulate_gaussian <- function(n = 1500L, d = 3L, sigma = 0.3, seed = 1L) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  colnames(X) <- paste0("x", seq_len(d))
  truth <- true_surface_nd(X)
  list(X = X, truth = truth, y = truth + rnorm(n, 0, sigma), family = "gaussian")
}

simulate_poisson <- function(n = 1500L, d = 3L, seed = 2L) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  colnames(X) <- paste0("x", seq_len(d))
  raw <- true_surface_nd(X, amp = 1)
  shift <- mean(raw)
  eta <- raw - shift + log(3)
  mu <- exp(eta)
  list(
    X = X, truth_eta = eta, truth_mu = mu, y = rpois(n, mu),
    family = "poisson", shift = shift
  )
}

simulate_bernoulli <- function(n = 2000L, d = 3L, seed = 3L) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  colnames(X) <- paste0("x", seq_len(d))
  raw <- true_surface_nd(X, amp = 1.8)
  shift <- mean(raw)
  eta <- raw - shift
  p <- plogis(eta)
  list(
    X = X, truth_eta = eta, truth_p = p, y = rbinom(n, 1L, p),
    family = "bernoulli", shift = shift
  )
}

holdout_poisson <- function(n_te = 4000L, d = 3L, shift, seed = 99L) {
  set.seed(seed)
  X <- matrix(runif(n_te * d), n_te, d)
  eta <- true_surface_nd(X) - shift + log(3)
  list(X = X, truth_eta = eta, truth_mu = exp(eta), shift = shift)
}

holdout_bernoulli <- function(n_te = 4000L, d = 3L, shift, seed = 99L) {
  set.seed(seed)
  X <- matrix(runif(n_te * d), n_te, d)
  eta <- true_surface_nd(X, amp = 1.8) - shift
  list(X = X, truth_eta = eta, truth_p = plogis(eta), shift = shift)
}

# Hold-out from same DGP
holdout_gaussian <- function(n_te = 4000L, d = 3L, seed = 99L) {
  set.seed(seed)
  X <- matrix(runif(n_te * d), n_te, d)
  list(X = X, truth = true_surface_nd(X))
}

rmse <- function(a, b) sqrt(mean((as.numeric(a) - as.numeric(b))^2))

.logloss <- function(y, p) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

.write_bench_csv <- function(df, name, out_dir) {
  path <- file.path(out_dir, name)
  write.csv(df, path, row.names = FALSE)
  message("Wrote ", path)
  invisible(path)
}

.default_control <- function(backend = "auto", ...) {
  dots <- list(...)
  base <- list(
    max_sweeps = 15L,
    pirls_maxit = 20L,
    als_sweeps_per_pirls = 3L,
    tol = 1e-6,
    tol_lambda = 1e-3,
    lambda_bounds = c(1e-2, 1e2),
    backend = backend,
    trace = FALSE,
    seed = 1L
  )
  do.call(tt_control, utils::modifyList(base, dots))
}
