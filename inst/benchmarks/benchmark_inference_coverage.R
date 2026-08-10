# Small coverage pilot for conditional Bayesian / frequentist TT intervals.
# Gate 1: Gaussian, fixed rank, fixed lambda. Pointwise only.
#
# Run: Rscript inst/benchmarks/benchmark_inference_coverage.R

suppressPackageStartupMessages({
  pkg_root <- if (file.exists("DESCRIPTION")) {
    "."
  } else if (file.exists("../../DESCRIPTION")) {
    "../.."
  } else {
    stop("Run from package root: Rscript inst/benchmarks/benchmark_inference_coverage.R")
  }
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    library(TTPsplines)
  }
})

.near_null_f <- function(X) {
  # Nearly in the 2nd-order difference null space: affine surface
  0.5 + 1.2 * X[, 1] - 0.8 * X[, 2] + if (ncol(X) >= 3) 0.3 * X[, 3] else 0
}

.run_one <- function(dgp, seed, n = 800L, rank = NULL, k = NULL, lambda = NULL,
                     sigma = NULL, n_test = 200L, level = 0.95) {
  # Defaults chosen so TT approximation bias is modest (coverage is not
  # meaningful for a grossly under-ranked smoother).
  if (identical(dgp, "ishigami")) {
    rank <- rank %||% 5L; k <- k %||% 8L; lambda <- lambda %||% 0.3; sigma <- sigma %||% 0.35
    sim <- simulate_ishigami(n = n, sigma = sigma, seed = seed)
    X <- sim$X; y <- sim$y
    f_true_fun <- function(Xnew) f_ishigami(Xnew)
  } else if (identical(dgp, "friedman")) {
    rank <- rank %||% 5L; k <- k %||% 8L; lambda <- lambda %||% 0.5; sigma <- sigma %||% 1
    sim <- simulate_friedman(n = n, sigma = sigma, seed = seed)
    X <- sim$X; y <- sim$y
    f_true_fun <- function(Xnew) f_friedman(Xnew)
  } else if (identical(dgp, "near_null")) {
    # Mild penalty so we are near but not deep in the null-space undercoverage regime
    rank <- rank %||% 2L; k <- k %||% 6L; lambda <- lambda %||% 1; sigma <- sigma %||% 0.35
    d <- 3L
    set.seed(seed)
    X <- matrix(runif(n * d), n, d)
    f <- .near_null_f(X)
    y <- f + rnorm(n, 0, sigma)
    f_true_fun <- .near_null_f
  } else {
    stop("Unknown dgp")
  }

  t_fit0 <- proc.time()[[3L]]
  fit <- ttps(
    y, X, family = gaussian(), rank = rank, k = k, lambda = lambda,
    control = tt_control(max_sweeps = 12, compute_edf = TRUE, seed = seed)
  )
  t_fit <- proc.time()[[3L]] - t_fit0

  set.seed(seed + 1000L)
  Xte <- matrix(runif(n_test * ncol(X)), n_test, ncol(X))
  # Match training domain roughly
  for (j in seq_len(ncol(X))) {
    Xte[, j] <- stats::runif(n_test, min(X[, j]), max(X[, j]))
  }
  truth <- f_true_fun(Xte)

  t_cov0 <- proc.time()[[3L]]
  invisible(vcov(fit, type = "bayesian"))
  t_cov <- proc.time()[[3L]] - t_cov0

  out <- list()
  for (ctype in c("bayesian", "frequentist")) {
    t_se0 <- proc.time()[[3L]]
    pr <- predict(
      fit, Xte, type = "response", se.fit = TRUE,
      interval = "confidence", level = level, vcov_type = ctype
    )
    t_se <- proc.time()[[3L]] - t_se0
    cover <- mean(truth >= pr$lower & truth <= pr$upper)
    width <- mean(pr$upper - pr$lower)
    rmse <- sqrt(mean((pr$fit - truth)^2))
    out[[ctype]] <- data.frame(
      dgp = dgp, covariance = ctype, seed = seed, n = n, rank = rank,
      k = k, lambda = lambda, sigma = sigma, n_test = n_test,
      nominal = level, coverage = cover, mean_width = width, rmse = rmse,
      fit_time = t_fit, cov_setup_time = t_cov, pred_se_time = t_se,
      edf = fit$edf %||% NA_real_,
      sigma2_hat = fit$._inf$data$scale %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

seeds <- 1:5
dgps <- c("ishigami", "friedman", "near_null")
rows <- list()
idx <- 1L
for (dgp in dgps) {
  for (s in seeds) {
    cat(sprintf("Running %s seed=%d ...\n", dgp, s))
    rows[[idx]] <- .run_one(dgp, seed = s)
    idx <- idx + 1L
  }
}
tab <- do.call(rbind, rows)
agg <- aggregate(
  cbind(coverage, mean_width, rmse, fit_time, cov_setup_time, pred_se_time) ~ dgp + covariance,
  data = tab, FUN = mean
)
agg$nominal <- 0.95

out_dir <- file.path("..", "..", "inst", "benchmarks", "results")
# script may run from package root
if (!dir.exists("inst/benchmarks")) {
  # already in package via load_all path
}
res_dir <- "inst/benchmarks/results"
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
utils::write.csv(tab, file.path(res_dir, "inference_coverage_pilot_raw.csv"),
                 row.names = FALSE)
utils::write.csv(agg, file.path(res_dir, "inference_coverage_pilot_summary.csv"),
                 row.names = FALSE)

cat("\n=== Coverage pilot summary (mean over seeds) ===\n")
print(agg)
cat("\nWrote:\n",
    file.path(res_dir, "inference_coverage_pilot_summary.csv"), "\n")
