# Chicago mortality: 5-way interaction without materializing the full tensor.
#
# Contrast with mgcv::bam:
#   bam(death ~ te(tmpd, o3median, pm10median, so2median, time,
#                  bs = "ps", k = c(5,5,5,5,8)),  # 5*5*5*5*8 = 5000 coeffs
#       family = poisson, method = "fREML", discrete = TRUE)
#
# TTPsplines stores O(d k r^2) coefficients instead of prod(k).
# Note: v0 uses a common basis size k per margin (bam's anisotropic k for
# time is not yet supported as a vector); k = 5 matches the pollutant margins.
#
# Requires Suggests: gamair
#   Rscript inst/examples/example_chicago_tt.R

root <- (function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), "../..")))
  }
  normalizePath(".")
})()

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required (or install TTPsplines and skip load_all).", call. = FALSE)
}
if (!requireNamespace("gamair", quietly = TRUE)) {
  stop("Install gamair for the chicago data: install.packages(\"gamair\")",
       call. = FALSE)
}

devtools::load_all(root, quiet = TRUE)
library(gamair)
data(chicago)

d <- stats::na.omit(chicago[, c(
  "death", "tmpd", "o3median", "pm10median", "so2median", "time"
)])
y <- d$death
X <- as.matrix(d[, c("tmpd", "o3median", "pm10median", "so2median", "time")])

# bam-style resolution on pollutants; isotropic k (package constraint)
k <- 5L
rank <- 3L
cx <- tt_complexity(d = ncol(X), p = k, rank = rank)

cat("Chicago deaths (gamair) — TT Poisson vs full te() storage\n")
cat(sprintf("  n = %d, d = %d\n", length(y), ncol(X)))
cat(sprintf("  Full te coeffs (k^d isotropic):     %s\n",
            format(cx$n_full, big.mark = ",")))
cat(sprintf("  bam k=c(5,5,5,5,8) product:         5,000\n"))
cat(sprintf("  TT stored (rank=%d, k=%d):          %s  (CR %.0fx)\n",
            rank, k, format(cx$n_tt_stored, big.mark = ","),
            cx$compression_ratio))
cat(sprintf("  TT intrinsic (minus gauge):         %s\n",
            format(cx$n_tt_intrinsic, big.mark = ",")))
cat("\nFitting ttps(..., family = poisson) — PIRLS-ALS, fixed lambda ...\n")

t0 <- proc.time()[[3L]]
fit <- ttps(
  y, X,
  family = poisson(),
  rank = rank,
  k = k,
  lambda = 1,                 # fixed; use "cGCV" later if desired
  optimizer = "auto",         # -> PIRLS-ALS
  control = tt_control(
    pirls_maxit = 20L,
    als_sweeps_per_pirls = 3L,
    max_sweeps = 12L,
    compute_edf = TRUE,
    edf_max_npar = 2500L,
    seed = 1L,
    trace = TRUE
  )
)
elapsed <- proc.time()[[3L]] - t0

print(summary(fit))
cat(sprintf("\nWall time: %.2fs\n", elapsed))
cat("Mean fitted deaths:", mean(fitted(fit)),
    "| mean observed:", mean(y), "\n")
cat("predict response [1:5]:\n")
print(predict(fit, X[1:5, , drop = FALSE], type = "response"))

# Optional: tiny additive bam for sanity (not the 5-way te)
if (requireNamespace("mgcv", quietly = TRUE)) {
  cat("\n--- Optional: mgcv additive bam (not the exploding te) ---\n")
  m_add <- mgcv::bam(
    death ~ s(tmpd, bs = "ps", k = 5) + s(o3median, bs = "ps", k = 5) +
      s(pm10median, bs = "ps", k = 5) + s(so2median, bs = "ps", k = 5) +
      s(time, bs = "ps", k = 8),
    family = poisson(), data = d, method = "fREML", discrete = TRUE
  )
  cat(sprintf("Additive bam deviance: %.4g | TT deviance: %.4g\n",
              m_add$deviance, deviance(fit)))
  cat("(Additive bam is a different model; shown only as a cheap reference.)\n")
}

cat("\nDone. To attempt the full te() (often very slow / memory-heavy):\n")
cat("  # bam(death ~ te(tmpd,o3median,pm10median,so2median,time,\n")
cat("  #               bs='ps', k=c(5,5,5,5,8)), family=poisson, ...)\n")
