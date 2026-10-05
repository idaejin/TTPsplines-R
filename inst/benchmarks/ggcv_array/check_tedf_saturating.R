## check_tedf_saturating.R — package T-EDF (fit$edf) versus the exact EDF at
## saturating TT rank, where the TT fit equals the dense tensor P-spline (GLAM).
## Exact EDF = tr((B'B + P)^{-1} B'B) with B = B2 (x) B1 and P = glam_penalty().
## Found 2026-10-01 while validating tt_ggcv_array(): T-EDF exceeds the exact
## EDF by about 9-10 df at moderate/strong lambda (d = 2, K = r = 5); same
## value at commit 34f06cf (2026-08-17), so not a regression.
## Run from the package root: Rscript inst/benchmarks/ggcv_array/check_tedf_saturating.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
ng <- c(8L, 7L); ax <- lapply(ng, function(m) seq(0, 1, length.out = m)); names(ax) <- c("x1", "x2")
K <- 5L; bb <- glam_grid_bases(ax, k = K); B <- kronecker(bb$B[[2]], bb$B[[1]])
exact_edf <- function(l) { P <- glam_penalty(c(K, K), l); sum(diag(solve(crossprod(B) + P, crossprod(B)))) }
set.seed(1); Y <- array(rnorm(prod(ng)), ng)
ce <- tt_control(max_sweeps = 30L, seed = 1L, compute_edf = TRUE, tol = 0)
rows <- lapply(list(c(1e-3, 1e-3), c(0.3, 1), c(1, 1), c(10, 10), c(1e3, 1e3)), function(lam) {
  f <- ttps(Y, axes = ax, array = TRUE, rank = K, k = K, lambda = lam, optimizer = "ALS", control = ce)
  u <- tt_gdf_array(Y, lambda = lam, axes = ax, rank = K, k = K, probes = "unit",
                    budget = "fixed", control = ce)  # the fit's own 30-sweep budget
  data.frame(lambda1 = lam[1], lambda2 = lam[2], edf_exact = exact_edf(lam),
             gdf_unit = u$gdf, tedf = as.numeric(f$edf)[1])
})
res <- do.call(rbind, rows)
print(res, digits = 4)
utils::write.csv(res, file.path("inst", "benchmarks", "ggcv_array", "check_tedf_saturating.csv"), row.names = FALSE)
