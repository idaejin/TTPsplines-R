# GLAM / full-tensor vs TT (FIXED λ, same basis / data)
# Compression + surface recovery — not λ selection.
# source("inst/benchmarks/benchmark_glam.R")

local({
  bench_dir <- (function() {
    of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(of)) return(dirname(normalizePath(of)))
    if (file.exists("helpers.R")) return(normalizePath("."))
    if (file.exists("inst/benchmarks/helpers.R"))
      return(normalizePath("inst/benchmarks"))
    system.file("benchmarks", package = "TTPsplines")
  })()
  source(file.path(bench_dir, "helpers.R"), local = TRUE)
  .ttps_bench_ensure_pkg()
  out_dir <- .ttps_bench_outdir()

  message("=== TTPsplines benchmark: GLAM grid vs TT (fixed λ) ===")

  # Modest 3D grid so dense p^3 stays feasible
  n_grid <- c(20L, 18L, 16L)
  p <- 8L
  lambda <- c(1, 1, 1)
  degree <- 3L

  g <- lapply(n_grid, function(nk) seq(0, 1, length.out = nk))
  # Truth array
  truth <- array(0, n_grid)
  for (i in seq_along(g[[1]])) {
    for (j in seq_along(g[[2]])) {
      for (k in seq_along(g[[3]])) {
        xx <- cbind(g[[1]][i], g[[2]][j], g[[3]][k])
        truth[i, j, k] <- true_surface_nd(xx)
      }
    }
  }
  set.seed(51)
  Y <- truth + array(rnorm(length(truth), 0, 0.25), dim(truth))

  # Marginal bases on grid axes
  knots <- lapply(g, function(gi) {
    # use package internal via ::: if exported helpers unavailable
    TTPsplines:::make_knots(gi, k = p, degree = degree)
  })
  B <- lapply(seq_along(g), function(j) {
    TTPsplines:::bspline_basis(g[[j]], knots[[j]], degree = degree)
  })

  # --- GLAM ---
  t0 <- proc.time()[["elapsed"]]
  glam <- glam_fit_gaussian(Y, B, lambda = lambda)
  time_glam <- proc.time()[["elapsed"]] - t0

  # --- TT on vectorized grid (same knots) ---
  # Expand observations = all cells, bases = row of B_k at cell index
  idx <- expand.grid(lapply(n_grid, seq_len), KEEP.OUT.ATTRS = FALSE)
  X_sc <- cbind(g[[1]][idx[[1]]], g[[2]][idx[[2]]], g[[3]][idx[[3]]])
  y_sc <- as.numeric(Y)
  truth_sc <- as.numeric(truth)

  rows <- list()
  for (r in c(1L, 2L, 3L, 4L, 6L)) {
    ctrl <- .default_control(backend = "auto", max_sweeps = 15L, seed = 3L)
    t1 <- proc.time()[["elapsed"]]
    fit <- ttps(
      y_sc, X_sc,
      family = gaussian(),
      rank = r,
      k = p,
      lambda = lambda,
      knots = knots,
      control = ctrl
    )
    elapsed <- proc.time()[["elapsed"]] - t1
    mu_tt <- fitted(fit)
    rows[[length(rows) + 1L]] <- data.frame(
      method = sprintf("TT-r%d", r),
      rank = r,
      npar = fit$npar_tt,
      npar_dense = fit$npar_dense,
      compression = fit$compression_ratio,
      rmse_y = rmse(mu_tt, y_sc),
      rmse_truth = rmse(mu_tt, truth_sc),
      rmse_vs_glam = rmse(mu_tt, as.numeric(glam$mu)),
      time_s = elapsed,
      backend = fit$backend,
      stringsAsFactors = FALSE
    )
    message(sprintf(
      "  TT r=%d | RMSE_truth=%.4f | vs GLAM=%.4f | npar=%d | %.2fs",
      r, rmse(mu_tt, truth_sc), rmse(mu_tt, as.numeric(glam$mu)),
      fit$npar_tt, elapsed
    ))
  }

  glam_row <- data.frame(
    method = "GLAM-fixed",
    rank = NA_integer_,
    npar = glam$npar,
    npar_dense = glam$npar,
    compression = 1,
    rmse_y = rmse(glam$mu, Y),
    rmse_truth = rmse(glam$mu, truth),
    rmse_vs_glam = 0,
    time_s = time_glam,
    backend = "array",
    stringsAsFactors = FALSE
  )

  tab <- rbind(glam_row, do.call(rbind, rows))
  print(tab, row.names = FALSE)
  .write_bench_csv(tab, "benchmark_glam.csv", out_dir)

  fig <- file.path(out_dir, "benchmark_glam.png")
  png(fig, width = 900, height = 400)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  tt <- tab[!is.na(tab$rank), ]
  plot(tt$rank, tt$rmse_truth, type = "b", pch = 19,
       xlab = "TT rank", ylab = "RMSE vs truth",
       main = "TT vs truth (fixed λ)")
  abline(h = glam_row$rmse_truth, col = 2, lty = 2)
  legend("topright", c("TT", "GLAM"), col = c(1, 2), lty = c(1, 2), bty = "n")
  plot(tt$rank, tt$npar, type = "b", pch = 19,
       xlab = "TT rank", ylab = "parameters",
       main = "Compression", ylim = range(c(tt$npar, glam_row$npar)))
  abline(h = glam_row$npar, col = 2, lty = 2)
  par(op)
  dev.off()
  message("Wrote ", fig)
  message(sprintf("GLAM time=%.3fs | cells=%d | p^3=%d",
                  time_glam, prod(n_grid), p^3))
})
