# Phase 1 — cGCV basins: starts x margin orders; re-score with TT-gGCV.
# Usage: TT_GGCV_QUICK=true Rscript inst/benchmarks/global_gcv/phase1_basins.R

source(file.path("inst", "benchmarks", "global_gcv", "phase1_config.R"))
.phase1_load_pkg()
cfg <- .phase1_config()
.phase1_ensure_dirs(cfg)

cache <- TTPsplines:::.tt_lab_new_cache()
lambda_starts <- if (cfg$quick) {
  list(low = c(1e-3, 1e-3), mid = c(1, 1), high = c(1e3, 1e3))
} else {
  list(
    low = c(1e-4, 1e-4), mid = c(1, 1), high = c(1e4, 1e4),
    aniso_lh = c(1e-2, 1e2), aniso_hl = c(1e2, 1e-2)
  )
}
orders <- list(
  `12` = c(1L, 2L),
  `21` = c(2L, 1L)
)

rows <- list()
traj_rows <- list()
ri <- 0L
ti <- 0L

for (sc in cfg$scenarios) {
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    scenario = sc, n = cfg$n, n_test = cfg$n_test, k = cfg$k, seed = cfg$seed0
  )
  probes <- .phase1_probe_bank(des$n, cfg$M_bank, cfg$seed0 + 99L)
  rk <- .phase1_ranks(des$k)$sufficient
  common_init <- TTPsplines:::.tt_with_preserved_seed({
    tt_initialize(d = 2L, rank = rk, k = des$k, seed = des$seed)
  })

  for (st_name in names(lambda_starts)) {
    for (ord_name in names(orders)) {
      lam0 <- lambda_starts[[st_name]]
      mord <- orders[[ord_name]]
      ctrl <- tt_control(
        max_sweeps = cfg$max_sweeps,
        outer_maxit = if (cfg$quick) 6L else 12L,
        compute_edf = FALSE,
        seed = des$seed,
        backend = "R",
        cgcv_update = "outer_simultaneous",
        cgcv_trace = TRUE,
        cgcv_margin_order = mord,
        lambda_start = lam0,
        warn_lambda_boundary = FALSE
      )
      t0 <- proc.time()[["elapsed"]]
      init_cgcv <- TTPsplines:::.tt_clone_cores(common_init)
      fit <- TTPsplines:::.tt_with_preserved_seed({
        ttps(
          des$y, des$X, family = gaussian(), rank = rk, k = des$k,
          degree = des$degree, penalty_order = des$penalty_order,
          lambda = "cGCV", optimizer = "ALS", backend = "R",
          init = init_cgcv,
          control = ctrl
        )
      })
      elapsed <- proc.time()[["elapsed"]] - t0

      # Re-score with global TT-gGCV (fixed lambda = cGCV solution)
      ev <- TTPsplines:::.tt_lab_eval_ggcv(
        lambda = as.numeric(fit$lambda), design = des, rank = rk,
        probes = probes, M = cfg$M_final,
        init_cores = common_init, init_policy = "cold_common",
        mc_bank_id = "bank0", cache = cache,
        control = TTPsplines:::.tt_lab_phase1_control(
          seed = des$seed, max_sweeps = cfg$max_sweeps
        ),
        epsilon_rel = cfg$epsilon_rel
      )
      mu_te <- tryCatch(
        as.numeric(predict(fit, newdata = des$X_test, type = "response")),
        error = function(e) rep(NA_real_, des$n_test)
      )
      rmse_test <- sqrt(mean((mu_te - des$y_test)^2))
      ise_true <- mean((mu_te - des$f_test)^2)

      ri <- ri + 1L
      rows[[ri]] <- data.frame(
        scenario = sc,
        start = st_name,
        margin_order = ord_name,
        lambda1 = fit$lambda[1],
        lambda2 = fit$lambda[2],
        log10_l1 = log10(fit$lambda[1]),
        log10_l2 = log10(fit$lambda[2]),
        lambda_boundary = paste(fit$lambda_boundary, collapse = ","),
        at_boundary = isTRUE(fit$lambda_at_boundary),
        n_outer = fit$n_outer %||% NA_integer_,
        n_sweeps = fit$n_sweeps,
        converged = fit$converged,
        global_gcv = ev$global_gcv,
        gdf = ev$gdf,
        gdf_se = ev$gdf_mc_se,
        rss = ev$rss,
        rmse_test = rmse_test,
        ise_true = ise_true,
        elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      )

      # Trajectory from cgcv trace if present
      tr <- fit$cgcv$trace
      if (is.data.frame(tr) && nrow(tr)) {
        # best-effort: record available numeric columns
        keep <- intersect(names(tr), c(
          "sweep", "outer", "margin", "lambda", "lambda_new",
          "gcv", "ed", "k"
        ))
        if (length(keep)) {
          tmp <- tr[, keep, drop = FALSE]
          tmp$scenario <- sc
          tmp$start <- st_name
          tmp$margin_order <- ord_name
          ti <- ti + 1L
          traj_rows[[ti]] <- tmp
        }
      }
      # Also store proposals if list of lambdas in history
      if (!is.null(fit$history) && length(fit$history)) {
        for (h in fit$history) {
          if (!is.null(h$lambda)) {
            ti <- ti + 1L
            traj_rows[[ti]] <- data.frame(
              scenario = sc, start = st_name, margin_order = ord_name,
              sweep = h$sweep %||% NA_integer_,
              lambda1 = h$lambda[1], lambda2 = h$lambda[2],
              stringsAsFactors = FALSE
            )
          }
        }
      }

      message(sprintf(
        "[%s|%s|%s] lam=(%.3g,%.3g) gGCV=%.4g rmse=%.3g",
        sc, st_name, ord_name, fit$lambda[1], fit$lambda[2],
        ev$global_gcv, rmse_test
      ))
    }
  }
}

tab <- do.call(rbind, rows)
.phase1_checkpoint_write(tab, file.path(cfg$out_dir, "phase1_basins.csv"))

if (length(traj_rows)) {
  # rbind may fail if schemas differ — bind carefully
  traj <- tryCatch(do.call(rbind, traj_rows), error = function(e) {
    # fallback: write separately as list summary
    NULL
  })
  if (!is.null(traj)) {
    .phase1_checkpoint_write(traj, file.path(cfg$out_dir, "phase1_cgcv_trajectories.csv"))
  }
}

# Trajectory / basin scatter figure
png(file.path(cfg$fig_dir, "phase1_cgcv_basins.png"), width = 900, height = 700)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (sc in cfg$scenarios) {
  sub <- tab[tab$scenario == sc, ]
  plot(sub$log10_l1, sub$log10_l2, pch = 16,
       col = ifelse(sub$margin_order == "12", "black", "darkred"),
       xlab = "log10(l1)", ylab = "log10(l2)",
       main = paste("cGCV basins:", sc),
       xlim = c(cfg$coarse_lo, cfg$coarse_hi),
       ylim = c(cfg$coarse_lo, cfg$coarse_hi))
  text(sub$log10_l1, sub$log10_l2, labels = sub$start, pos = 3, cex = 0.7)
  legend("topright", c("order 1-2", "order 2-1"),
         col = c("black", "darkred"), pch = 16, bty = "n", cex = 0.8)
}
par(op)
dev.off()

print(tab)
message("Basins done.")
