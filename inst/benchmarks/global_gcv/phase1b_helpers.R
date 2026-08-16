# Phase 1B replication helpers (cold_common primary; CRN; resume).

`%||%` <- function(a, b) if (!is.null(a)) a else b

.phase1b_append_csv <- function(df, path) {
  if (is.null(df) || !nrow(df)) return(invisible(path))
  if (!file.exists(path)) {
    utils::write.csv(df, path, row.names = FALSE)
  } else {
    utils::write.table(
      df, path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
    )
  }
  invisible(path)
}

.phase1b_done_keys <- function(path, key_cols = c("scenario", "rep")) {
  if (!file.exists(path)) return(character(0))
  tab <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || !nrow(tab)) return(character(0))
  if (!all(key_cols %in% names(tab))) return(character(0))
  apply(tab[, key_cols, drop = FALSE], 1L, function(r) paste(r, collapse = "|"))
}

.phase1b_paired_stats <- function(delta, label = "delta") {
  delta <- as.numeric(delta)
  delta <- delta[is.finite(delta)]
  n <- length(delta)
  if (!n) {
    return(data.frame(
      label = label, n = 0L, mean = NA_real_, median = NA_real_, sd = NA_real_,
      q10 = NA_real_, q90 = NA_real_,
      mean_lo = NA_real_, mean_hi = NA_real_,
      prop_pos = NA_real_, prop_neg = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  # Normal approx MC interval for the mean
  se <- stats::sd(delta) / sqrt(n)
  data.frame(
    label = label,
    n = n,
    mean = mean(delta),
    median = stats::median(delta),
    sd = stats::sd(delta),
    q10 = as.numeric(stats::quantile(delta, 0.10)),
    q90 = as.numeric(stats::quantile(delta, 0.90)),
    mean_lo = mean(delta) - 1.96 * se,
    mean_hi = mean(delta) + 1.96 * se,
    prop_pos = mean(delta > 0),
    prop_neg = mean(delta < 0),
    stringsAsFactors = FALSE
  )
}

.phase1b_run_cgcv_starts <- function(des, rank, common_init, ctrl_base,
                                     starts, probes, M_final, cache, cfg) {
  rows <- list()
  for (nm in names(starts)) {
    lam0 <- starts[[nm]]
    ctrl <- ctrl_base
    ctrl$lambda_start <- lam0
    ctrl$cgcv_update <- "outer_simultaneous"
    ctrl$cgcv_trace <- TRUE
    ctrl$outer_maxit <- if (cfg$smoke) 6L else 12L
    ctrl$warn_lambda_boundary <- FALSE
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      TTPsplines:::.tt_with_preserved_seed({
        ttps(
          des$y, des$X, family = gaussian(), rank = rank, k = des$k,
          degree = des$degree, penalty_order = des$penalty_order,
          lambda = "cGCV", optimizer = "ALS", backend = "R",
          init = TTPsplines:::.tt_clone_cores(common_init),
          control = ctrl
        )
      }),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    if (inherits(fit, "error")) {
      rows[[nm]] <- list(
        start = nm, ok = FALSE, reason = conditionMessage(fit),
        lambda = c(NA_real_, NA_real_), elapsed = elapsed, fit = NULL, ev = NULL
      )
      next
    }
    ev <- TTPsplines:::.tt_lab_eval_ggcv(
      lambda = as.numeric(fit$lambda), design = des, rank = rank,
      probes = probes, M = M_final,
      init_cores = common_init, init_policy = "cold_common",
      mc_bank_id = "bank0", cache = cache, control = ctrl_base,
      epsilon_rel = cfg$epsilon_rel
    )
    mu_te <- tryCatch(
      as.numeric(predict(fit, newdata = des$X_test, type = "response")),
      error = function(e) rep(NA_real_, des$n_test)
    )
    rows[[nm]] <- list(
      start = nm, ok = TRUE, reason = NA_character_,
      lambda = as.numeric(fit$lambda),
      log10_l1 = log10(fit$lambda[1]), log10_l2 = log10(fit$lambda[2]),
      log_ratio = log10(fit$lambda[1] / fit$lambda[2]),
      global_gcv = ev$global_gcv, gdf = ev$gdf, gdf_se = ev$gdf_mc_se,
      rss = ev$rss,
      test_mse = mean((mu_te - des$y_test)^2),
      ise = mean((mu_te - des$f_test)^2),
      converged = isTRUE(fit$converged),
      n_sweeps = fit$n_sweeps,
      at_boundary = isTRUE(fit$lambda_at_boundary),
      elapsed = elapsed,
      fit = fit, ev = ev
    )
  }
  rows
}

.phase1b_grid_oracle <- function(des, rank, common_init, probes, cfg, cache,
                                 ctrl, center) {
  box <- c(-5, 5)
  gr <- .phase1b_refine_grid(center, cfg$refine_halfwidth, cfg$refine_by, box)
  # Evaluate rectangular grid (may be anisotropic extents)
  coords <- expand.grid(log10_l1 = gr$g1, log10_l2 = gr$g2, KEEP.OUT.ATTRS = FALSE)
  rows <- vector("list", nrow(coords))
  for (i in seq_len(nrow(coords))) {
    lam <- 10^c(coords$log10_l1[i], coords$log10_l2[i])
    t0 <- proc.time()[["elapsed"]]
    ev <- TTPsplines:::.tt_lab_eval_ggcv(
      lambda = lam, design = des, rank = rank, probes = probes,
      M = cfg$M_search, init_cores = common_init, init_policy = "cold_common",
      mc_bank_id = "bank0", cache = cache, control = ctrl,
      epsilon_rel = cfg$epsilon_rel
    )
    full <- tryCatch(
      TTPsplines:::.tt_lab_eval_full_gcv(des, lam),
      error = function(e) list(gcv = NA_real_, edf = NA_real_, rss = NA_real_)
    )
    rows[[i]] <- data.frame(
      log10_l1 = coords$log10_l1[i], log10_l2 = coords$log10_l2[i],
      lambda1 = lam[1], lambda2 = lam[2],
      global_gcv = ev$global_gcv, gdf = ev$gdf, rss = ev$rss,
      valid = ev$valid, converged = ev$converged, n_sweeps = ev$n_sweeps,
      test_mse = ev$rmse_test^2, ise = ev$ise_true,
      gcv_full = full$gcv, edf_full = full$edf,
      elapsed = proc.time()[["elapsed"]] - t0,
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, rows)
  tab$flat_1pct <- TTPsplines:::.tt_lab_flatness_mask(tab$global_gcv, 0.01)

  ok <- is.finite(tab$global_gcv) & isTRUE(tab$valid %in% TRUE)
  if (!any(ok)) ok <- is.finite(tab$global_gcv)
  i_best <- which(ok)[which.min(tab$global_gcv[ok])]
  on_edge <- tab$log10_l1[i_best] %in% range(gr$g1) ||
    tab$log10_l2[i_best] %in% range(gr$g2)

  widened <- FALSE
  if (isTRUE(on_edge) && isTRUE(cfg$widen_once)) {
    widened <- TRUE
    gr2 <- .phase1b_refine_grid(
      center, cfg$refine_halfwidth + cfg$refine_by, cfg$refine_by, box
    )
    coords2 <- expand.grid(log10_l1 = gr2$g1, log10_l2 = gr2$g2, KEEP.OUT.ATTRS = FALSE)
    # Only evaluate new points
    key_old <- paste(tab$log10_l1, tab$log10_l2)
    key_new <- paste(coords2$log10_l1, coords2$log10_l2)
    add <- coords2[!key_new %in% key_old, , drop = FALSE]
    extra <- vector("list", nrow(add))
    for (i in seq_len(nrow(add))) {
      lam <- 10^c(add$log10_l1[i], add$log10_l2[i])
      ev <- TTPsplines:::.tt_lab_eval_ggcv(
        lambda = lam, design = des, rank = rank, probes = probes,
        M = cfg$M_search, init_cores = common_init, init_policy = "cold_common",
        mc_bank_id = "bank0", cache = cache, control = ctrl,
        epsilon_rel = cfg$epsilon_rel
      )
      full <- tryCatch(
        TTPsplines:::.tt_lab_eval_full_gcv(des, lam),
        error = function(e) list(gcv = NA_real_, edf = NA_real_)
      )
      extra[[i]] <- data.frame(
        log10_l1 = add$log10_l1[i], log10_l2 = add$log10_l2[i],
        lambda1 = lam[1], lambda2 = lam[2],
        global_gcv = ev$global_gcv, gdf = ev$gdf, rss = ev$rss,
        valid = ev$valid, converged = ev$converged, n_sweeps = ev$n_sweeps,
        test_mse = ev$rmse_test^2, ise = ev$ise_true,
        gcv_full = full$gcv, edf_full = full$edf,
        elapsed = NA_real_, stringsAsFactors = FALSE
      )
    }
    if (length(extra)) {
      add_tab <- do.call(rbind, extra)
      # Align columns with tab (extra may omit flat_1pct)
      for (nm in setdiff(names(tab), names(add_tab))) {
        add_tab[[nm]] <- NA
      }
      add_tab <- add_tab[, names(tab), drop = FALSE]
      tab <- rbind(tab, add_tab)
      tab$flat_1pct <- TTPsplines:::.tt_lab_flatness_mask(tab$global_gcv, 0.01)
      ok <- is.finite(tab$global_gcv)
      i_best <- which(ok)[which.min(tab$global_gcv[ok])]
      on_edge <- tab$log10_l1[i_best] %in% range(c(gr$g1, gr2$g1)) ||
        tab$log10_l2[i_best] %in% range(c(gr$g2, gr2$g2))
    }
  }

  # Final re-eval with M_final
  lam_best <- c(tab$lambda1[i_best], tab$lambda2[i_best])
  evf <- TTPsplines:::.tt_lab_eval_ggcv(
    lambda = lam_best, design = des, rank = rank, probes = probes,
    M = cfg$M_final, init_cores = common_init, init_policy = "cold_common",
    mc_bank_id = "bank0", cache = cache, control = ctrl,
    epsilon_rel = cfg$epsilon_rel
  )
  list(
    surface = tab,
    i_best = i_best,
    on_edge = on_edge,
    widened = widened,
    lambda = lam_best,
    log10 = c(tab$log10_l1[i_best], tab$log10_l2[i_best]),
    ev_final = evf,
    flat_n = sum(tab$flat_1pct, na.rm = TRUE),
    flat_frac = mean(tab$flat_1pct, na.rm = TRUE)
  )
}

#' One Phase-1B replicate for a single scenario.
.phase1b_one_rep <- function(scenario, rep, cfg) {
  seed <- as.integer(cfg$seed0 + (rep - 1L) * 97L +
                       match(scenario, cfg$scenarios) * 13L)
  t_rep0 <- proc.time()[["elapsed"]]
  des <- TTPsplines:::.tt_lab_phase1_make_design(
    scenario = scenario, n = cfg$n, n_test = cfg$n_test, k = cfg$k, seed = seed
  )
  probes <- .phase1_probe_bank(des$n, cfg$M_bank, seed + 99L)
  ranks <- .phase1_ranks(des$k)
  rk <- ranks$sufficient
  ctrl <- TTPsplines:::.tt_lab_phase1_control(seed = seed, max_sweeps = cfg$max_sweeps)
  common_init <- TTPsplines:::.tt_with_preserved_seed({
    tt_initialize(d = 2L, rank = rk, k = des$k, seed = seed)
  })
  cache <- TTPsplines:::.tt_lab_new_cache()

  starts <- list(
    low = 10^c(-3, -3),
    central = c(1, 1),
    high = 10^c(3, 3)
  )

  # --- cGCV starts ---
  cg <- .phase1b_run_cgcv_starts(
    des, rk, common_init, ctrl, starts, probes, cfg$M_final, cache, cfg
  )
  cg_ok <- cg[vapply(cg, function(z) isTRUE(z$ok), logical(1))]
  # cGCV_default = central (documented default-like)
  cg_default <- cg[["central"]]
  # cGCV_best_global = min TT-gGCV among starts (no test set)
  if (length(cg_ok)) {
    gcv_vals <- vapply(cg_ok, function(z) as.numeric(z$global_gcv), numeric(1))
    cg_best <- cg_ok[[which.min(gcv_vals)]]
  } else {
    cg_best <- cg_default
  }

  # --- Grid oracle ---
  center <- cfg$centers[[scenario]]
  if (is.null(center)) center <- c(0, 0)
  oracle <- .phase1b_grid_oracle(
    des, rk, common_init, probes, cfg, cache, ctrl, center
  )

  # --- Full-TP min on same refined surface (geometry check) ---
  surf <- oracle$surface
  i_full <- which.min(surf$gcv_full)
  full_th <- c(surf$log10_l1[i_full], surf$log10_l2[i_full])

  # --- Optional Nelder ---
  nelder <- NULL
  if (isTRUE(cfg$do_nelder) && rep <= cfg$nelder_max_reps) {
    starts_nm <- list(
      center = c(0, 0),
      cgcv = c(cg_default$log10_l1 %||% 0, cg_default$log10_l2 %||% 0),
      grid = oracle$log10,
      low = c(-3, -3),
      high = c(3, 3)
    )
    nelder <- tryCatch(
      TTPsplines:::.tt_lab_optimize_ggcv(
        design = des, rank = rk, probes = probes,
        M_search = cfg$M_search, M_final = cfg$M_final,
        starts = starts_nm, box = c(-5, 5), common_init = common_init,
        mc_bank_id = "bank0", cache = cache, control = ctrl,
        epsilon_rel = cfg$epsilon_rel,
        optim_control = list(maxit = if (cfg$smoke) 15L else 40L, reltol = 1e-3)
      ),
      error = function(e) NULL
    )
  }

  # --- Optional restricted rank (secondary) ---
  restr <- NULL
  if (isTRUE(cfg$do_restricted) && rep <= cfg$restricted_max_reps) {
    rk_r <- ranks$restricted
    init_r <- TTPsplines:::.tt_with_preserved_seed({
      tt_initialize(d = 2L, rank = rk_r, k = des$k, seed = seed)
    })
    restr <- tryCatch(
      .phase1b_grid_oracle(
        des, rk_r, init_r, probes, cfg, TTPsplines:::.tt_lab_new_cache(),
        ctrl, center
      ),
      error = function(e) list(error = conditionMessage(e))
    )
  }

  # Build method rows (primary methods)
  method_row <- function(method, log10, lam, ev, test_mse, ise,
                         elapsed_sec = NA_real_, cgcv_start = NA_character_) {
    in_flat <- NA
    if (!is.null(oracle$surface) && all(is.finite(log10))) {
      d2 <- (oracle$surface$log10_l1 - log10[1])^2 +
        (oracle$surface$log10_l2 - log10[2])^2
      j <- which.min(d2)
      in_flat <- isTRUE(oracle$surface$flat_1pct[j])
    }
    data.frame(
      mode = cfg$mode,
      scenario = scenario,
      rep = as.integer(rep),
      seed = seed,
      method = method,
      rank = as.integer(rk)[1L],
      log10_l1 = as.numeric(log10[1]), log10_l2 = as.numeric(log10[2]),
      lambda1 = as.numeric(lam[1]), lambda2 = as.numeric(lam[2]),
      log_ratio = as.numeric(log10(lam[1] / lam[2])),
      global_gcv = as.numeric(ev$global_gcv %||% NA_real_)[1],
      gdf = as.numeric(ev$gdf %||% NA_real_)[1],
      gdf_se = as.numeric(ev$gdf_mc_se %||% NA_real_)[1],
      rss = as.numeric(ev$rss %||% NA_real_)[1],
      test_mse = as.numeric(test_mse)[1],
      ise = as.numeric(ise)[1],
      dist_to_oracle = as.numeric(sqrt(sum((log10 - oracle$log10)^2)))[1],
      in_flat_1pct = in_flat,
      oracle_on_edge = isTRUE(oracle$on_edge),
      oracle_widened = isTRUE(oracle$widened),
      flat_frac = as.numeric(oracle$flat_frac)[1],
      converged = isTRUE(ev$converged %||% TRUE),
      n_sweeps = as.integer(ev$n_sweeps %||% NA_integer_)[1],
      elapsed_sec = as.numeric(elapsed_sec)[1],
      cgcv_start = as.character(cgcv_start)[1],
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  rows[[1]] <- method_row(
    "tt_grid_oracle", oracle$log10, oracle$lambda, oracle$ev_final,
    oracle$ev_final$rmse_test^2, oracle$ev_final$ise_true,
    elapsed_sec = sum(oracle$surface$elapsed, na.rm = TRUE)
  )
  ev_full <- TTPsplines:::.tt_lab_eval_ggcv(
    lambda = 10^full_th, design = des, rank = rk, probes = probes,
    M = cfg$M_final, init_cores = common_init, init_policy = "cold_common",
    mc_bank_id = "bank0", cache = cache, control = ctrl,
    epsilon_rel = cfg$epsilon_rel
  )
  rows[[2]] <- method_row(
    "full_tp_gcv", full_th, 10^full_th, ev_full,
    ev_full$rmse_test^2, ev_full$ise_true
  )
  if (isTRUE(cg_default$ok)) {
    rows[[length(rows) + 1L]] <- method_row(
      "cGCV_default", c(cg_default$log10_l1, cg_default$log10_l2),
      cg_default$lambda, cg_default$ev,
      cg_default$test_mse, cg_default$ise,
      elapsed_sec = cg_default$elapsed, cgcv_start = "central"
    )
  }
  if (isTRUE(cg_best$ok)) {
    rows[[length(rows) + 1L]] <- method_row(
      "cGCV_best_global", c(cg_best$log10_l1, cg_best$log10_l2),
      cg_best$lambda, cg_best$ev,
      cg_best$test_mse, cg_best$ise,
      elapsed_sec = cg_best$elapsed, cgcv_start = cg_best$start
    )
  }
  for (nm in names(cg)) {
    z <- cg[[nm]]
    if (!isTRUE(z$ok)) next
    rows[[length(rows) + 1L]] <- method_row(
      paste0("cGCV_start_", nm), c(z$log10_l1, z$log10_l2), z$lambda, z$ev,
      z$test_mse, z$ise,
      elapsed_sec = z$elapsed, cgcv_start = nm
    )
  }
  if (!is.null(nelder) && nrow(nelder)) {
    okn <- nelder[isTRUE(nelder$ok) | nelder$ok == TRUE, , drop = FALSE]
    if (nrow(okn)) {
      j <- which.min(okn$global_gcv_final)
      evn <- TTPsplines:::.tt_lab_eval_ggcv(
        lambda = c(okn$lambda1[j], okn$lambda2[j]), design = des, rank = rk,
        probes = probes, M = cfg$M_final, init_cores = common_init,
        init_policy = "cold_common", mc_bank_id = "bank0", cache = cache,
        control = ctrl, epsilon_rel = cfg$epsilon_rel
      )
      rows[[length(rows) + 1L]] <- method_row(
        "tt_nelder", c(okn$log10_l1[j], okn$log10_l2[j]),
        c(okn$lambda1[j], okn$lambda2[j]), evn,
        evn$rmse_test^2, evn$ise_true,
        elapsed_sec = sum(nelder$elapsed_sec, na.rm = TRUE)
      )
    }
  }

  rep_tab <- do.call(rbind, rows)
  # Align columns
  # Stability table for this rep
  stab <- do.call(rbind, lapply(names(cg), function(nm) {
    z <- cg[[nm]]
    if (!isTRUE(z$ok)) {
      return(data.frame(
        mode = cfg$mode, scenario = scenario, rep = rep, seed = seed,
        start = nm, ok = FALSE, reason = z$reason %||% NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      mode = cfg$mode, scenario = scenario, rep = rep, seed = seed,
      start = nm, ok = TRUE, reason = NA_character_,
      log10_l1 = z$log10_l1, log10_l2 = z$log10_l2, log_ratio = z$log_ratio,
      global_gcv = z$global_gcv, gdf = z$gdf, gdf_se = z$gdf_se,
      test_mse = z$test_mse, ise = z$ise,
      delta_gcv_vs_oracle = z$global_gcv - oracle$ev_final$global_gcv,
      delta_mse_vs_oracle = z$test_mse - (oracle$ev_final$rmse_test^2),
      delta_ise_vs_oracle = z$ise - oracle$ev_final$ise_true,
      dist_to_oracle = sqrt((z$log10_l1 - oracle$log10[1])^2 +
                              (z$log10_l2 - oracle$log10[2])^2),
      at_boundary = z$at_boundary,
      elapsed_sec = z$elapsed,
      stringsAsFactors = FALSE
    )
  }))

  aniso <- data.frame(
    mode = cfg$mode, scenario = scenario, rep = rep, seed = seed,
    method = c("tt_grid_oracle", "cGCV_default", "cGCV_best_global", "full_tp_gcv"),
    log_ratio = c(
      log10(oracle$lambda[1] / oracle$lambda[2]),
      if (isTRUE(cg_default$ok)) cg_default$log_ratio else NA_real_,
      if (isTRUE(cg_best$ok)) cg_best$log_ratio else NA_real_,
      log10((10^full_th[1]) / (10^full_th[2]))
    ),
    sign_pos = c(
      log10(oracle$lambda[1] / oracle$lambda[2]) > 0,
      if (isTRUE(cg_default$ok)) cg_default$log_ratio > 0 else NA,
      if (isTRUE(cg_best$ok)) cg_best$log_ratio > 0 else NA,
      (full_th[1] - full_th[2]) > 0
    ),
    truth_sign = cfg$aniso_sign_truth[[scenario]] %||% 0,
    stringsAsFactors = FALSE
  )
  aniso$sign_match_truth <- if (scenario == "strong_aniso") {
    aniso$sign_pos == TRUE
  } else {
    NA
  }

  list(
    replication = rep_tab,
    stability = stab,
    anisotropy = aniso,
    surface = cbind(
      mode = cfg$mode, scenario = scenario, rep = rep, seed = seed,
      oracle$surface, stringsAsFactors = FALSE
    ),
    restricted = restr,
    elapsed_sec = proc.time()[["elapsed"]] - t_rep0,
    failure = NULL
  )
}
