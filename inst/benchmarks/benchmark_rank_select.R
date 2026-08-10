# Simulation audit: oracle vs min-CV vs 1-SE ranks
#   Rscript inst/benchmarks/benchmark_rank_select.R
#
# Oracle uses truth ONLY for evaluation — never inside tt_rank_select().

root <- if (file.exists("DESCRIPTION") &&
             identical(unname(read.dcf("DESCRIPTION")[, "Package"]), "TTPsplines")) {
  normalizePath(".")
} else {
  NULL
}
if (!is.null(root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, quiet = TRUE)
} else {
  library(TTPsplines)
}

audit_one <- function(name, y, X, f_true, ranks = 1:5, k = 6,
                      folds = 5, seed = 1, max_sweeps = 8) {
  ctrl <- tt_control(max_sweeps = max_sweeps, compute_edf = FALSE)
  t0 <- proc.time()[["elapsed"]]
  sel <- tt_rank_select(
    y, X, ranks = ranks, k = k, lambda = 1, folds = folds,
    rule = "1se", seed = seed, control = ctrl
  )
  cv_time <- proc.time()[["elapsed"]] - t0

  oracle_rmse <- vapply(ranks, function(r) {
    fit <- ttps(y, X, rank = r, k = k, lambda = 1, control = ctrl)
    sqrt(mean((fitted(fit) - f_true)^2))
  }, numeric(1))
  r_oracle <- ranks[which.min(oracle_rmse)]

  fit_sel <- tt_rank_refit(sel)
  fit_min <- tt_rank_refit(sel, rank = sel$rank_min)
  fit_ora <- ttps(y, X, rank = r_oracle, k = k, lambda = 1, control = ctrl)

  rmse <- function(fit) sqrt(mean((fitted(fit) - f_true)^2))
  data.frame(
    surface = name,
    oracle_rank = r_oracle,
    min_cv_rank = sel$rank_min,
    se1_rank = sel$rank_1se,
    oracle_rmse = rmse(fit_ora),
    min_cv_rmse = rmse(fit_min),
    se1_rmse = rmse(fit_sel),
    se1_npar = fit_sel$npar_tt,
    se1_CR = fit_sel$compression_ratio,
    cv_time_s = cv_time,
    stringsAsFactors = FALSE
  )
}

data(ishigami)
data(sobol_g)
data(friedman)

tab <- rbind(
  audit_one("Ishigami", ishigami$y,
            as.matrix(ishigami[, c("x1", "x2", "x3")]), ishigami$f, k = 8),
  audit_one("Sobol-g", sobol_g$y,
            as.matrix(sobol_g[, paste0("x", 1:4)]), sobol_g$f, k = 6),
  audit_one("Friedman", friedman$y,
            as.matrix(friedman[, paste0("x", 1:5)]), friedman$f, k = 5)
)
print(tab, digits = 4, row.names = FALSE)
