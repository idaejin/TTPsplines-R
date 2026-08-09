# TTPsplines benchmarks
#
# Reproducible examples comparing rank, λ mode (fixed / cGCV), and GLAM.
# **Not** executed by `R CMD check`.
#
# ## Quick start (development)
#
# ```r
# setwd("path/to/ttpsplines-pkg")
# devtools::load_all()
# source("inst/benchmarks/run_all.R")
# ```
#
# Or selectively:
#
# ```bash
# TTPSPLINES_BENCH_WHICH=gaussian,glam Rscript inst/benchmarks/run_all.R
# ```
#
# Outputs (CSV + PNG) go to `inst/benchmarks/results/` by default
# (override with `TTPSPLINES_BENCH_OUT`).
#
# | Script | What |
# |---|---|
# | `benchmark_gaussian.R` | Scattered Gaussian; ranks × {fixed, cGCV} |
# | `benchmark_poisson.R` | Scattered Poisson; ranks × {fixed, cGCV} |
# | `benchmark_bernoulli.R` | Scattered Bernoulli; ranks × {fixed, cGCV} |
# | `benchmark_rank.R` | `tt_rank_profile` + test RMSE |
# | `benchmark_glam.R` | Grid GLAM vs TT, **same fixed λ** (compression) |
# | `benchmark_optimizers_fixed.R` | ALS vs LBFGS vs Adam stub; fixed λ; same init |
# | `benchmark_optimizers_cgcv.R` | ALS vs LBFGS vs Adam stub; cGCV; same init |
# | `diagnose_als_vs_lbfgs.R` | Objective / multi-init / cross warm-start / λ paths |
# | `benchmark_bernoulli_audit.R` | Bernoulli PIRLS vs LBFGS root-cause audit (`--quick`/`--full`) |
# | `python/benchmark_keras_tt.py` | Optional autodiff comparator (no TF dependency) |
#
# Shared DGP / metrics: `helpers.R`.
#
# Optimizer fair protocol: same data, bases, knots, rank, λ (or same cGCV),
# same `tt_initialize()` cores. Compare η/μ/deviance/time — not TT cores.
