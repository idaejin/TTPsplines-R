# Experimental global TT-gGCV laboratory (Monte Carlo GDF)

> **CLOSED — KEEP cGCV** (product). Experimental joint optimizer: `tt_global_lambda_optimize()`.
>
> **v0 GATE PASS** (\(d=2\), \(d=3\)) — CLOSED.  
> Product KEEP cGCV · TT-gGCV = internal oracle · v1 = IDEA parked (not pending).  
> Three-way separation locked. See `GATE_REPORT.md`.

Status: **product path closed (KEEP cGCV)**; **experimental joint λ optimizer** under construction (SAA + Sobol + nlminb + final multistart). Not classical exact GCV. Does not modify `lambda = "cGCV"` defaults.

## Matched-budget comparison (NEXT algorithmic gate)

Sobol vs grid / random / LHS at matched \(N\in\{16d,32d,64d\}\) — optimization of \(Q_r\), not prediction vs `cGCV`.

Protocol: [`PROTOCOL_MATCHED_BUDGET.md`](PROTOCOL_MATCHED_BUDGET.md)  
v0/v1 status: [`GATE_REPORT.md`](GATE_REPORT.md)

```bash
TT_GGCV_MB_MODE=SMOKE Rscript inst/benchmarks/global_gcv/run_matched_budget_gate.R
TT_GGCV_MB_MODE=GATE  Rscript inst/benchmarks/global_gcv/run_matched_budget_gate.R
```

Outputs: `results/matched_budget_{smoke,gate}/`

## Phase 0

```bash
Rscript inst/benchmarks/global_gcv/run_full_tensor_validation.R
Rscript inst/benchmarks/global_gcv/run_epsilon_M_sensitivity.R
Rscript inst/benchmarks/global_gcv/run_lambda_surface_d2.R
```

Report: `REPORT.md`

## Phase 1 QUICK

```bash
TT_GGCV_QUICK=true Rscript inst/benchmarks/global_gcv/phase1_run_all.R
```

Report: `PHASE1_REPORT.md` · outputs: `results/phase1/`

## Phase 1B FULL selective

```bash
# Smoke (separate directory)
TT_GGCV_MODE=SMOKE TT_GGCV_R=2 \
TT_GGCV_SCENARIOS=smooth_smooth,strong_aniso \
Rscript inst/benchmarks/global_gcv/phase1b_run_full.R

# Full
TT_GGCV_MODE=FULL TT_GGCV_R=50 \
TT_GGCV_SCENARIOS=smooth_smooth,strong_aniso \
Rscript inst/benchmarks/global_gcv/phase1b_run_full.R
```

Report: `PHASE1B_FULL_REPORT.md` · outputs: `results/phase1b_full/` (never writes into `results/phase1/`)

## Internal API

```r
TTPsplines:::tt_global_gdf_mc(...)
TTPsplines:::tt_global_gcv(...)
TTPsplines:::.tt_lab_eval_ggcv(...)
```
