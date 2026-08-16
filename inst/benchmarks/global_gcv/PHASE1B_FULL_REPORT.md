# Phase 1B FULL report: selective validation of TT-gGCV vs cGCV

> **CLOSED — KEEP cGCV** (2026-08-16). Track closed. Do not run further Monte Carlo of this type. Do not put TT-gGCV in the current manuscript.
>
> Canonical conclusion: *cGCV remains the operational smoothing-parameter selector. Although the joint TT-gGCV oracle can attain lower prediction error in a strongly anisotropic scenario, this advantage is associated with flat objective regions and boundary minima and does not translate into better recovery of the underlying anisotropic smoothing structure. Multistart initialization provides no material improvement. TT-gGCV is therefore retained only as an internal diagnostic benchmark.*

**Status:** internal laboratory. Not a public method.  
**Date:** 2026-08-16  
**Mode:** `FULL` (not SMOKE / not QUICK)  
**Artifacts:** `results/phase1b_full/` · `figures/phase1b_full/`  
**QUICK Phase 1 preserved at:** `results/phase1/` (untouched)

---

## Configuration (exact)

| Item | Value |
|------|------:|
| Scenarios | `smooth_smooth`, `strong_aniso` |
| Replicates \(R\) | 50 per scenario (100 completed keys) |
| \(n\) / \(n_{\mathrm{test}}\) | 200 / 400 |
| \(k\) | 8 |
| Rank (primary) | sufficient \(=k\) |
| \(M_{\mathrm{search}}\) / \(M_{\mathrm{final}}\) | 5 / 20 |
| MC bank | 40 (CRN within replicate) |
| Init policy | `cold_common` |
| Refined grid | half-width 1.5, step 0.5 about QUICK centers \((-1,-3)\) and \((3,3)\); widen-once if edge |
| cGCV starts | low / central / high (`cGCV_default` = central) |
| Nelder / restricted | first 10 reps only |
| Wall time | **1037 s (~17.3 min)** |
| Failures | **0** |

Reproduce:

```bash
TT_GGCV_MODE=FULL TT_GGCV_R=50 \
TT_GGCV_SCENARIOS=smooth_smooth,strong_aniso \
TT_GGCV_DO_NELDER=true TT_GGCV_DO_RESTRICTED=true \
Rscript inst/benchmarks/global_gcv/phase1b_run_full.R
```

Smoke (already run; separate dir):

```bash
TT_GGCV_MODE=SMOKE TT_GGCV_R=2 \
TT_GGCV_SCENARIOS=smooth_smooth,strong_aniso \
Rscript inst/benchmarks/global_gcv/phase1b_run_full.R
```

---

## FULL vs QUICK

| Aspect | QUICK (Phase 1) | FULL (Phase 1B) |
|--------|-----------------|-----------------|
| Scenarios | 4 | 2 interior candidates |
| \(R\) | 2 | 50 |
| Grid | coarse \(\Delta=2\) | refined around QUICK centers |
| Primary metrics | descriptive | paired CI / win rates |
| `smooth_smooth` | cGCV often better MSE | **equivalent** (CI includes 0) |
| `strong_aniso` | similar MSE | oracle lower MSE on average, but **edge/flat artefacts** |

---

## Paired results (cGCV_default − tt_grid_oracle)

Convention: \(\Delta < 0\) ⇒ cGCV better.

### `smooth_smooth` (\(n=50\))

| Metric | mean | 95% CI | prop(cGCV better) |
|--------|-----:|--------|------------------:|
| \(\Delta\) test-MSE | −0.00056 | [−0.00124, 0.00012] | 0.52 |
| \(\Delta\) ISE | −0.00067 | [−0.00123, −0.00011] | 0.56 |
| \(\Delta\) TT-gGCV | −0.00037 | [−0.00046, −0.00027] | 0.86 |
| \(\|\Delta\log_{10}\lambda\|_2\) | 0.26 | — | — |

Relative \(\Delta\)MSE ≈ **−0.37%** of oracle MSE. Lambdas close. Surfaces mostly interior (edge rate 0.10).  
**Interpretation:** apparent QUICK “cGCV predicts better” **does not persist** as a material paired gain; methods are predictively equivalent. Mild ISE edge for cGCV; GCV gaps tiny vs MC noise.

Why QUICK looked stronger: small \(R\), coarser grid, smaller \(n_{\mathrm{test}}\) — sampling noise + oracle overfit to noisy GCV, not a stable effect.

### `strong_aniso` (\(n=50\))

| Metric | mean | 95% CI | prop(cGCV better) |
|--------|-----:|--------|------------------:|
| \(\Delta\) test-MSE | **+0.025** | [0.002, 0.048] | 0.30 |
| \(\Delta\) ISE | **+0.026** | [0.002, 0.050] | 0.28 |
| \(\Delta\) TT-gGCV | −0.019 | [−0.024, −0.014] | 0.88 |
| \(\|\Delta\log_{10}\lambda\|_2\) | **3.90** | — | — |

Relative \(\Delta\)MSE ≈ **+5.5%** (oracle lower MSE). But:

- oracle **edge rate 0.86**, mean flat 1% fraction **0.59**;
- mean \(\log_{10}(\lambda_1/\lambda_2)\): cGCV ≈ **+2.55** (correct sign 86%), TT-grid oracle ≈ **0.01** (correct sign only **32%**).

**Interpretation:** large λ disagreement. Oracle can win MSE by sitting on refined-grid boundaries / flat plateaus, while **failing anisotropy recovery**. cGCV better matches the intended directional roughness. Do **not** read the MSE gap as a clean win for joint TT-gGCV selection.

---

## Initialization stability

Across low/central/high cGCV starts:

| Scenario | median sd(test-MSE) across starts | mean range |
|----------|----------------------------------:|-----------:|
| smooth_smooth | 0.00089 | 0.0020 |
| strong_aniso | 0.00112 | 0.010 |

`cGCV_best_global` (min TT-gGCV among starts, **no test leakage**) vs `cGCV_default`:

- `smooth_smooth`: mean \(\Delta\)MSE ≈ 0 (CI includes 0)
- `strong_aniso`: mean \(\Delta\)MSE ≈ 0.00085 (CI includes 0)

**Multistart does not deliver material predictive gains** at this design size.

---

## Anisotropy (`strong_aniso`)

Truth: margin 2 more oscillatory ⇒ expect \(\lambda_1/\lambda_2 > 1\) (\(\log_{10}\) ratio \(>0\)).

| Method | mean \(\log_{10}(\lambda_1/\lambda_2)\) | P(sign correct) |
|--------|----------------------------------------:|----------------:|
| cGCV_default | 2.55 | **0.86** |
| cGCV_best_global | 2.61 | **0.88** |
| tt_grid_oracle | 0.01 | 0.32 |
| full_tp_gcv | 0.38 | 0.44 |

cGCV identifies the **direction** of anisotropy far more reliably than the refined TT-gGCV grid in this FULL run.

---

## Rank (secondary)

Restricted-rank probes (first 10 reps) remain diagnostic only; primary comparisons use sufficient rank. Inf/edge behaviour under restriction was already flagged in Phase 1 — do not equate full-TP EDF with TT GDF.

---

## Cost

- FULL wall ≈ **17 min** for 100 replicate-keys on this machine.  
- Per-rep cost dominated by refined grid × \((1+M_{\mathrm{search}})\) ALS fits + \(M_{\mathrm{final}}\) re-scores.  
- Still far heavier than a single outer-simultaneous cGCV fit → acceptable as **oracle/diagnostic**, not as default selector.

---

## Failures / nonconvergence

- `phase1b_failures.csv`: **not created** (zero hard failures).  
- Some GCV cells / ratios non-finite on flat anisotropic surfaces (contributes to Inf in |ΔGCV| summaries); recorded via `valid` / edge flags rather than silent drops.

---

## Limitations

1. \(d=2\) only; Gaussian only.  
2. Refined grid still discrete; edge-heavy on `strong_aniso`.  
3. MC GDF noise (\(M_{\mathrm{final}}=20\)) remains non-negligible for tiny GCV gaps.  
4. Nelder benchmark only on first 10 reps.  
5. Exit code 1 on orchestrator once due to `tee` before `mkdir` (run still completed; CSVs written).

---

## Decision

### **CLOSED — KEEP cGCV** (2026-08-16)

Product path locked. No further Monte Carlo of this type. TT-gGCV remains internal diagnostic only (not exported; not for current manuscript).

| Option | Verdict |
|--------|---------|
| **KEEP cGCV** | **CLOSED YES** — predictively equivalent on `smooth_smooth`; on `strong_aniso` recovers anisotropy and is stable across starts |
| **IMPROVE INITIALIZATION** | **No** — multistart TT-gGCV re-ranking does not materially improve test-MSE |
| **DEVELOP GLOBAL METHOD** | **No** — MSE gains of the grid “oracle” on `strong_aniso` coincide with edge/flat artefacts and **worse** anisotropy recovery |
| **STOP GLOBAL TRACK** (as replacement) | **CLOSED** as product replacement; **keep** TT-gGCV MC code + 59 tests as **internal diagnostic** for ALS / rank / init / parameterization audits |

Aligns with Phase 1 QUICK **CONDITIONAL GO**, now closed by \(R=50\): differences in \(\lambda\) are real; predictive superiority of joint TT-gGCV is **not** established as a reason to replace cGCV.

---

## Tests (pre-FULL)

```text
test-global-gcv-lab.R:     FAIL 0 · PASS 23 · SKIP 0
test-global-gcv-phase1.R:  FAIL 0 · PASS 17 · SKIP 0
test-global-gcv-phase1b.R: FAIL 0 · PASS 19 · SKIP 0
```

Total relevant: **FAIL 0 · PASS 59 · SKIP 0**

---

## Files created / modified

**Added**

- `inst/benchmarks/global_gcv/phase1b_config.R`
- `inst/benchmarks/global_gcv/phase1b_helpers.R`
- `inst/benchmarks/global_gcv/phase1b_run_full.R`
- `inst/benchmarks/global_gcv/phase1b_analyze.R`
- `inst/benchmarks/global_gcv/PHASE1B_FULL_REPORT.md`
- `tests/testthat/test-global-gcv-phase1b.R`
- `results/phase1b_full/*`, `results/phase1b_smoke/*`
- `figures/phase1b_full/*`, `figures/phase1b_smoke/*`

**Modified**

- `R/global_gcv_lab.R` — no further change required for 1B beyond Phase 1 helpers already present
- `inst/benchmarks/global_gcv/README.md` — document Phase 1B commands

QUICK `results/phase1/` **not overwritten**.
