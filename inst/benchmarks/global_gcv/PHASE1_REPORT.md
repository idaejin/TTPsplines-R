# Phase 1 report: validating global TT-gGCV as an oracle vs cGCV

**Status:** internal laboratory (QUICK mode unless noted). Not a public method.  
**Date:** 2026-08-16  
**Code:** `R/global_gcv_lab.R` (Phase-1 helpers) · `inst/benchmarks/global_gcv/phase1_*.R`  
**Artifacts:** `results/phase1/` · `figures/phase1/`

Terminology unchanged from Fase 0: **global TT-gGCV**, **Monte Carlo GDF**,  
\(\mathrm{GDF}_{TT}\neq\sum_k\mathrm{ed}_k\).

---

## Question

Can a fixed-rank TT fitter support a **joint** smoothing criterion (TT-gGCV via MC-GDF) that:

1. recovers **interior** \(\boldsymbol\lambda\) minima when they exist;
2. serves as an oracle to judge **outer simultaneous cGCV**;
3. reveals dependence on **cGCV starts / margin order**;
4. separates **rank** vs **\(\lambda\)** effects;
5. admits a **derivative-free joint optimiser** as a benchmark (not a product path)?

No TMB. No public API / vignette / manuscript changes.

---

## Design (QUICK)

| Knob | QUICK value |
|------|-------------|
| \(n=n_{\mathrm{test}}\) | 120 |
| \(k\) | 6 |
| Coarse \(\log_{10}\lambda\) | \(\{-5,-3,-1,1,3,5\}\) |
| \(M_{\mathrm{search}}/M_{\mathrm{final}}\) | 3 / 8 |
| MC bank | 20 shared Rademacher columns (CRN) |
| Replicates (compact MC) | 2 |
| Init policy (primary) | `cold_common` (cloned `tt_initialize`) |
| Continuation | one serpentine pass (sufficient rank) |

Scenarios (\(d=2\)):

| ID | Truth sketch | \(\sigma\) | Full-TP min (coarse) | Interior? |
|----|--------------|------------|----------------------|-----------|
| `smooth_smooth` | \(\sin(2\pi x_1)+\cos(2\pi x_2)\) | 0.35 | \((-1,-3)\) | **yes** |
| `smooth_rough` | low-freq \(x_1\) + high-freq \(x_2\) | 0.40 | \((-5,5)\) | no (edge) |
| `strong_aniso` | linear \(x_1\) + fast \(x_2\) | 0.45 | \((3,3)\) | **yes** |
| `weak_signal` | tiny signal | 1.00 | \((5,-1)\) | no (edge) |

Boundary conclusions were checked by widening the box by one coarse step (design scout).

---

## Methods compared

1. **full-TP GCV** (exact EDF; scattered design)  
2. **TT-gGCV grid oracle** (cold_common, sufficient / restricted rank)  
3. **outer simultaneous cGCV** (public; several \(\lambda^{(0)}\) × margin orders)  
4. **Nelder–Mead on \(\theta=\log_{10}\boldsymbol\lambda\)** (multi-start; \(M_{\mathrm{search}}\) then \(M_{\mathrm{final}}\))

Objective evaluations use a **cache** keyed by  
`dataset + rank + theta + init_policy + mc_bank + M`, with **common random numbers**.

---

## Results (QUICK)

### Interior minima

- Full-TP and TT-grid (sufficient) agree on **smooth_smooth** at \((-1,-3)\).  
- **strong_aniso** has an interior full-TP min at \((3,3)\); TT-grid (coarse) preferred \((1,5)\) — same high-\(\lambda\) plateau (flat 1% region size 8/36).  
- Edge minima remain for rough / weak-signal designs (not forced to be interior).

### cGCV vs TT-gGCV oracle (compact MC, mean over 2 reps)

| Scenario | mean \(\Delta\) TT-gGCV (cGCV − oracle) | mean \(\|\Delta\log_{10}\lambda\|_2\) | ISE(cGCV) vs ISE(oracle) |
|----------|----------------------------------------:|--------------------------------------:|--------------------------|
| smooth_smooth | **−0.005** (cGCV slightly better on noisy GCV) | 2.60 | **0.036 vs 0.122** (cGCV better predictively) |
| smooth_rough | −0.002 | 1.19 | similar (~0.75 vs 0.74) |
| strong_aniso | −0.001 | 3.66 | similar (~0.33 vs 0.34) |
| weak_signal | −0.015 | 2.25 | similar / slightly better cGCV |

**Lambda distance is large; predictive / GCV gaps are small** on flat surfaces.  
Primary metrics (prediction + \(\Delta\)GCV) do **not** show systematic superiority of TT-gGCV over cGCV in QUICK.

### Basins (starts × orders)

- Outer simultaneous cGCV: **order `1-2` vs `2-1` usually near-identical** for low/mid starts.  
- **Start dependence is real** on `smooth_smooth`: high \(\lambda^{(0)}\) → larger \(\hat\lambda\), worse TT-gGCV and RMSE than low/mid starts.  
- `strong_aniso` / `weak_signal`: large spread in \(\log\lambda\) but tiny spread in GCV/ISE (flatness).  
- Trajectories recorded from `fit$history` / `cgcv$trace` when present (`phase1_cgcv_trajectories.csv` may be partial).

### Rank

- Restricted \(r=1\) shifts grid minima and can yield **Inf / unstable** TT-gGCV cells (e.g. `smooth_rough`).  
- Do **not** read full-TP EDF as TT GDF under rank restriction (confirmed again).

### Cold vs continuation

- Serpentine continuation vs cold_common: deltas exist; some non-finite cells under continuation (ALS basin / GDF failures).  
- Cold_common remains the **fair** \(\lambda\) comparison policy.

### Joint optimisation (Nelder–Mead)

- Viable as a **benchmark**: multi-start required; some starts hit box faces.  
- Best starts often recover near-cGCV / near-grid quality, not a clear predictive win over cGCV in QUICK.  
- Cost: tens of TT-gGCV evaluations per start (each \(\approx 1+M\) ALS fits).

---

## Cost

| Step | QUICK wall (approx) |
|------|---------------------|
| designs | ~0.6s |
| surfaces (after fix) | ~21s |
| basins | ~2.8s |
| optimize | ~7.7s |
| compact MC | ~13s |
| **First orchestrated run** | ~32s (surfaces initially failed; fixed and re-run) |

TT-gGCV remains an **oracle / diagnostic**, not a public default.

---

## Limitations

1. QUICK only (\(n=120\), coarse grid, \(M=3/8\), 2 MC reps).  
2. MC noise in the objective; CRN + cache mitigate but do not remove it.  
3. Fixed-\(\lambda\) ALS `converged=TRUE` remains weak.  
4. Some restricted-rank / continuation cells non-finite.  
5. No \(d=3\) in this phase.  
6. Surfaces script initially hit a `data.frame` length-0 bug on NULL `gdf_mc_se` from failed early returns — **fixed** in lab helpers (scalar coercion + complete early-return fields). Minimal productive fix; no public API change.

---

## Decision

### **CONDITIONAL GO**

| Criterion | Verdict |
|-----------|---------|
| GO (material reproducible gain vs cGCV) | **Not met** in QUICK |
| CONDITIONAL GO (λ/basin diffs; flat surfaces; small predictive impact) | **Met** |
| NO-GO (cGCV always equivalent + cost useless) | **Not met** — oracle still informative for basins / starts / rank |

**Recommendation for next phase**

1. Run `TT_GGCV_QUICK=false` (finer grid, larger \(M\), more reps) on the two interior scenarios.  
2. Keep cGCV as the public selector; use TT-gGCV as **oracle scoring** of cGCV starts and as a Paper-2 diagnostic.  
3. Consider TMB / implicit GDF **only if** FULL mode shows systematic predictive gains or if oracle cost becomes the bottleneck for theory checks — not to replace cGCV yet.  
4. Harden ALS convergence reporting if restricted-rank Inf cells persist (documented; not a broad refactor).

---

## Tests

```text
NOT_CRAN=true test-global-gcv-lab.R      → FAIL 0 | PASS 23 | SKIP 0
NOT_CRAN=true test-global-gcv-phase1.R  → FAIL 0 | PASS 17 | SKIP 0
```

## How to reproduce

```bash
# QUICK (default)
TT_GGCV_QUICK=true Rscript inst/benchmarks/global_gcv/phase1_run_all.R

# FULL (heavier)
TT_GGCV_QUICK=false TT_GGCV_R=5 Rscript inst/benchmarks/global_gcv/phase1_run_all.R
```

## Files created / modified (Phase 1)

**Modified**

- `R/global_gcv_lab.R` — Phase-1 helpers (designs, cache, grid, optim, flatness); early-return / scalar hardening  
- `DESCRIPTION` — already listed `global_gcv_lab.R` from Fase 0  

**Added**

- `inst/benchmarks/global_gcv/phase1_config.R`  
- `inst/benchmarks/global_gcv/phase1_designs.R`  
- `inst/benchmarks/global_gcv/phase1_surfaces.R`  
- `inst/benchmarks/global_gcv/phase1_basins.R`  
- `inst/benchmarks/global_gcv/phase1_optimize.R`  
- `inst/benchmarks/global_gcv/phase1_mc.R`  
- `inst/benchmarks/global_gcv/phase1_run_all.R`  
- `inst/benchmarks/global_gcv/PHASE1_REPORT.md`  
- `tests/testthat/test-global-gcv-phase1.R`  
- `results/phase1/phase1_*.csv`  
- `figures/phase1/phase1_*.png`  

Fase 0 artifacts under `results/` / `figures/` (non-phase1) preserved.
