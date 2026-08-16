# Evidence report — Margin Activity Path

**Date:** 2026-08-17  
**Status:** PRELIMINARY (small MC + demos; not Paper-1 claim)  
**Code:** `run_evidence_mc.R` · API `tt_margin_activity_path()` · `tt_margin_drop_test()`  
**Protocol:** `PROTOCOL_MAP.md`

---

## Verdict

| Claim | Evidence | Status |
|-------|----------|--------|
| Path + CV 1-SE recovers clear additive active set and improves holdout MSE vs full \(d\) cGCV | **Strong on favorable DGPs** (B=20 additive; reprex; demos) | PRELIMINARY → leaning ESTABLISHED *for this DGP class* |
| Same for pure product interaction \(\sin(2\pi x_1 x_2)\) | Good TPR; occasional extra margin; MSE gain smaller | PRELIMINARY |
| Nested drop test p-values are a reliable selector | High TPR but **inflated FPR** vs path | Diagnostic only |
| Method is ready as Paper-1 product claim | No — scope, designs, and baselines incomplete | REJECTED for Paper-1 core |

**Practical reading:** use **Margin Activity Path (1-SE)** as screening; treat drop-test p-values as secondary.

---

## 1. Sources of evidence

| Source | What | Strength |
|--------|------|----------|
| Unit tests | API smoke; top-2 scores / keep x1 | Very weak (single seed) |
| `reprex_margin_activity_path.R` | Fixed seeds; additive2, \(d=6\) | Weak (1 seed) |
| `standalone_margin_activity_path.R` | Same DGP, user-facing | Weak (1 seed) |
| `run_map_example.R` / early MAP figs | \(n_{\mathrm{signal}}=2\) or \(5\), \(d=10\) | Weak–moderate |
| **`run_evidence_mc.R` (this report)** | B=20 × 2 DGPs | Moderate gate |

---

## 2. Mini MC (primary new evidence)

**Design:** \(n=200\), \(d=6\), \(r=2\), \(k=5\), \(\sigma=0.3\), backend R, B=20 seeds per DGP.

| DGP | Truth |
|-----|--------|
| `additive2` | \(y=\sin(2\pi x_1)+\sin(2\pi x_2)+\varepsilon\) |
| `interaction` | \(y=\sin(2\pi x_1 x_2)+\varepsilon\) |

### Summary (`evidence_mc_summary.csv`)

| dgp | exact_set | TPR path | FPR path | TPR drop | FPR drop | MSE full | MSE path | MSE path / full |
|-----|----------:|---------:|---------:|---------:|---------:|---------:|---------:|----------------:|
| additive2 | **1.00** | **1.00** | **0.00** | 1.00 | 0.16 | 0.386 | **0.113** | **0.30** |
| interaction | 0.90 | **1.00** | 0.06 | 1.00 | 0.12 | 0.187 | 0.126 | 0.67 |

### Interpretation

1. **Additive (favorable):** path recovers exact `{x1,x2}` in **20/20** seeds; FPR=0; holdout MSE ~**3.4×** better than full cGCV on average.  
2. **Interaction:** still TPR=1; exact set 18/20 (2 seeds keep an extra null); MSE gain real but milder (~1.5×). Partial-range screening still “sees” both factors.  
3. **Drop test:** never misses signals here (TPR=1) but keeps ~12–16% of nulls on average → **liberal**; do not use alone for sparse selection.

Raw rows: `results/evidence_mc_raw.csv`.

---

## 3. Single-seed demos (consistent with MC)

| Run | Path keep | MSE full → path |
|-----|-----------|-----------------|
| Reprex (seeds 20260817/18) | x2,x1 | 0.429 → 0.123 |
| Standalone (same) | x2,x1 | 0.429 → 0.123 |
| Earlier \(d=10\), \(n_{\mathrm{signal}}=5\) demo | top signal block via CV | large MSE drop (order-of-magnitude in that script) |

Linear lasso baseline (benchmark script): can select the right *names* but **predicts poorly** under nonlinear DGP — selection ≠ TT fit quality.

---

## 4. What is *not* evidenced

- Correlated covariates / confounding  
- Weak SNR, many weak signals  
- High \(d\) (20–50) with sparse truth  
- Non-Gaussian families  
- Fair competitor: group-lasso on B-spline bases, gam `select=TRUE`, etc.  
- Calibration of nested F p-values (type-I under global null)  
- Stability of \(S_j\) under reference point (median vs mean)  
- Joint choice of \(r\) after screening  

---

## 5. Epistemic update (vs PROTOCOL §9)

| Claim | Before | After mini-MC |
|-------|--------|----------------|
| Additive path separates signal | PRELIMINARY (1 seed) | **PRELIMINARY+** (B=20 exact recovery) |
| Interaction path | untested | PRELIMINARY (good TPR, some FPR) |
| Drop test as formal cut | diagnostic | **Confirmed liberal** — keep as diagnostic |
| Multi-seed MC | NEXT | **Started** (gate B=20); enlarge before ESTABLISHED |

---

## 6. Recommendation

1. **Ship / use:** Margin Activity Path + 1-SE as **practical screening** for additive-ish sparse margins.  
2. **Do not claim:** calibrated inference from drop-test p-values; universality; Paper-1 centrality.  
3. **Next evidence gate:** B≥100; add correlated design + weak-signal; report CI on TPR/FPR/MSE ratio; optional group-lasso baseline.

---

## Reproduce

```bash
Rscript inst/benchmarks/margin_path/run_evidence_mc.R
Rscript inst/examples/reprex_margin_activity_path.R
```
