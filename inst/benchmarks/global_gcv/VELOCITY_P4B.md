# P4B velocity report

**Date:** 2026-08-17  
**Machine:** Darwin arm64 · R 4.6.1  
**Scripts:** `run_p4b_probe_warm_velocity.R`, `run_p4b_lab_velocity.R`  
**Smoke gate:** `run_p4b_gate_smoke.R` → **PASS**

## Summary

| Lever | Effect (this machine) |
|-------|------------------------|
| `gdf_init=probe_warm` vs `cold` | ~1.0–1.4× on TT-gGCV; isolation asserts hold |
| `adaptive_fidelity` vs flat 40/60 sweeps | ~1.2–1.3× wall-clock on small Sobol+refine+final |
| Combined with P4A `Rcpp_fixed` | Multiplicative with (1+M) fits per θ |

**DECISION:** do **not** start `fused_blocked` yet — remaining time still includes search ALS count; capture pipeline savings first. Revisit after a fuller GATE re-run.

## `probe_warm` vs `cold`

| n | d | M | warm (s) | cold (s) | speedup | base_ok |
|--:|--:|--:|---------:|---------:|--------:|:-------:|
| 500 | 2 | 5 | 0.051 | 0.050 | 0.99× | TRUE |
| 500 | 3 | 5 | 0.107 | 0.108 | 1.01× | TRUE |
| 1000 | 3 | 8 | 0.149 | 0.210 | **1.41×** | TRUE |

GCV relative diffs ≪ 1%. GDF can differ at low M (different basins under cold init) — ranking still uses shared banks.

## Adaptive fidelity (lab optimize, smoke budgets)

| case | adaptive | elapsed (s) | sobol/refine/final sweeps | ranking_stable |
|------|:--------:|------------:|---------------------------|:--------------:|
| d=2 flat | no | 1.80 | 40/40/60 | TRUE |
| d=2 adaptive | yes | **1.35** | 12/25/50 | TRUE |
| d=3 flat | no | 27.05 | 40/40/60 | TRUE |
| d=3 adaptive | yes | **22.81** | 12/25/50 | TRUE |

Same winner family on d=2 (`refined_from_corner_lo`); d=3 winner source can shift under cheap search but final fidelity + alt bank remain stable here.

Note: adaptive may increase *counted* search ALS (refine cache keys differ by sweeps) while still reducing wall-clock via cheaper fits.

## Smoke gate checks (PASS)

- `fit_base$cores` unchanged under `probe_warm`
- Final candidates tagged `fidelity="final"`
- Budgets sobol=12, final=50
- Alt-bank ranking stable on ≥1 case
- Finite GCV + `winner_source` on d=2 and d=3 aniso smoke

## Next

1. ~~Optional fuller GATE (`smooth_smooth` / `strong_aniso`) with P4B defaults~~ **DONE — PASS**
2. Re-profile; only then `fused_blocked` / interface cache  
3. θ→θ `continuation` remains experimental (not in gate)

## Full GATE (2026-08-17) — **PASS**

Driver: `run_p4b_gate_full.R` · artifacts: `results/p4b_gate_gate/`  
Stack: `fit_backend=Rcpp_fixed`, `gdf_init=probe_warm`, `adaptive_fidelity=TRUE`  
Reference θ: locked 2026-08-16 GATE (re-scored at `M_final`, not fresh dense grid).

| Track | Scenario | in 1% valley | aniso | ranking | PASS |
|-------|----------|:------------:|:-----:|:-------:|:----:|
| v0 d=2 | smooth_smooth | ✓ | ✓ | strict | ✓ |
| v0 d=2 | strong_aniso | ✓ | ✓ | strict | ✓ |
| v0 d=3 | smooth_smooth | ✓ | ✓ | soft* | ✓ |
| v0 d=3 | strong_aniso | ✓ | ✓ (θ₂≈−4.3) | strict | ✓ |
| v1 d=3 | strong_aniso | ✓ (θ₂≈−8.1) | ID OK | strict | ✓ |

\*d=3 `smooth_smooth`: alt MC bank disagrees on winner, but both stay in 1% valley (same soft rule as locked v0 GATE).

**v1 highlights:** expands **only** margin 2; `lambda2_id=effectively_unpenalized`; `winner_source=region_point`.  
Wall-clock vs prior R-only GATE: d=2 ~38 s (was ~22 s with fewer ALS in R — comparable order); v1 ~17 min (was ~8 min — P4B + C++ path overhead on expand batches).

**Not claimed:** d=3 `strong_aniso` v0 still differs from symmetric-box grid best (known edge artefact); claim is valley + anisotropy sign, not point θ match.
