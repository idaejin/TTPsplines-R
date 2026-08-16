# Global λ optimizer — GATE report

**Date:** 2026-08-16  
**Product λ path:** **KEEP cGCV** (unchanged)  
**Experimental APIs (not exported):**
- v0: `tt_global_lambda_optimize()`
- v1: `tt_global_lambda_optimize_v1()` — adaptive box, near-optimal region, λ-ID labels

**Artifacts:** `results/global_opt_{gate,seeds,ablation,seed31,d3_gate,d3_aniso_expand,v1_d3_aniso}/`

---

## Canonical claim (locked)

> **Global-λ optimizer v0: GATE PASS on \(d=2\) and \(d=3\).**  
> **Global-λ optimizer v1 (adaptive): GATE PASS on \(d=3\) `strong_aniso`.**  
> Sobol + refinement recovers near-optimal TT-gGCV regions. v1 auto-expands only the truncated margin via edge probes, reports \(\mathcal R_{0.01}\), and labels flat low-λ plateaus (`effectively_unpenalized` / `weakly_identified`). Public product remains **KEEP cGCV**.

| Layer | Status |
|-------|--------|
| Productive method | **KEEP cGCV** |
| v0 experimental optimizer | **GATE PASS** (\(d=2\), \(d=3\)) — CLOSED |
| v1 adaptive region diagnostics | **GATE PASS** on \(d=3\) `strong_aniso` (edge probe + ID labels) — PRELIMINARY beyond that case |

**Three-way separation (locked):** productive `cGCV` selector ≠ validated experimental global-λ optimizer ≠ methodological line (v1).

---

## v1 GATE — \(d=3\) `strong_aniso` (PASS)

Start box: \([-5,5]^3\). Driver: `run_v1_gate_d3_aniso.R`.

| Check | Result |
|-------|--------|
| Expands **only** margin 2 | YES (2 expansions → \(\theta_2\) lo \(-5\to-8\to-11\)) |
| Does **not** expand margins 1, 3 | YES (`force_lo=0,1,0`) |
| Recovers valley \(\theta_2\lesssim-4.5\) | YES (\(\hat\theta_2\approx-5.29\); \(I_2=[-11,-3.32]\)) |
| λ₂ label | `effectively_unpenalized` (flat profile for \(\theta_2\le-4\)) |
| Near-optimal region | 20 points; outer_search_stable = TRUE |
| Batches used | 3 (64+64+128); no spurious extra margins |

**Mechanism that unlocked expand:** Sobol+refine best sat at \(\theta_2\approx-3.3\) (interior). Edge probe holding \(\theta_{-2}\) fixed showed \(Q(\theta_2=-5)\) still within the 1% band → directed expand of margin 2 only.

**Interpretation (unchanged from v0 lock):** report \(I_2\) / plateau, not point precision of \(\hat\theta_2\).

---

## Interpretation lock — \(d=3\) `strong_aniso`

Numerical minimizer (v0 directed expand):

\[
\widehat{\boldsymbol\theta}=(0.04,-5.82,-1.94).
\]

**Do not** treat \(\widehat\theta_2\) as a precise estimate. The scientifically relevant statement is

\[
\theta_2\lesssim -4.5,
\]

a region of essentially null penalty on margin 2 where GCV, GDF, and ISE change little.

| Fact | Status |
|------|--------|
| Numerical interior minimum exists | YES (in expanded box) |
| Evidence of asymptotic descent to \(-\infty\) | **NO** (flat plateau) |
| \(\lambda_2\) identification | **Weak / effectively unpenalized** |
| Near-optimal region > point minimizer | **YES** |

---

## Evidence summary

### \(d=2\) (v0)

- Valley 1% + anisotropy: **PASS** (`smooth_smooth`, `strong_aniso`).
- Ablation `no_cgcv`: **PASS** (`refined_from_sobol`).
- Seed 31: strict ranking fail → **SOFT_PASS** (both candidates in 1% valley).

### \(d=3\) (v0)

- `smooth_smooth`: valley PASS; explore **138 vs 343**; ALS **2189 vs 4459**.
- `strong_aniso` (symmetric box): edge artefact at \(\theta_2=-5\).
- Directed expand \(\theta_2\in[-8,5]\): **INTERIOR_MIN** in flat valley; winner `refined_from_sobol`.

### \(d=3\) (v1 adaptive)

- Auto-expand from \([-5,5]^3\): **PASS** (only \(\theta_2\); label `effectively_unpenalized`).
- Tests: `tests/testthat/test-global-lambda-optimize-v1.R` (19 PASS).

---

## Metrics (locked)

| Metric | Role |
|--------|------|
| 1% valley | Primary |
| Anisotropy | Primary |
| Soft MC ranking / tie rule | Flip OK if both in valley + aniso |
| Directed margin expand | Expand only affected coords (edge probe + soft tol) |
| Report region on flats | Required for honest λ diagnostics |
| λ identifiability labels | `well_identified` / `weakly_identified` / `boundary_unresolved` / `effectively_unpenalized` |
| \(d_\theta\) / point \(\hat\theta\) alone | Secondary on plateaus |

---

## After v1 (priority)

1. Validate adaptive path on \(d=4\)–\(6\) (known designs first).  
2. Matched-budget Sobol vs random/LHS — [`PROTOCOL_MATCHED_BUDGET.md`](PROTOCOL_MATCHED_BUDGET.md) (algorithmic efficiency; does not reopen KEEP cGCV).  
3. Then IFT / cheaper GDF (computational route).  
4. Surrogate / adaptive Sobol for \(d=5\)–\(10\) only after (1)–(2).

### \(d=4\)–\(6\) validation plan

| Step | Design | Pass criteria (v1) |
|------|--------|--------------------|
| d=4 SMOKE | `.tt_lab_phase1_make_design_dn` `strong_aniso` + `smooth_smooth` | Finite; region; aniso: expand only margin 2 if needed + weak/unpen ID |
| d=4 GATE | Same, full batches | Same as SMOKE at d=3-comparable budgets |
| d=5 | Package `friedman` (known) | Finite; region; no spurious expand of clearly interior margins |
| d=6 | `design_dn(d=6)` known scenarios | Same structural checks; cost reported |

Driver: `run_v1_gate_d4.R` (`TT_GGCV_V1_MODE=SMOKE|GATE`).

#### d=4 SMOKE `strong_aniso` (2026-08-16) — FAIL (diagnostic)

Artifacts: `results/global_opt_v1_d4_smoke/strong_aniso/`.

| Check | Result |
|-------|--------|
| Finite / has region | TRUE |
| Search recovers low-\(\theta_2\) point | YES (\(I_2\approx\{-6.63\}\) after expand) |
| Expand includes margin 2 | YES |
| Expand only margin 2 | **NO** (also margin 4) → `expand_clean=FALSE` |
| λ₂ label weak/unpen | **NO** (`well_identified` on singleton region) |
| Final winner | `cgcv_anchor` at \(\theta_2\approx+4\) (not the search valley) |
| Wall time | ~39 min (SMOKE budgets) |

**Interpretation (PRELIMINARY):** at \(d=4\) the adaptive search can reach a truncated low-\(\theta_2\) slab, but (i) final MC re-eval still prefers the cGCV anchor, and (ii) soft edge rules still expand an extra margin. Next diagnostics before GATE: keep valley candidates in the final shortlist; tighten expand to edge-probe-only when \(d\ge4\); optionally drop `include_cgcv_anchor` from winner selection when a probed face is competitive.

Do **not** treat “more fixed Sobol points” as the next scientific step.

Public product remains **KEEP cGCV**; v0/v1 stay experimental / diagnostic.

Brain: [[TTPsplines-global-gGCV-architecture]]

---

## Reproduce

```bash
TT_GGCV_OPT_MODE=GATE Rscript inst/benchmarks/global_gcv/run_global_lambda_optimize.R
TT_GGCV_OPT_MODE=GATE Rscript inst/benchmarks/global_gcv/run_global_lambda_optimize_d3.R
Rscript inst/benchmarks/global_gcv/run_d3_aniso_expand.R
Rscript inst/benchmarks/global_gcv/run_v1_gate_d3_aniso.R
TT_GGCV_V1_MODE=SMOKE Rscript inst/benchmarks/global_gcv/run_v1_gate_d4.R
TT_GGCV_V1_MODE=GATE  Rscript inst/benchmarks/global_gcv/run_v1_gate_d4.R
TT_GGCV_MB_MODE=SMOKE Rscript inst/benchmarks/global_gcv/run_matched_budget_gate.R
TT_GGCV_MB_MODE=GATE  Rscript inst/benchmarks/global_gcv/run_matched_budget_gate.R
NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-global-lambda-optimize-v1.R")'
```
