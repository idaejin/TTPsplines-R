# Paper 1 — Simulation section draft (rank selection + GLAM vs TT)

**Status:** PRELIMINARY (2026-08-10).  
**Provenance:** Factual package results; interpretation for manuscript.

Depends on closed pipeline:

1. GLAM-vs-TT with `TT-minCV` / `TT-1SE` (+ oracle reference only)
2. Multi-seed (\(N=30\)) rank-selection audit
3. Rank × \(\lambda\) (fixed grid, then cGCV)

## Proposed structure

1. **Recovery / compression** — TT vs dense counts, CR.
2. **Rank selection** — CV-min vs 1-SE vs oracle (simulation only).
3. **Rank × smoothing** — \(r\) = structural capacity, \(\lambda\) = smoothness.
4. **Scalability** — GLAM feasible vs \(k^d\) cliff; TT continues.
5. **Families** — Gaussian primary; Poisson / Bernoulli as extensions.

## Central comparison table

\[
\boxed{
\text{Oracle TT}
\;|\;
\text{CV-min TT}
\;|\;
\text{1-SE TT}
\;|\;
\text{GLAM}
}
\]

Metrics per cell: RMSE (vs truth, simulation), `npar`, compression, time
(+ `select_time` for CV).

API: `compare_glam_tt_*()` → `summarize_glam_tt_compare()`.

### Snapshot (compact design, single seed)

| d | k | GLAM RMSE | Oracle r/RMSE | minCV r/RMSE | 1SE r/RMSE |
|---|---|-----------|---------------|--------------|------------|
| 3 | 8 | 0.240 | 3 / 0.075 | 3 / 0.075 | 3 / 0.075 |
| 5 | 6 | infeas. | 3 / 0.161 | 3 / 0.161 | 3 / 0.161 |
| 7 | 3 | 0.079 | 2 / 0.032 | 2 / 0.032 | 1 / 0.055 |

(Full table: `Rscript inst/examples/example_glam_gaussian_vs_tt.R`.)

## Multi-seed (\(N=30\), n=400)

| Surface | mean oracle r | mean 1SE r | mean RMSE oracle | mean RMSE 1SE | penalty | 1SE mode |
|---------|---------------|------------|------------------|---------------|---------|----------|
| Ishigami | 5.0 | 3.1 | 0.325 | 0.354 | +0.029 | r=3 (90%) |
| Sobol-g | 5.0 | 3.0 | 0.166 | 0.208 | +0.042 | r=3 (100%) |
| Friedman | 3.1 | 2.7 | 0.644 | 0.810 | +0.166 | r=3 (73%) |

**GLAM-vs-TT** cell \(d=3\), \(k=8\), 30 seeds:

- mean RMSE: GLAM 0.240 | Oracle 0.081 | minCV 0.084 | **1SE 0.086**
- 1SE penalty vs oracle ≈ **0.0045** (small)
- mean ranks: oracle 2.80 | minCV 2.80 | 1SE 2.63

**Claim (safe):** On the grid GLAM-vs-TT design, 1-SE preserves nearly all of the
oracle TT advantage over GLAM without using \(f_{\mathrm{true}}\).

**Caveat:** On scattered Friedman with fixed \(\lambda=1\), 1-SE is more
conservative and the truth-RMSE penalty is larger; cGCV changes the picture
(next section).

## Rank × \(\lambda\) (Friedman, n=400)

| Setting | Selected r (1SE) | RMSE vs truth |
|---------|------------------|---------------|
| \(\lambda=0.1\) fixed | 2 | 1.23 |
| \(\lambda=1\) fixed | 2 | 1.23 |
| \(\lambda=10\) fixed | 2 | 1.23 |
| `lambda="cGCV"` | **3** | **0.53** |

Fixed-\(\lambda\) grid also shows: at strong smoothing (\(\lambda=10\)), rank 3
RMSE worsens vs \(\lambda=1\) (1.06 vs 0.65), while rank 4 remains competitive —
higher \(r\) can **interact** with \(\lambda\), not substitute for it.

**Message for Paper 1:**

\[
r \neq \lambda \neq \mathrm{EDF}.
\]

Outer CV chooses structural capacity; inner cGCV (training-only) chooses
smoothness. Do not collapse them into one “complexity” knob.

## Wording for the manuscript

- Use **1-SE** (or min-CV) as the practical TT entry in GLAM-vs-TT tables.
- Keep **oracle** in a grey / footnote column: “simulation oracle; not a method.”
- Avoid: “1-SE proves the true rank is …”
- Prefer: “1-SE selects the smallest candidate whose CV error is within one
  SE of the minimum.”

## Scripts

- `inst/examples/example_glam_gaussian_vs_tt.R`
- `inst/benchmarks/benchmark_rank_select_multiseed.R` (`N_SEEDS=30` or `50`)
- `inst/benchmarks/benchmark_rank_lambda.R`
- Outputs under `inst/benchmarks/output/`

## Still open before locking the section

- [ ] Optional \(N=50\) confirmation (`N_SEEDS=50`)
- [ ] Poisson / Bernoulli analogues of the central table
- [ ] Figure panels: CV curves; RMSE vs npar for Oracle|minCV|1SE|GLAM
