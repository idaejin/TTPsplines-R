# Experimental report: global TT-gGCV via Monte Carlo GDF

**Status:** exploratory laboratory only. Not a validated public method.
**Date:** 2026-08-16  
**Scope:** Gaussian, fixed TT rank, fixed numeric \(\boldsymbol\lambda\), no TMB, no API / manuscript changes.

Terminology used here:

- **global TT-gGCV** — criterion \(n\,\mathrm{RSS}/(n-\widehat{\mathrm{GDF}})^2\) for the converged fixed-rank TT fitter
- **Monte Carlo GDF** — Hutchinson / Rademacher estimate of \(\mathrm{tr}(\partial\hat y/\partial y^\top)\)
- **Not** classical exact GCV; **not** \(\sum_k\mathrm{ed}_k\)

Code: `R/global_gcv_lab.R` (internal). Scripts: `inst/benchmarks/global_gcv/`.

---

## Phase 0 — API inspection (summary)

| Item | Finding |
|------|---------|
| Public fit | `ttps(..., lambda = numeric, optimizer = "ALS", family = gaussian())` |
| Internal engine | `tt_als_fit` → `tt_als_fit_sequential` when `lambda` is fixed |
| Rank | `rank` (scalar / length \(d-1\) / full chain); frozen across GDF refits |
| Init / warm start | `init =` list of cores from `tt_initialize()` or `fit$cores` (**supported**) |
| Seed | `tt_control(seed=)`; `initialize_tt_cores()` calls `set.seed` **without restoring** |
| ALS order | `tt_control(cgcv_margin_order=)` (also used under fixed \(\lambda\)) |
| Tol / sweeps / backend | `tol`, `max_sweeps`, `backend` / `ttps(backend=)` |
| Extract | `fitted()`, `$deviance` (Gaussian RSS), `$cores`, `$converged`, `$n_sweeps`, `$lambda` |
| Canonicalization | Left-orthogonal gauge only in **inference**; ALS does **not** re-orthogonalize between sweeps |
| Nondeterminism | Random init if `init=NULL`; OpenMP/`gram_threads`; seed side effects of init |

Plan executed: internal lab functions + `inst/benchmarks/global_gcv/` scripts + light `testthat` tests. No NAMESPACE exports.

---

## Answers to the scientific questions

### 1. Can global GDF be estimated stably?

**Mostly yes in this toy regime (CONDITIONAL).**

With warm-started ALS, shared Rademacher probes, and fixed \((r,\boldsymbol\lambda)\):

- Estimates are **reproducible** given the same probes.
- **Insensitive to `epsilon_rel`** over \(\{10^{-2},10^{-3},10^{-4}\}\) when rank is sufficient (map behaves linearly under warm start).
- **Sensitive to \(M\)** through Monte Carlo variance (expected): error vs full EDF falls from \(\approx 2.3\) (\(M=3\)) to \(\approx 0.22\) (\(M=20\)), with MC SE \(\approx 0.70\) at \(M=20\).

Caveat: fixed-\(\lambda\) ALS reports `converged=TRUE` even at `max_sweeps`. The lab trusts that flag and records `hit_max_sweeps` / last relative RSS as diagnostics.

### 2. Do forward and central agree?

**Yes here.** On the rank-sufficient case they matched to machine precision across the \(\epsilon\times M\) grid — consistent with an essentially linear local map. Prefer **forward** for cost in linear regimes; keep central as a nonlinearity check.

### 3. Does GDF match full-tensor EDF when rank is sufficient?

**Within Monte Carlo error, yes.**

| Quantity | Value |
|----------|------:|
| \(\mathrm{EDF}_{\mathrm{full}}\) | 6.945 |
| \(\widehat{\mathrm{GDF}}\) forward/central (\(M=20\)) | 6.725 |
| Abs. error | 0.221 (inside MC SE \(\approx 0.70\)) |
| \(\mathrm{RSS}_{\mathrm{TT}}\) vs \(\mathrm{RSS}_{\mathrm{full}}\) | identical |
| \(\|\hat y_{\mathrm{TT}}-\hat y_{\mathrm{full}}\|_2\) | \(\sim 10^{-13}\) |

So the **fit** matches the full TP; the **divergence estimate** is consistent up to Hutchinson noise.

### 4. What happens under rank restriction?

With \(r=1\) (deliberately too small for \(k=5\), \(d=2\)):

- RSS rises slightly; fitted surface departs from full TP (\(\mathrm{RMSE}\approx 0.022\), \(\mathrm{cor}\approx 0.994\)).
- \(\widehat{\mathrm{GDF}}\approx 4.38 < \mathrm{EDF}_{\mathrm{full}}\) — fewer effective degrees of freedom, as expected for a restricted estimator.
- Do **not** interpret this as failure of MC-GDF; the map \(\mathcal A_r\) changed.

### 5. Does the global surface have an interpretable interior minimum?

**Not on this \(n=80\) toy grid.** Both global TT-gGCV and full-TP GCV minimize at the **corner** \((\log_{10}\lambda_1,\log_{10}\lambda_2)=(-3,-3)\). The surface is smooth and qualitatively aligned with full-TP GCV, but this design prefers the lowest smoothing on the grid (undersmooth GCV behaviour / interactive truth). **No interior minimum to trust yet.**

### 6. Where does outer simultaneous cGCV land?

Outer cGCV \(\hat{\boldsymbol\lambda}\approx(10^{-2.35},\,10^{-2.62})\), near the low-\(\lambda\) region preferred by the surfaces, but **not** identical to the grid corner minimum. Treat as “same qualitative basin,” not as agreement of selectors.

### 7. Is the cost only a benchmark, or a practical method?

**Benchmark / oracle first.**

- This toy (\(n=80\), \(k=5\), \(d=2\), warm start) finished a 49-point surface with \(M=5\) forward probes in a few seconds.
- Scaling: each criterion evaluation costs \(\approx 1 + M\) (forward) or \(1 + 2M\) (central) full ALS fits. For production \(\lambda\) search in higher \(d\) / larger \(n\), that is heavy versus cGCV’s \(d\) conditional 1D searches.
- Recommendation: keep MC global GDF as a **diagnostic / validation oracle**, not as the default public \(\lambda\) engine.

### 8. Enough evidence to justify TMB / implicit differentiation later?

**CONDITIONAL GO** for a *secondary* research track:

| Evidence for | Evidence against / missing |
|--------------|----------------------------|
| GDF is estimable and probe-reproducible | Need cases with **interior** \(\lambda\) minima |
| Matches full EDF within MC noise when \(r\) sufficient | Need larger \(n\), \(d=3\), and rank-restricted surfaces |
| Forward≈central when map is linear | Need basin-jump stress tests (cold start, multi-start) |
| Warm start works with public `init=` | Fixed-ALS convergence flag is weak |
| Cost OK as oracle | Not yet competitive as primary selector |

TMB/implicit AD would be justified if (a) MC-GDF remains the scientific target, and (b) we need cheaper \(\partial\hat y/\partial y\) for optimization — **after** showing an interior, stable global surface on harder designs.

---

## Provisional lab settings (not public defaults)

From nested-probe sensitivity (still only exploratory):

- `scheme = "forward"` for cost (or `"central"` as check)
- `epsilon_rel = 1e-3`
- `M = 10`–`20` when comparing to full EDF

---

## Recommendation

### **CONDITIONAL GO** (toward a later TMB / implicit-diff study)

**Go if:** continue as an internal oracle / Paper-2 diagnostic, expand designs (interior minima, \(d=3\), cold-start basin audits), keep cGCV as the public \(\lambda\) method.

**No-go if:** the goal is to replace cGCV in the public API or manuscript with global TT-gGCV on present evidence alone.

---

## Files produced

- `R/global_gcv_lab.R` — internal estimators
- `tests/testthat/test-global-gcv-lab.R` — 23 unit tests (`NOT_CRAN=true`)
- `inst/benchmarks/global_gcv/run_*.R` — validation / sensitivity / surface
- `inst/benchmarks/global_gcv/results/*.csv`
- `inst/benchmarks/global_gcv/figures/*.png`

## Commands run

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-global-gcv-lab.R")'
# → FAIL 0 | PASS 23

Rscript inst/benchmarks/global_gcv/run_full_tensor_validation.R
Rscript inst/benchmarks/global_gcv/run_epsilon_M_sensitivity.R
Rscript inst/benchmarks/global_gcv/run_lambda_surface_d2.R
```

## Limitations (explicit)

1. Gaussian only; rank fixed by construction.
2. ALS `converged` flag is weak; basin changes must be diagnosed separately.
3. Toy surface has **boundary** GCV minimum — do not overclaim selector quality.
4. MC SE at moderate \(M\) still large relative to tiny EDF gaps.
5. No manuscript / vignette / public default changes (by design).
