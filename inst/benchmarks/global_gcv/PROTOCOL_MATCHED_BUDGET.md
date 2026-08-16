# Matched-budget protocol — Sobol vs grid / random / LHS

**Status:** READY TO RUN (algorithmic gate first; statistical MC later)  
**Date:** 2026-08-16  
**Product lock:** KEEP `cGCV` (unchanged). This protocol does **not** reopen product λ selection.  
**Related:** `GATE_REPORT.md` (v0/v1 PASS) · Brain `TTPsplines-global-gGCV-architecture`

---

## Research question

> Does Sobol exploration locate the near-optimal region of the **joint TT-gGCV** criterion with fewer evaluations than a Cartesian grid, independent random search, or Latin hypercube sampling?

**Not** the first question: whether Sobol improves predictive performance vs `cGCV`.  
**First** validate that it improves **optimization of the same criterion** \(Q_r(\boldsymbol\theta)\).

Three separable questions:

1. **Coverage** — Does low-discrepancy sampling enter the 1% valley more often / sooner?
2. **Refinement / adaptation** — Does `nlminb` (and edge expand) reduce regret beyond pure exploration?
3. **Statistical consequence** — Do algorithmic gains survive data variability? (only after algorithmic GO)

---

## 1. Objective (fixed rank)

With TT rank \(r\) fixed:

\[
\widehat{\boldsymbol\theta}
=
\arg\min_{\boldsymbol\theta\in\Theta}
Q_r(\boldsymbol\theta),
\qquad
\theta_j=\log_{10}\lambda_j,
\]

\[
Q_r(\boldsymbol\theta)
=
\frac{
n\,\mathrm{RSS}_r(\boldsymbol\theta)
}{
\bigl[n-\widehat{\mathrm{GDF}}_r(\boldsymbol\theta)\bigr]^2
}.
\]

Initial box:

\[
\Theta=[-5,5]^d.
\]

Expand a face by margin only if that face remains in the near-optimal region (M6).  
**Do not** let DMRG / rank adaptation change \(r\) during search. Rank comparison is a later experiment.

---

## 2. Methods

| ID | Method | Role |
|----|--------|------|
| M0 | `cGCV` | Operational baseline (cost / anchor / predictive reference). Not the same eval type as global designs. |
| M1 | Cartesian grid | Reference only for \(d=2,3\) (e.g. \(11^d\)). Not scalable. |
| M2 | Random search | IID uniform in \([-5,5]^d\). Tests “many points alone”. |
| M3 | Latin hypercube | Same \(N\) as Sobol. Marginal stratification vs joint low discrepancy. |
| M4 | Sobol pure | \(\widehat\theta=\arg\min_{m=1,\ldots,N} Q(\theta^{(m)})\). |
| M5 | Sobol + refine | Sobol → diverse elites → bounded `nlminb` → final re-eval. **Primary candidate.** |
| M6 | Sobol adaptive | Batched Sobol → edge probes → directed expand → regional enrich → refine → \(\mathcal R_{0.01}\). Full stack. |

---

## 3. Matched budget

Global methods receive the **same number of initial exploration points**:

\[
N\in\{16d,\,32d,\,64d\}.
\]

| \(d\) | Low | Mid | High |
|------:|----:|----:|-----:|
| 2 | 32 | 64 | 128 |
| 3 | 48 | 96 | 192 |
| 5 | 80 | 160 | 320 |
| 6 | 96 | 192 | 384 |

Log **separately**, compare on **total** cost:

- `n_explore`
- `n_refine`
- `n_edge`
- `n_expand_strip`
- `n_reeval`
- `n_als_total`
- `n_failed` (non-convergent / non-finite \(Q\))
- wall time

Primary efficiency metric uses total evaluations (and ALS fits), not initial \(N\) alone.

---

## 4. Scenarios

| ID | Setting | Goal |
|----|---------|------|
| A | \(d=2\) `smooth_smooth` | Interior min; surface viz; grid equivalence |
| B | \(d=2\) anisotropic | Recover \(\lambda_1\gg\lambda_2\); tilted valleys; λ-distance ≠ GCV equivalence |
| C | \(d=3\) strong anisotropy | Edge probes; directed expand; \(\theta_j\lesssim c\); weak ID |
| D | \(d=5\) gradual anisotropy \(\theta_1<\cdots<\theta_5\) | Scalability; seed stability; no grid |
| E | Weak / flat GCV | Return stable near-optimal **region**; avoid false precision |

Reuse existing lab designs where available (`.tt_lab_phase1_make_design_dn`, etc.).

---

## 5. Fair evaluation at each \(\boldsymbol\theta\)

```r
lambda <- 10^theta
init <- clone(common_initial_cores)
fit <- ttps(..., rank = r, lambda = lambda, init = init)
```

Policy for the **main** comparison:

- `cold_common` (same initial cores for all methods / points)
- same MC bank within a replicate
- same ALS tolerance / `max_sweeps`
- RNG preserved; per-point cache; convergence logged

Continuation / warm-start: separate ablation only — do **not** mix into the primary head-to-head.

---

## 6. Monte Carlo GDF control

\[
M_{\mathrm{search}} < M_{\mathrm{final}}
\qquad
\text{(e.g. } 20 < 200\text{)}.
\]

During search: one MC bank for all points (CRN); store SE(GDF); do not refresh the bank mid-replicate.  
Finalists: \(M_{\mathrm{final}}\); common bank; alternate bank; ALS multistart.  
Do not promote winners on \(M_{\mathrm{search}}\) alone.

---

## 7. Sobol design (R)

```r
U <- qrng::sobol(n = N_max, d = d, randomize = "Owen", seed = seed)
theta <- sweep(U, 2, upper - lower, `*`)
theta <- sweep(theta, 2, lower, `+`)
```

Generate \(N_{\max}\) once; consume prefixes (`1:N`) for convergence curves without changing the sequence.

---

## 8. Diverse elites (before refine)

Do **not** take the raw top-\(K\) by \(Q\) (they may collapse to one pocket).

1. Sort by \(Q\).
2. Take best; exclude points within distance \(\delta\) (in \(\theta\) / \(\log_{10}\lambda\) space).
3. Take next; repeat until \(K\in\{5,\ldots,10\}\).
4. Run bounded `nlminb` from each elite.

**Hyperparameter to report:** \(\delta\) (default suggestion: \(0.5\)–\(1\) on the \(\theta\) scale).

---

## 9. Adaptive expansion (M6)

For \(\widehat\theta\), probe faces:

\[
Q(L_j,\widehat\theta_{-j}),
\qquad
Q(U_j,\widehat\theta_{-j}).
\]

Expand margin \(j\) if

\[
Q_{\mathrm{face}} \le 1.01\, Q_{\min}.
\]

Enrich only the new strip (Sobol or LHS). Stop when the face is no longer competitive, an interior min appears, profile stabilizes (`effectively_unpenalized`), or `max_expansions` is hit.

---

## 10. Reference optimum

| \(d\) | Reference |
|------:|-----------|
| 2 | Dense grid + local refine |
| 3 | Moderate grid + profiles + high-budget Sobol |
| 5–6 | Multi-seed high-budget Sobol + multistart refine + alternate MC banks |

Name:

> **high-budget global TT-gGCV reference**

Never “exact global optimum”.

---

## 11. Primary metrics

| Metric | Definition / role |
|--------|-------------------|
| 1% valley | \(Q(\widehat\theta_m)\le 1.01\,Q_{\mathrm{ref}}\); success rate \(P(\text{1% valley})\) |
| Relative regret | \(R_Q=(Q(\widehat\theta_m)-Q_{\mathrm{ref}})/Q_{\mathrm{ref}}\) |
| Efficiency | Evals / ALS / time until valley entry; failed evals |
| Anisotropy | \(\mathrm{sign}(\theta_j-\theta_k)\) and \(\theta_j-\theta_k=\log_{10}(\lambda_j/\lambda_k)\) |
| Stability | Across seeds; alternate MC bank; \(\mathcal R_{0.01}\); ID labels |

On flats, region / labels beat point \(\hat\theta\).

---

## 12. Secondary metrics

\(\|\widehat\theta-\theta_{\mathrm{ref}}\|\), test-MSE, ISE, RSS, GDF, ALS convergence, fitted-value distance, true-\(f\) recovery.

λ-distance is **not** primary (plateaus contain distant λ with similar \(Q\)).

---

## 13. Replications

### Algorithmic gate (first)

One fixed dataset per scenario; design seeds:

```text
1, 11, 21, 31, 41
```

Studies optimization only.

### Statistical Monte Carlo (later)

\(R=50\) datasets per scenario at mid budget — only after algorithmic GO / CONDITIONAL GO.

Do **not** start with \(R=50\).

---

## 14. Ablations

| Config | Question |
|--------|----------|
| Sobol pure | Is coverage enough? |
| Sobol + refine | Does `nlminb` help? |
| Sobol adaptive | Do edge probes / expand help? |
| Sobol without `cGCV` | Independent of anchor? |
| `cGCV` + refine | Is a local solution enough? |
| Sobol + `cGCV` anchor | Cost reduction? |
| Random + refine | Sobol > random? |
| LHS + refine | Sobol > marginal stratification? |

---

## 15. Hypotheses

- **H1** Sobol enters the 1% valley more often than random at matched \(N\).
- **H2** Advantage grows with \(d\).
- **H3** Sobol+`nlminb` has lower regret than Sobol pure.
- **H4** Edge probes correct artificial boundary minima.
- **H5** Sobol without `cGCV` finds the same valley (higher cost OK).
- **H6** GCV gains vs `cGCV` are smaller than λ differences (plateaus).

---

## 16. Decision gate

| Verdict | Criteria |
|---------|----------|
| **GO** | Adaptive / Sobol+refine: ≥4/5 seeds in valley; beats random/LHS on evals or stability; works without `cGCV`; recovers anisotropy; scales to \(d=5\); labels flats/boundaries honestly |
| **CONDITIONAL GO** | Better coverage but not cost; or advantage only under strong anisotropy |
| **NO-GO** | Random/LHS equivalent; and `cGCV` already in the same valley at much lower cost |

**Defendable conclusion (if GO):**

> Adaptive low-discrepancy exploration recovers near-optimal regions of the joint TT-GCV criterion more efficiently than Cartesian search and more reproducibly than random designs, while `cGCV` remains an efficient operational selector.

**Not** the target claim: “Sobol selects better λ than `cGCV`.”

---

## 17. Execution order (cost-aware)

Already covered by v0/v1 GATE (do not redo as primary work): scenarios A–C mechanics, M5 viability, `no_cgcv`, directed expand on \(d=3\) `strong_aniso`.

**Remaining priority for this protocol:**

```text
1. Matched-budget gate (d=2,3): Random / LHS / Sobol / Sobol+refine
   N in {16d, 32d, 64d}, seeds {1,11,21,31,41}, cold_common
2. Ablations (cGCV+refine, Sobol+/-anchor) if (1) is not NO-GO
3. Scenario E (flat) on d=2 — region stability
4. d=5 (+ d=4 smoke) vs high-budget reference — H2
5. M6 only if M5 leaves boundary residues
6. R=50 statistical MC only after GO / CONDITIONAL GO
7. Final algorithmic claim (product stays KEEP cGCV)
```

Parallel track (from `GATE_REPORT.md`): v1 structural validation on \(d=4\)–\(6\), then IFT.  
Do **not** treat “more fixed Sobol points alone” as the next scientific step.

---

## 18. Artifacts

| Artifact | Purpose |
|----------|---------|
| `run_matched_budget_gate.R` | Driver (implemented) |
| `results/matched_budget_{smoke,gate}/` | Per-seed / per-\(N\) ledgers + aggregate |
| `MATCHED_BUDGET_REPORT.md` | Inside results dir — verdict GO / CONDITIONAL / NO-GO |

**Driver notes (v0):**

- \(d=2\) only (grid reference). Higher \(d\) later.
- Sobol = package Joe–Kuo + Cranley–Patterson rotation (seeded); no `qrng` dependency.
- `seed_data` fixes dataset, common cores, and MC banks; `seed_design` only scrambles the design.
- Methods: `random`, `lhs`, `sobol_pure`, `sobol_refine` (+ `cgcv` baseline row).
- Optional env: `TT_GGCV_MB_METHODS`, `TT_GGCV_MB_SCENARIOS`, `TT_GGCV_MB_N_MULT`.
- **Bugfix (2026-08-16):** evaluator cache keys were `%.6f`, which collapsed `nlminb` finite-difference steps → refine never moved. Keys are now `%.12f`; refine objectives use `use_cache = FALSE`.
- After each run, \(Q_{\mathrm{ref}}\) is lifted to \(\min(Q_{\mathrm{grid}}, \min Q_{\mathrm{final}})\) so a coarse grid cannot create false FAIL when a design beats it.

```bash
TT_GGCV_MB_MODE=SMOKE Rscript inst/benchmarks/global_gcv/run_matched_budget_gate.R
TT_GGCV_MB_MODE=GATE  Rscript inst/benchmarks/global_gcv/run_matched_budget_gate.R
```
