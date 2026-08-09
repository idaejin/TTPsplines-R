# Bernoulli PIRLS audit — ALS vs L-BFGS (fixed λ)

**Status:** diagnostic complete. Minimal FIX1+FIX2 + gate → `docs/BERNOULLI_PIRLS_STABILIZATION.md` (Outcome B).  
**Script:** `inst/benchmarks/benchmark_bernoulli_audit.R` (`--quick` / `--full`)  
**Results:** `inst/benchmarks/results/bernoulli_audit/` (do not overwrite)  
**Run:** 2026-08-09 (`--quick`: n=500, d=3, k=6, ranks 2–3)

Final technical question answered at §14–16.  
**Paper framing:** §−1 below; hub [[TTPsplines]] / [[TTPsplines-Paper-1]].

---

## −1. Epistemic framing (**DECISION** — how to write this)

**Do not sell:** “We discover an instability of Bernoulli TT models.”

**ESTABLISHED (logistic / IRLS literature):** extreme fitted probabilities can collapse PIRLS weights \(W_i=p_i(1-p_i)\to 0\) (complete / quasi-complete separation, near-separation, overconfident fits). Classical remedies: step-halving, trust regions, stronger regularization, Firth / bias-reduced likelihood.

**What is TT-specific (PRELIMINARY → interesting HYPOTHESIS):** that known mechanism interacts with flexible coefficient-TT structure and P-spline smoothing:

\[
\boxed{r \times \boldsymbol\lambda \times \text{logistic likelihood}}.
\]

In TT-P-splines this need not be classical separation in the original design; it can be **near-separation / overconfident fitting** from a flexible \(\Theta_{\mathrm{TT}}\). Higher \(r\) raises capacity; lower \(\lambda\) raises local freedom → \(|\eta|\gg 0\), \(p\simeq 0/1\), \(W\simeq 0\).

Two superimposed layers:

1. separation / extreme \(p\) / IRLS instability — **known**;  
2. nonconvex TT factorization + rank + directional smoothing — **TTPsplines-specific**, shapes *which* extreme-logit regions are reachable.

**Optimizer implication of this audit:** L-BFGS does **not** eliminate the statistical extreme-\(p\) problem; it follows a **different trajectory** in TT-core space and may generalize better even when ALS train \(L\) is similar or slightly better. FIX1 (consistent \(W/z\)) + FIX2 (true-objective step-halving) are **natural PIRLS stabilizations**, not TT hacks.

**POSITIONING prose:**

> As is well known for logistic regression, extreme fitted probabilities can lead to degenerating IRLS weights. In the TT-P-spline setting, this phenomenon interacts with the TT rank and the strength of smoothness regularization, requiring appropriate stabilization of the PIRLS updates.

**Worth pursuing later (not a headline claim yet):** \(r\uparrow,\lambda\downarrow\) ⇒ greater propensity to extreme-logit solutions — characterize *modulation*, do not rediscover IRLS pathology.

---

## 0. Setup

Shared across ALS and L-BFGS within each experiment:

- family = `binomial(link = "logit")`
- same `tt_initialize()` cores
- same bases / knots / difference penalties
- fixed isotropic `lambda` (default 1 unless noted)

Penalized objective (both methods):

\[
L(G,\alpha)
=
-\sum_i\bigl[y_i\log p_i+(1-y_i)\log(1-p_i)\bigr]
+\tfrac12\sum_k\lambda_k\,g_k^\top P_k g_k,
\qquad
p_i=\mathrm{logistic}(\alpha+f_i).
\]

ALS path: PIRLS working LS + weighted TT-ALS.  
L-BFGS path: direct minimization of \(L\) over packed cores (intercept fixed at GLM init unless warm-started).

---

## 1. Objective parity (Q1)

| rank | `tt_objective` | L-BFGS internal | rel discrepancy | η RMSE |
|------|---------------:|----------------:|----------------:|-------:|
| 2 | 290.5 | 290.5 | **0** | **0** |
| 3 | 280.3 | 280.3 | **0** | **0** |

**PASS** (tolerance \(10^{-8}\)).

**Decision:** **not A**. Both optimizers evaluate the **same** Bernoulli + TT penalty model. Benchmark gaps are algorithmic / path-dependent, not a silent objective mismatch.

---

## 2. PIRLS equations / code audit (Q2)

Package `glm_working()` for Bernoulli:

```r
mu <- plogis(eta)
mu <- pmin(pmax(mu, 1e-5), 1 - 1e-5)
var <- mu * (1 - mu)
w   <- pmax(var, 1e-4)
z   <- eta + (y - mu) / var   # divides by var, not w
```

Findings:

1. **μ clipping** to \([10^{-5},1-10^{-5}]\) — standard.
2. **W floor** at \(10^{-4}\) while **z uses unclipped `var`** — mild inconsistency when `var < 10^{-4}`.
3. Soft “damping” in `tt_pirls_fit` is **not** step-halving: if deviance rises by >25%, PIRLS **stops** and keeps the previous iterate.
4. When `init` is supplied, cores are **not** scaled by 0.05 (fair vs L-BFGS). When `init=NULL`, ALS shrinks init by 0.05 (asymmetric).

**Decision:** minor **B** (W/z inconsistency) present, but **not** sufficient alone to explain the RMSE gap (see §4).

---

## 3. Objective trajectories (Q3)

Traced PIRLS (same updates as package, logging true \(L\)):

- **0 / 40** iterations increased the true penalized Bernoulli objective (r=2 and r=3).
- Train objectives of finished ALS vs L-BFGS are **close**, while test RMSE is not.

| rank | method | train obj | test RMSE(η) | max\|η\| |
|------|--------|----------:|-------------:|---------:|
| 2 | ALS | 285.7 | **1.83** | 5.27 |
| 2 | LBFGS | 284.1 | **0.54** | 3.60 |
| 3 | ALS | 276.8 | **1.52** | 7.93 |
| 3 | LBFGS | 272.8 | **0.80** | 3.80 |

Plots: `q3_traj_r{2,3}.png` (objective, deviance, max\|η\|, median W).

**Interpretation:** ALS is not “failing to decrease \(L\)” on this path. It reaches a **different region** of the nonconvex landscape: similar (sometimes even better) train \(L\), much wilder η, much worse test RMSE.

---

## 4. Weight / probability diagnostics (Q5)

From `q5_separation.csv` (λ=1, shared init):

| rank | method | min p | max p | prop p&lt;1e-6 | prop p&gt;1−1e-6 |
|------|--------|------:|------:|--------------:|-----------------:|
| 2 | ALS | 0.006 | 0.995 | 0 | 0 |
| 2 | LBFGS | 0.027 | 0.943 | 0 | 0 |
| 3 | ALS | 0.0004 | 0.999 | 0 | 0 |
| 3 | LBFGS | 0.022 | 0.949 | 0 | 0 |

Hard separation counts are zero under the 1e-6 rule, but ALS is **systematically more extreme** (larger max\|η\|, p closer to 0/1). Rank 3 amplifies this.

Multi-init (Q10) shows **catastrophic** ALS failures at r=3 (2/5 seeds): max\|η\| up to hundreds, objective \(10^4\)–\(10^5\), soft-damping abort at PIRLS≈12.

**Decision:** **F** (quasi-separation / extreme probabilities) supported, especially as rank grows.

---

## 5. Conditioning (Q6)

Core-1 condition numbers logged per PIRLS iteration (`q6_cond_r*.csv`). Bernoulli weighted Gram + λP remains factorizable with ridge in these runs; no systematic Cholesky collapse on the quick grid.

Conditioning alone does **not** explain the gap better than the extreme-η / path story. Still worth monitoring when W floors bind.

**Decision:** **G** secondary, not primary.

---

## 6. Inner ALS convergence (Q7)

`als_sweeps_per_pirls ∈ {1,2,5,10,20}`, `pirls_maxit=40`, λ=1:

| rank | sweeps | train obj | test RMSE | max\|η\| |
|------|-------:|----------:|----------:|---------:|
| 2 | 1–20 | ≈285.7–286 | ≈1.69–1.85 | ≈5.0–5.3 |
| 2 | LBFGS | **284.1** | **0.54** | 3.6 |
| 3 | 1 | 278.6 | 1.27 | **11.1** |
| 3 | 10 | **259.5** | 1.39 | 6.8 |
| 3 | 20 | 262.7 | 1.19 | 6.1 |
| 3 | LBFGS | 272.8 | **0.80** | 3.8 |

More inner sweeps can **lower train \(L\) below L-BFGS** (r=3, sweeps=10) while **test RMSE stays worse**. So “inexact PIRLS” is **not** the main explanation of underperformance on test.

**Decision:** **C** rejected as primary cause of the RMSE gap.

---

## 7. Outer PIRLS convergence (Q8)

Varying `pirls_maxit` with 10 inner sweeps does not close the ALS–LBFGS test gap on the quick grid (see `q8_outer_pirls.csv`). Soft damping can terminate early after deviance spikes.

Package currently exposes a single `converged = is.finite(deviance)` — too coarse (noted for a later API fix; not changed here).

**Decision:** **D** secondary (early abort via soft damping can freeze bad states).

---

## 8. Step-halving experiment (Q4)

True-objective step-halving on the outer PIRLS step:

| rank | raw #obj↑ | step-halving rejects | median α |
|------|----------:|---------------------:|---------:|
| 2 | 0 | 0 | 1 |
| 3 | 0 | 0 | 1 |

On this DGP/init, full PIRLS steps already decreased \(L\). **Missing step-halving is not why raw PIRLS “fails to descend” here.**

However, the package’s soft deviance guard is a **different** (and harsher) mechanism: it **aborts** instead of damping. That matters for the catastrophic multi-init failures.

**Decision:** **E** partially — not needed for monotone \(L\) on the main seed, but the current soft abort is the wrong safeguard.

---

## 9. Cross warm-start (Q9)

Intercept preserved when polishing ALS → LBFGS via internal API.

| rank | A (ALS) → A2 (LBFGS) | B (LBFGS) → B2 (ALS) |
|------|----------------------|----------------------|
| 2 | Δobj ≈ −0.03; RMSE **1.85 → 1.85** | Δobj ≈ −5.3; RMSE **0.54 → 0.65** |
| 3 | Δobj ≈ −0.71; RMSE **1.32 → 1.51** | Δobj ≈ −13.7; RMSE **0.80 → 1.19** |

Reading:

1. **LBFGS from ALS barely moves** — ALS is near a stationary point of \(L\), but a **bad-predicting** one (large η).
2. **ALS from LBFGS lowers train \(L\) further but worsens test RMSE and increases max\|η\|** — PIRLS+ALS **leaves** a good L-BFGS region toward overconfident fits.

This is the strongest single diagnostic against “just run ALS longer.”

---

## 10. Multiple initializations (Q10)

5 seeds, λ=1, medians:

| rank | method | med obj | med RMSE_te | med max\|η\| |
|------|--------|--------:|------------:|-------------:|
| 2 | ALS | 285.7 | **1.84** | 5.27 |
| 2 | LBFGS | 284.5 | **0.54** | 3.47 |
| 3 | ALS | 273 (IQR up to **2e4**) | **3.22** | **10.2** |
| 3 | LBFGS | 274 | **0.79** | 3.88 |

ALS–LBFGS gap is **systematic**, not one unlucky `tt_initialize()`. At r=3 ALS has heavy right-tail instability.

**Decision:** **H** (initialization / basin) yes for catastrophic failures; even “good” ALS basins generalize poorly vs LBFGS.

---

## 11. Lambda sensitivity (Q11)

| rank | λ | ALS RMSE_te | LBFGS RMSE_te | ALS max\|η\| |
|------|---|------------:|-------------:|-------------:|
| 2 | 0.1 | 1.97 | 0.71 | 5.3 |
| 2 | 1 | 1.85 | 0.54 | 5.2 |
| 2 | 10 | **0.85** | 1.10* | 3.7 |
| 2 | 100 | **0.82** | 1.10* | 3.5 |
| 3 | 0.1 | **202** (blow-up) | 1.75 | 9.3 |
| 3 | 1 | 1.39 | **0.80** | 6.8 |
| 3 | 10 | 1.01 | 0.86 | 5.1 |
| 3 | 100 | 0.92 | 1.10* | 4.6 |

\*LBFGS at λ≥10 often collapses to near-null fits on this quick grid (max\|η\|≈0.06) — separate LBFGS / scaling issue under strong penalty; still, **ALS clearly stabilizes as λ↑**.

**Decision:** **I** (rank × smoothing × logistic) strongly supported.

---

## 12. Dense / small-d check (Q12)

d=2, k=5, n=300 (`q12_d2_check.csv`): ALS again shows larger max\|η\| / worse η RMSE than LBFGS at moderate rank — consistent with a **Bernoulli+PIRLS path** issue, not only high-d TT compression.

(Full Kronecker logistic P-spline baseline left for `--full` / denser follow-up.)

---

## 13. EDF audit (Q13)

Conditional working EDF vs λ (core-wise, final ALS weights) decreases with λ and stays finite/positive (`q13_edf.csv`). No evidence that bogus EDF caused the optimizer gap. EDF is **not** used as an explanation here.

---

## 14. Root cause (classification)

| Code | Verdict |
|------|---------|
| A objective mismatch | **Ruled out** (Q1) |
| B incorrect PIRLS formulas | **Minor** W/z floor inconsistency |
| C insufficient inner ALS | **Ruled out as primary** (Q7) |
| D insufficient outer PIRLS | **Secondary** (soft abort) |
| E missing step-halving | **Partial** — true-obj damping unused; soft abort wrong |
| F quasi-separation / extreme η | **Primary** |
| G poor conditioning | Secondary |
| H init / local basins | **Primary** (esp. r=3 catastrophics) |
| I rank × smoothing × logit | **Primary** |
| J limitation of ALS-PIRLS path for Bernoulli TT | **Primary** (path leaves good LBFGS regions) |
| K other | Soft deviance abort masquerading as “damping” |

### One-sentence root cause

**Bernoulli TT–PIRLS+ALS and L-BFGS optimize the same \(L\), but PIRLS’s local WLS path drives TT fits toward overconfident / large-η stationary regions (worse as rank↑ / λ↓), while L-BFGS consistently finds better-generalizing points; more ALS sweeps can even improve train \(L\) while harming test RMSE.**

---

## 15. Recommended minimal fixes

**Applied** in package (see `BERNOULLI_PIRLS_STABILIZATION.md`): (1) true-objective step-halving; (2) consistent \(W/z\); richer `convergence`. Gate = Outcome B. Package default: Bernoulli → LBFGS via `optimizer="auto"`.

Still optional / not forced as headline:

3. Bernoulli-specific η **diagnostic** logging (done) vs hard η cap (**not** first-line; deferred).
4. Hybrid ALS→LBFGS (**experimental**; did not close gap in gate).
5. Dedicated study of propensity for extreme-logit solutions as \(r\uparrow,\lambda\downarrow\) (**NEXT** scientific angle — §−1).

Do **not** “fix” Bernoulli rhetorically by switching to LBFGS without documenting the PIRLS × rank × λ interaction.

---

## 16. Answer to the final question

> Why does Bernoulli TT-P-spline fitting with PIRLS+ALS currently underperform global L-BFGS, particularly as TT rank increases, and what is the smallest statistically and numerically justified correction?

Because **the objectives match**, but **PIRLS+ALS follows a path that favors extreme linear predictors** in a nonconvex TT logistic landscape. Higher TT rank adds capacity for those extremes; weak λ removes the only brake. L-BFGS stays in smoother basins with better test η. The smallest justified correction is **true-objective (and η-aware) outer damping / step-halving plus consistent W/z weights** — not more inner ALS sweeps, and not pretending the models differ.
