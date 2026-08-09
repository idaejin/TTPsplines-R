# Bernoulli PIRLS stabilization — FIX 1 + FIX 2 gate

**Status:** implementation + decision gate complete (2026-08-09).  
**Audit (unchanged baseline):** `docs/BERNOULLI_PIRLS_AUDIT.md`  
**Script:** `inst/benchmarks/benchmark_bernoulli_stabilization.R` (`--quick` / `--full`)  
**Results:** `inst/benchmarks/results/bernoulli_stabilization/`  
**Does not overwrite:** `inst/benchmarks/results/bernoulli_audit/`

This note does **not** re-diagnose the audit. It records the minimal fixes, tests, benchmarks, and the A/B/C gate.

---

## 1. Changes implemented

### FIX 1 — consistent Bernoulli working variance

In `glm_working()` (`R/families.R`):

\[
\mu \leftarrow \mathrm{clip}(\mathrm{logistic}(\eta),\varepsilon_\mu),\quad
v_{\mathrm{raw}}=\mu(1-\mu),\quad
v_{\mathrm{work}}=\max(v_{\mathrm{raw}},\,v_{\mathrm{floor}}),
\]

\[
W = v_{\mathrm{work}},\qquad
z = \eta + (y-\mu)/v_{\mathrm{work}}.
\]

Controls (`tt_control()`):

- `binomial_mu_eps = 1e-5`
- `binomial_weight_floor = 1e-4`

(Not claimed optimal — exposed for sensitivity.)

### FIX 2 — true-objective outer step-halving

Removed soft deviance abort (deviance ↑ > 25% → stop).

Bernoulli outer step (when `pirls_step_halving = TRUE`):

1. Run weighted ALS on \((W,z)\) → candidate cores / intercept.
2. Accept step \(\alpha\in\{1,1/2,\ldots\}\) if
   \[
   L(\alpha)\le L_{\mathrm{old}}+\tau,
   \]
   with \(\tau=\max(\texttt{objective\_tol},\,\texttt{tol}\cdot\max(1,|L_{\mathrm{old}}|))\).
3. Blend in **parameter space**:
   \[
   G_k(\alpha)=G_k^{\mathrm{old}}+\alpha(G_k^{\mathrm{cand}}-G_k^{\mathrm{old}}),
   \]
   same for intercept; recompute \(\eta(\alpha)\) by TT contraction.
4. If no \(\alpha\ge\texttt{step\_min}\) improves \(L\):
   - restore last valid iterate;
   - if near a numerical plateau → `reason = "stationary (no improving PIRLS step)"`;
   - else → `pirls = FALSE`, `reason = "PIRLS line search failed"`.

**Acceptance uses the true penalized Bernoulli objective** (same as L-BFGS parity in the audit), not deviance / working RSS alone.

Controls: `pirls_step_halving`, `step_factor=0.5`, `step_min=1/128`, `objective_tol=1e-10`.

### Convergence / history

`fit$convergence` now carries at least:

`overall`, `pirls`, `als`, `reason`, `n_pirls`, `n_als_sweeps`, `n_step_halvings`.

PIRLS history rows (always stored on the R Bernoulli path): objective, nll, penalty, deviance, `max_abs_eta`, min/max probability, weight summaries, ALS sweeps, proposed/accepted step, halvings.

### Explicitly **not** done in this gate

- No hard \(|\eta|\) cap.
- No Adam / Keras.
- No cFS / cREML / penalty change.
- Default optimizer remains ALS (hybrid is experimental only).

### Outcome-B follow-up (experimental)

`optimizer = "hybrid"`: PIRLS+ALS warm-start → short L-BFGS polish (`hybrid_lbfgs_maxit`, default 50) on the **same** penalized objective.

---

## 2. Mathematical justification

PIRLS solves a sequence of **weighted least-squares** TT-ALS problems that approximate Newton steps for the GLM. Those directions need not be descent directions for the true nonlinear objective

\[
L(G,\alpha)=-\ell_{\mathrm{Bern}}(G,\alpha)+\tfrac12\sum_k\lambda_k g_k^\top P_k g_k.
\]

Consistent \((W,z)\) removes a small first-order inconsistency when \(v_{\mathrm{raw}}<v_{\mathrm{floor}}\).  
True-objective step-halving enforces monotone (within tolerance) outer progress in \(L\) along a **parameter-space** segment between consecutive PIRLS iterates.

Neither change alters the statistical model; both alter only the **numerical path**.

---

## 3. W/z parity test

`tests/testthat/test-bernoulli-stabilization.R`:

- Extreme \(\eta\) force flooring.
- Asserts `weight == var_work` and `z == eta + (y-mu)/var_work`.

**PASS.**

---

## 4. Objective monotonicity test

With `pirls_step_halving = TRUE`, accepted outer history satisfies

\[
L_{t+1}\le L_t + \texttt{tolerance}.
\]

Also checks \(\eta\) recomputed from cores+intercept equals `linear.predictors`.

**PASS.**

---

## 5. Gauge / interpolation handling

Current R ALS **does not** re-orthogonalize / canonicalize cores between outer iterates. Old and candidate cores lie on a continuous ALS path, so parameter-space blending is statistically meaningful along that path.

`fit$convergence$gauge_note` documents: if future ALS adds gauge fixing, align representations **before** blending. Gauge non-identifiability is **not** hidden.

No alternative line-search (e.g. η-only) was used as a substitute for TT parameters.

---

## 6. Original benchmark (Exp1)

Full mode: \(n=800\), \(d=3\), \(k=6\), \(\lambda=1\), shared `tt_initialize(seed=123)`.

| rank | method | train \(L\) | test RMSE(\(\eta\)) | max\|\(\eta\)\| | PIRLS | halvings | catastrophic\* |
|-----:|--------|------------:|--------------------:|---------------:|------:|---------:|:--------------:|
| 2 | ALS_old | 483.8 | 0.935 | 5.22 | 40 | 0 | no |
| 2 | ALS_new | 484.1 | 0.933 | 5.48 | 40 | 0 | no |
| 2 | LBFGS | 452.2 | **0.510** | 4.43 | — | 0 | no |
| 3 | ALS_old | 463.8 | 1.368 | **20.13** | 40 | 0 | yes (blowup) |
| 3 | ALS_new | 469.9 | 1.139 | 7.81 | 4 | 8 | yes (line search) |
| 3 | LBFGS | 438.5 | **0.738** | 5.90 | — | 0 | no |

\*Catastrophic (diagnostic only): non-finite \(L\), or max\|\(\eta\)\| > 20, or early `PIRLS line search failed`. **Not** used to alter fits.

**Reading:** FIX1+FIX2 do **not** close the ALS–LBFGS test-RMSE gap. At rank 3 they reduce max\|\(\eta\)\| vs old ALS but often stop via line-search failure after a few outer steps.

---

## 7. Multi-init (Exp2) — 20 seeds

Median (IQR) of test RMSE(\(\eta\)), max\|\(\eta\)\|, time; catastrophic / blowup / line-search rates.

| rank | method | med RMSE | IQR RMSE | med max\|\(\eta\)\| | med \(t\) (s) | blowup | line-search fail |
|-----:|--------|---------:|---------:|-------------------:|--------------:|-------:|-----------------:|
| 2 | ALS_old | 0.935 | [0.933, 0.936] | 5.23 | 0.63 | 0% | 0% |
| 2 | ALS_new | 0.933 | [0.931, 0.934] | 5.52 | 0.68 | 0% | 0% |
| 2 | LBFGS | **0.498** | [0.495, 0.501] | 4.13 | 0.55 | 0% | — |
| 3 | ALS_old | 1.273 | [1.272, 1.295] | 8.70 | 1.21 | **20%** | 0% |
| 3 | ALS_new | 1.139 | [1.133, 1.148] | 7.70 | 0.18 | 0% | **90%** |
| 3 | LBFGS | **0.731** | [0.721, 0.740] | 5.37 | 0.87 | 0% | — |

Files: `exp2_multiinit.csv`, `exp2_multiinit_summary.csv`, `exp2_failure_modes.csv`.

---

## 8. Rank × λ (Exp3)

Ranks \(2,3,4\); \(\lambda\in\{0.1,1,10,100\}\); shared init.

- **Blowups** (max\|\(\eta\)\| > 20 or non-finite): ALS_old rises with rank (0% / 50% / 75%); ALS_new **0%** blowups in this grid; LBFGS 0% / 0% / 25%.
- ALS_new often improves RMSE vs ALS_old at fixed \(\lambda\), but **LBFGS remains clearly better** at moderate \(\lambda\) (e.g. \(\lambda=1\)).
- Strong **rank × λ** interaction remains (small \(\lambda\), large rank → harder landscape).

File: `exp3_rank_lambda.csv`.

---

## 9. Catastrophic failure rates

Pre-specified diagnostic definition (see §6). Separating modes:

| | ALS_old | ALS_new | LBFGS |
|--|--------:|--------:|------:|
| Rank-2 blowup (Exp2) | 0% | 0% | 0% |
| Rank-3 blowup (Exp2) | 20% | 0% | 0% |
| Rank-3 line-search fail (Exp2) | 0% | 90% | — |

**Interpretation:** FIX1+FIX2 reduce **η-blowups** relative to legacy soft-abort ALS, but replace many trajectories with **early line-search stalls** that still leave ALS in a worse-predicting region than LBFGS.

---

## 10. Computational cost

On the full Exp2 protocol (n=800, λ=1, 20 inits):

- ALS_new ≈ ALS_old at rank 2 (~0.7 s).
- At rank 3, ALS_new is *faster* (~0.18 s) because line search often aborts early — not because it reaches LBFGS quality.
- Cold LBFGS ≈ 0.55–0.90 s and dominates statistically.

---

## 11. Decision gate

| Outcome | Criterion | Result |
|---------|-----------|--------|
| **A** | Eliminate catastrophic failures **and** substantially close ALS–LBFGS RMSE gap | **FAIL** |
| **B** | Numerical path more controlled, but ALS still reaches worse-predicting stationary regions | **PASS** |
| **C** | FIX1+FIX2 worsen / destabilize | **FAIL** (mild improvement vs old ALS; model unchanged) |

**Gate result: OUTCOME B.**

---

## 12. Hybrid follow-up (Exp4) — only because Outcome B

`optimizer = "hybrid"`: ALS/PIRLS → L-BFGS polish (`hybrid_lbfgs_maxit = 50`).

| rank | method | med RMSE | med max\|\(\eta\)\| | med \(t\) |
|-----:|--------|---------:|-------------------:|----------:|
| 2 | ALS_new | 0.933 | 5.52 | 0.69 |
| 2 | hybrid | 0.933 | 5.47 | 0.78 |
| 2 | LBFGS | **0.498** | 4.13 | 0.58 |
| 3 | ALS_new | 1.139 | 7.70 | 0.19 |
| 3 | hybrid | 1.147 | 8.03 | 0.36 |
| 3 | LBFGS | **0.731** | 5.37 | 0.90 |

Probe (same init, rank 2, λ=1): hybrid polish with `maxit ∈ {50,100,300}` stays at \(L\approx 484\), RMSE≈0.94; **cold LBFGS from the same init** reaches \(L\approx 451\), RMSE≈0.49.

**Conclusion:** short (or even long) L-BFGS polish **does not escape** the ALS basin. Hybrid is **not** a statistically meaningful shortcut to cold LBFGS on this Bernoulli fixed-λ protocol.

Files: `exp4_hybrid.csv`, `exp4_hybrid_summary.csv`.

---

## 13. Recommendation for the package

1. **Default optimizer is `auto`:** **LBFGS for Bernoulli**, **ALS** for Gaussian/Poisson.
2. **Ship FIX1 + FIX2** for Bernoulli PIRLS when `optimizer = "ALS"` is requested explicitly (consistent \(W/z\), true-objective step-halving, richer `convergence` / history).
3. **Document** that Bernoulli ALS remains available for research / cGCV alignment studies, but is not the recommended fixed-λ predictor.
4. Keep `optimizer = "hybrid"` as **experimental**; do **not** advertise it as closing the Bernoulli gap (evidence against).
5. Do **not** add a hard η cap yet; keep logging `max_abs_eta`.
6. Revisit Bernoulli **cGCV** with the LBFGS default path (this gate used fixed λ only).

---

## 14. Recommendation for Paper 1 (lab narrative)

### Framing — **DECISION** (do not overclaim)

Do **not** sell this as discovering a new pathology of Bernoulli TT models.

**ESTABLISHED (classical logistic / IRLS):** extreme fitted probabilities can degenerate IRLS weights \(W_i=p_i(1-p_i)\to 0\) (complete/quasi-complete separation, near-separation, overconfident fits). Standard remedies include step-halving, trust regions, stronger regularization, Firth / bias-reduced likelihood.

**PRELIMINARY (TT-P-spline specific):** that known mechanism **interacts** with TT capacity and P-spline smoothing:

\[
\underbrace{\text{extreme }p\text{ / IRLS weight collapse}}_{\text{known logistic phenomenon}}
\;+\;
\underbrace{\text{nonconvex TT factorization }\times\, r\times\lambda}_{\text{TTPsplines structure}}.
\]

Working narrative sentence (**POSITIONING**):

> As is well known for logistic regression, extreme fitted probabilities can lead to degenerating IRLS weights. In the TT-P-spline setting, this phenomenon interacts with the TT rank and the strength of smoothness regularization, requiring appropriate stabilization of the PIRLS updates.

The scientifically interesting claim to pursue (**HYPOTHESIS**, not yet a paper headline):

\[
r\uparrow,\quad \lambda\downarrow
\quad\Longrightarrow\quad
\text{greater propensity for extreme-logit / near-separated PIRLS solutions}.
\]

That is: **not** “IRLS fails with extreme \(p\)” (known), but **how rank and smoothing modulate that failure mode** inside logistic TT-P-splines.

### Optimizer path (audit + gate)

- ALS/PIRLS and L-BFGS evaluate the **same** penalized Bernoulli objective (**ESTABLISHED**).
- L-BFGS does **not** magically remove separation risk; it follows a **different trajectory** in the nonconvex TT-core landscape and can land in better-generalizing basins even when ALS train \(L\) is similar or slightly better (**PRELIMINARY**).
- FIX1 (consistent \(W/z\)) + FIX2 (true-objective step-halving) are **natural PIRLS stabilizations**, not TT-specific hacks (**DECISION** on interpretation).
- They improve path hygiene but do **not** close the ALS–LBFGS predictive gap alone (gate Outcome B).
- Hybrid ALS→LBFGS polish stays in the ALS basin here — further evidence of path-dependent local structure.

### Publishable threads (do not force one)

- **(C/D)** rank × smoothing × logistic likelihood interaction (main interesting angle).
- Optimizer-dependent basins under a shared objective.
- ALS+cGCV may still matter for **automatic smoothing** — keep separate from fixed-λ Bernoulli prediction.

---

## Final questions

### Do consistent PIRLS working quantities and true-objective step-halving remove the Bernoulli instability of TT-PIRLS+ALS without changing the underlying statistical model?

**No.**  
The statistical model is unchanged, and the path is somewhat better controlled (fewer η-blowups vs legacy soft-abort ALS; monotone \(L\) when steps are accepted). They do **not** remove the ALS–LBFGS test-RMSE gap or the rank×λ fragility.

### Does a short ALS/PIRLS → L-BFGS hybrid provide a statistically meaningful and computationally efficient alternative?

**No (on these benchmarks).**  
Hybrid stays near the ALS solution; cold LBFGS from the same initialization finds a better-predicting basin at comparable or lower cost. Do not treat hybrid as a Bernoulli fix without further evidence.
