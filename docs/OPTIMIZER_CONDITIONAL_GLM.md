# Conditional GLM solvers (Damped-Newton-ALS / LBFGS-ALS)

Status: **gate completed** (2026-08-09). Fixed \(\lambda\) only; **no cGCV**.

## 1. Conditional objective

With other TT cores fixed, the predictor is affine in the vectorized core \(g_k=\mathrm{vec}(G_k)\):

\[
\eta = o_k + X_k g_k,
\qquad
Q_k(g_k)
=
-\ell(y;o_k+X_kg_k)
+\tfrac{\lambda_k}{2}g_k^\top P_k g_k.
\]

For Bernoulli-logit, \(\mu=\mathrm{logit}^{-1}(\eta)\),

\[
\nabla Q_k = X_k^\top(\mu-y)+\lambda_k P_k g_k,
\qquad
H_k = X_k^\top W X_k + \lambda_k P_k,
\quad W_i=\mu_i(1-\mu_i).
\]

Global penalized objective (shared by all solvers at fixed \(\lambda\)):

\[
Q(G,\alpha)
=
-\ell(y;\eta_G)
+\tfrac12\sum_k\lambda_k J_k(G_k).
\]

Implementation: `R/optimizer_conditional.R` reuses TT interfaces / `tt_design_core`,
`tt_glm_penalized_objective`, and sparse `crossprod` / `solve_spd_ridge`.

**Intercept:** fixed at GLM `init_intercept` (same convention as global `LBFGS`).
PIRLS-ALS may refresh the intercept inside its outer loop — that difference is
documented, not “fixed away” by changing the model.

## 2. PIRLS-ALS

Constructs working weights / responses and solves **weighted least-squares**
core problems (optionally with true-objective step-halving). It does **not**
require Armijo on \(Q_k\) at every Newton proposal.

## 3. Damped-Newton-ALS (`optimizer = "Damped-Newton-ALS"`)

Per core: form true \(\nabla Q_k\) and \(H_k\), solve \(H_k\Delta=-\nabla Q_k\),
accept \(\alpha\Delta\) only under **Armijo** on the true \(Q_k\)
(\(c=10^{-4}\), \(\rho=0.5\) by default). Optional explicit `dn_ridge`
(default `0`; never silent).

This is **not** PIRLS under another name: direction may be Newton/IRLS-related
for canonical Bernoulli, but **step control is on the true penalized likelihood**.

## 4. Block-LBFGS-ALS (`optimizer = "LBFGS-ALS"`)

Per core: minimize \(Q_k\) by L-BFGS-B with the **analytical** conditional
gradient; then update the core and continue ALS. Distinct from global `LBFGS`.

## 5. Global LBFGS (`optimizer = "LBFGS"`)

Joint optimization over **all** packed TT cores with the same \(Q\).

| Name | Scope | Objective used for steps |
|------|-------|--------------------------|
| `PIRLS-ALS` | alternating cores | working WLS (+ optional outer true-obj halving) |
| `Damped-Newton-ALS` | alternating cores | true \(Q_k\) + Armijo |
| `LBFGS-ALS` | alternating cores | true \(Q_k\) via L-BFGS |
| `LBFGS` | all cores jointly | true \(Q\) |
| `GD` | all cores jointly | true \(Q\) (diagnostic; not in this gate) |

## 6. Gradient / parity checks

Automated in `tests/testthat/test-conditional-solvers.R`:

- Conditional \(\eta\) / NLL matches global evaluator at shared cores.
- Analytical \(\nabla Q_k\) vs central FD (Bernoulli + Poisson), \(\max_j|\cdot|<10^{-4}\).
- Damped-Newton history monotone in global \(Q\) (numerical tol).
- LBFGS-ALS smoke fit finite with probabilities in \((0,1)\).

**Result:** all checks **PASS**. Gate proceeded.

## 7. Bernoulli gate

Script: `inst/benchmarks/benchmark_conditional_glm.R --full`  
Results: `inst/benchmarks/results/conditional_glm/`  
Setup: same DGP helpers as prior Bernoulli benches; \(n=800\), \(n_{\mathrm{te}}=2500\),
\(k=6\), ranks \(2,3\), \(\lambda=1\), **10 seeds**; methods above (no GD).

### Median test metrics (10 seeds)

| method | rank | median RMSE(\(\eta\)) | median log-loss | median max\|\(\eta\)\| | median runtime | failure rate |
|--------|-----:|----------------------:|----------------:|-----------------------:|---------------:|-------------:|
| PIRLS-ALS | 2 | 0.933 | 0.654 | 5.52 | 0.68s | 0 |
| Damped-Newton-ALS | 2 | 0.818 | 0.640 | 4.29 | 0.09s | 0 |
| LBFGS-ALS | 2 | 0.815 | 0.640 | 4.34 | 0.36s | 0 |
| LBFGS | 2 | 0.496 | 0.603 | 4.10 | 0.55s | 0 |
| PIRLS-ALS | 3 | 1.143 | 0.683 | 7.70 | 0.20s | 0 |
| Damped-Newton-ALS | 3 | 1.162 | 0.679 | 8.17 | 0.18s | 0 |
| LBFGS-ALS | 3 | 0.932 | 0.657 | 5.74 | 0.66s | 0 |
| LBFGS | 3 | 0.727 | 0.631 | 5.22 | 0.86s | 0 |

Shared-init snapshot is consistent with the multi-init ordering
(global LBFGS best; ALS block methods do not close the gap).

## 8. Poisson sanity

Same ranks / \(\lambda=1\) / 10 seeds (no claim of a full Poisson study).

| method | rank | median RMSE(\(\eta\)) | median max\|\(\eta\)\| | median runtime | failure rate |
|--------|-----:|----------------------:|-----------------------:|---------------:|-------------:|
| PIRLS-ALS | 2 | 0.107 | 2.64 | 0.72s | 0 |
| Damped-Newton-ALS | 2 | 0.404 | 2.71 | 0.09s | 0 |
| LBFGS-ALS | 2 | 0.345 | 2.71 | 0.42s | 0 |
| LBFGS | 2 | 0.109 | 2.64 | 0.58s | 0 |
| PIRLS-ALS | 3 | 0.146 | 2.62 | 1.30s | 0 |
| Damped-Newton-ALS | 3 | 0.317 | 2.69 | 0.18s | 0 |
| LBFGS-ALS | 3 | 0.178 | 2.71 | 0.72s | 0 |
| LBFGS | 3 | 0.170 | 2.72 | 0.88s | 0 |

**Sanity:** new conditional solvers are **not** free wins. Where PIRLS already
works (Poisson), PIRLS \(\approx\) global LBFGS and **beats** Damped-Newton /
often beats LBFGS-ALS. Do not replace Poisson defaults with DN/LBFGS-ALS.

## 9. Runtime

On this gate, Damped-Newton-ALS is fast (often early relative-\(Q\) stop) but
lands in weaker basins. LBFGS-ALS is mid-cost. Global LBFGS remains competitive
in wall time for Bernoulli while delivering better test RMSE.

## 10. Interpretation (scenarios A–E)

Pre-registered:

- **A** \(LBFGS\text{-}ALS \approx\) global \(LBFGS \gg\) PIRLS → alternating OK, PIRLS path bad  
- **B** \(LBFGS\text{-}ALS \approx\) PIRLS \(\ll\) global → alternating itself hurts  
- **C** DampedNewton \(\approx\) global → true-obj Newton enough  
- **D** DN \(\approx\) PIRLS but LBFGS-ALS \(\approx\) global → need conditional direct likelihood  
- **E** all core-wise ALS \(\ll\) global LBFGS → Bernoulli favors joint direct likelihood  

**Best match: Scenario E**, with a mild **B/D hybrid nuance**:

- LBFGS-ALS is **somewhat better** than PIRLS-ALS on Bernoulli (esp. rank 3),
  so pure “ALS ≡ PIRLS” is too strong.
- Neither LBFGS-ALS nor Damped-Newton-ALS approaches global LBFGS.
- Damped-Newton \(\approx\) PIRLS (or slightly better at rank 2); **not** Scenario C.

**HYPOTHESIS (retained):** solver choice may depend on likelihood geometry under a
common TT–P-spline representation. This gate supports **Bernoulli → prefer
global direct-likelihood (`LBFGS`)**; **Poisson → keep PIRLS-ALS** as a strong
structure-aware default.

## 11. Recommendation for Paper 1

- Report Bernoulli as a **solver-path** issue under TT, not a new pathology claim.
- Headline Bernoulli optimizer: **global LBFGS** (or other joint direct \(\mathcal L\)).
- Present Damped-Newton-ALS / LBFGS-ALS as **negative / partial** conditional
  experiments: true-objective core control does **not** recover global quality.
- Keep Poisson PIRLS success in the narrative (this sanity check agrees).

## 12. Package API recommendation

Provisional API (unchanged scientific defaults):

```
optimizer = "PIRLS-ALS" | "Damped-Newton-ALS" | "LBFGS-ALS" | "LBFGS" | "GD" | "auto"
```

- `auto`: family-aware v1 defaults — Gaussian → `ALS`, Poisson → `PIRLS-ALS`,
  binomial → `LBFGS`. Transparent via `optimizer_requested` / `optimizer_used` /
  `optimizer_reason`. Do **not** map `auto` to Damped-Newton or LBFGS-ALS.
- Manual overrides always available (`PIRLS-ALS`, `LBFGS`, …).
- `GD`: experimental diagnostic baseline (see `OPTIMIZER_GD.md`).
- `Damped-Newton-ALS`, `LBFGS-ALS`: **EXPERIMENTAL** public options.

### Classification

| method | decision |
|--------|----------|
| `LBFGS` | **KEEP — production candidate** (Bernoulli default via `auto`) |
| `PIRLS-ALS` | **KEEP — production** (Poisson / Gaussian ALS path); Bernoulli: available, not preferred |
| `LBFGS-ALS` | **KEEP — benchmark/reference** (EXPERIMENTAL production) |
| `Damped-Newton-ALS` | **EXPERIMENTAL** (not a Bernoulli closer; weak on Poisson) |
| `GD` | **EXPERIMENTAL** diagnostic (completed earlier) |

### NEXT (not in this task)

1. rank \(\times\) \(\lambda\) validation for survivors  
2. surviving solvers \(\times\) cGCV  
3. scaling study  
4. real-data Poisson traffic application  
5. final package API freeze  
