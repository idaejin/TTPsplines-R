# ALS vs L-BFGS diagnostic findings (Gaussian, n=800, d=3, k=6)

**Script:** `inst/benchmarks/diagnose_als_vs_lbfgs.R`  
**Outputs:** `inst/benchmarks/results/diagnose_als_vs_lbfgs*.csv`

## Fixed λ = 1 — training penalized objective

Definition (shared by ALS and L-BFGS):

\[
\frac12\|y-\alpha-\hat f\|^2 + \frac12\sum_k \lambda_k\, g_k^\top P_k g_k.
\]

| Setting (r=2, shared init) | Objective | RSS | RMSE test |
|----------------------------|----------:|----:|----------:|
| ALS 20 sweeps (old default-ish) | 79.3 | 150.8 | 0.350 |
| ALS 200 sweeps, tol=1e-12 | **35.56** | 68.0 | **0.072** |
| L-BFGS from same init | 37.45 | 68.9 | 0.082 |
| L-BFGS warm-started from ALS | **35.55** | 68.0 | 0.073 |
| ALS warm-started from L-BFGS | 36.23 | 68.2 | 0.075 |

**Conclusion (fixed λ, good init):** the earlier ALS ≪ L-BFGS RMSE gap was largely **early stopping**, not a structural ALS failure. Long ALS reaches (and slightly beats) cold-start L-BFGS on the **same** penalized objective. L-BFGS from the ALS solution only shaves ~0.01 in objective — ALS was already near a good stationary point.

## Multi-init (12 seeds, ALS 80 sweeps)

| Rank | ALS better obj | median obj ALS | median obj LBFGS | median RMSE_te ALS | median RMSE_te LBFGS |
|------|----------------|---------------:|-----------------:|-------------------:|---------------------:|
| 2 | 2/12 | 63.6 | **37.7** | 0.246 | **0.082** |
| 3 | 0/12 | 38.0 | **36.3** | 0.118 | **0.099** |

**Conclusion:** L-BFGS is **more robust across random initializations**. ALS can match L-BFGS from a favorable start + enough sweeps, but often stalls in worse basins under moderate sweep budgets.

## cGCV

| Rank | Method | λ (approx) | RMSE test |
|------|--------|------------|----------:|
| 2 | ALS-cGCV | (0.01, 0.01, 7.0) anisotropic | **0.102** |
| 2 | LBFGS-outer-cGCV | (0.01, 0.01, 0.01) boundary | 0.235 |
| 3 | ALS-cGCV | (0.01, 0.01, 4.9) | 0.118 |
| 3 | LBFGS-outer-cGCV | (0.01, 0.01, 0.01) | 0.113 |

**Conclusion:** ALS-cGCV stays aligned with **conditional** λ updates and finds anisotropic λ; outer L-BFGS+cGCV tends to hit the lower bound uniformly. Structural alignment of core-wise ALS with conditional GCV remains a distinct algorithmic story — separate from the fixed-λ early-stopping issue.

## Package follow-ups

1. Default `max_sweeps` raised (20 → 50) and `tol` tightened (1e-6 → 1e-8).
2. Use `tt_objective(fit, X)` for fair training-objective comparisons.
3. Optimizer benchmarks should report **objective**, not only RMSE, and use adequate ALS sweeps.
4. Narrative candidate: ALS-cGCV as the natural conditional algorithm; L-BFGS as global contrast / multi-start insurance for fixed λ.
