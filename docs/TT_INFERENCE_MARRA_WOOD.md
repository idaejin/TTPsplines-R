# Uncertainty for Tensor-Train P-spline smooths

**Status:** Gate 1/2 implementation note (Gaussian conditional inference).  
**Package:** `TTPsplines` · see also vignette `uncertainty`.

## 1. Statistical motivation

For a penalized TT P-spline with objective

\[
Q(\theta)=-\ell(\theta)+\tfrac12\theta^\top S_\lambda\theta,
\]

Marra–Wood / Wahba / Nychka style inference targets **pointwise** uncertainty for the identifiable smooth

\[
\hat f(x),\qquad\text{not individual TT core entries}.
\]

Literature (mgcv; Marra & Wood 2012) distinguishes:

| Level | Object | v1 |
|---|---|---|
| 1 | \(\operatorname{Var}\{\hat f(x)\mid \hat r,\hat\lambda\}\) | **implemented** |
| 2 | smoothing-parameter uncertainty (`unconditional=TRUE`) | **not implemented** |
| 3 | rank-selection uncertainty | **not implemented** |

v1 intervals are therefore

\[
\mathrm{CI}\bigl[f(x)\bigm|\hat r,\hat\lambda\bigr]
\]

approximate, local-linearization, **pointwise** (not simultaneous bands).

Conceptual positioning (four distinct objects):

\[
\boxed{r\neq\lambda\neq\mathrm{EDF}\neq\text{uncertainty}}
\]

- \(r\): structural TT capacity  
- \(\lambda\): directional smoothness  
- EDF: effective complexity after penalization  
- \(\operatorname{SE}(\hat f)\): conditional inferential uncertainty  

## 2. Bayesian covariance (default)

Let \(J\) be the stacked Jacobian of the TT map \(f(\theta)\) w.r.t. packed cores (same construction as joint linearized EDF), and append an intercept column:

\[
X=\bigl[1\mid J\bigr],\qquad
H=X^\top X,\qquad
S=\mathrm{blkdiag}(0,S_\lambda),\qquad
H_p=H+S+\varepsilon I.
\]

Gaussian approximate Bayesian covariance:

\[
V_B=\hat\sigma^2 H_p^{-1}.
\]

Default for `predict(..., se.fit=TRUE)` and `vcov(fit)`.

Motivation: Bayesian/ penalized-spline covariance typically has better frequentist interval properties for smooth functions than the sandwich form because it accounts more appropriately for smoothing bias (Marra & Wood).

## 3. Frequentist sandwich (secondary)

\[
V_F=\hat\sigma^2 H_p^{-1} H H_p^{-1}.
\]

Available as `vcov(fit, type="frequentist")` and `predict(..., vcov_type="frequentist")` for research / benchmarking. **Not** the default.

## 4. TT gauge issue

TT cores are gauge-redundant:

\[
G_k\mapsto G_k A,\qquad G_{k+1}\mapsto A^{-1}G_{k+1}
\]

leaves the coefficient tensor unchanged. Consequences:

- raw \(H_p\) in packed coordinates is singular / nearly singular along gauge directions;
- core-entry covariance is not scientifically interpretable;
- inference must target \(f(x)\).

## 5. Parameterization used (Option A, light)

**Left-orthogonal QR gauge fix** (sequential) before forming \(J\) and \(H_p\):

1. Unfold each \(G_k\) (\(k=1,\ldots,d-1\)) as \((r_{k-1}p_k)\times r_k\).
2. Thin QR (non-pivoted); absorb \(R\) into \(G_{k+1}\).
3. Sign-normalize \(\mathrm{diag}(R)\).

Inference coordinates: `(intercept, packed left-orthogonal cores)`.

Dimension note: factorization size is \(N_{\mathrm{TT}}+1\) (stored parameters + intercept). Intrinsic manifold dimension is \(N_{\mathrm{TT}}-\sum r_k^2\); gauge null directions are controlled by the penalty plus a tiny ridge. Prediction SEs are validated to be **gauge-invariant**.

## 6. Prediction Jacobian

For newdata rows,

\[
J_x=\frac{\partial\eta(x)}{\partial\vartheta}\Big|_{\hat\vartheta}
=\bigl[1\mid \partial f/\partial\theta\bigr],
\]

built from the same conditional designs as EDF (`tt_design_core` / `tt_stacked_jacobian`). Then

\[
\widehat{\mathrm{Var}}\{\hat\eta(x)\}=J_x V_\vartheta J_x^\top,\qquad
\mathrm{SE}=\sqrt{\widehat{\mathrm{Var}}}.
\]

Implementation uses Cholesky solves of \(H_p\) and computes **diagonal-only** prediction variance (`full_cov=FALSE`). No \(m\times m\) prediction covariance is materialized.

## 7. Conditional interpretation

Always:

```text
Conditional on TT rank:   yes
Conditional on lambda:    yes
Smoothing uncertainty:    not included
Rank uncertainty:         not included
```

If `lambda="cGCV"`, \(\hat\lambda\) is still treated as **fixed** in \(V_B/V_F\).  
`unconditional=TRUE` errors clearly — it must not silently return the conditional covariance.

## 8. Gaussian derivation / scale

Objective matches ALS / L-BFGS Gaussian path (½ RSS + penalty). Gauss–Newton data Hessian \(H=X^\top X\).

Scale estimator (documented on `fit$inference$scale_estimator`):

\[
\hat\sigma^2=\frac{\mathrm{RSS}}{n-\mathrm{edf}_{\mathrm{aug}}},\qquad
\mathrm{edf}_{\mathrm{aug}}=\mathrm{tr}(H_p^{-1}H)
\]

(including the unpenalized intercept). Compatible with the joint linearized EDF infrastructure (same \(J\), same \(S_\lambda\), gauge-fixed cores).

## 9. GLM extension

**Not enabled in Gate 1/2.** Planned: Fisher / observed information on the link scale with PIRLS weights; intervals built on \(\eta\) then inverse-link transformed. Do not form \(\hat p\pm z\,\mathrm{SE}_p\) as default.

## 10. Coverage validation

Pilot (`inst/benchmarks/benchmark_inference_coverage.R`; 5 seeds; \(n=800\); pointwise 95%):

| DGP | Covariance | Nominal | Empirical (mean) | Mean width | RMSE |
|---|---|---:|---:|---:|---:|
| near_null (affine) | Bayesian | 0.95 | 0.900 | 0.246 | 0.075 |
| near_null | Frequentist | 0.95 | 0.872 | 0.233 | 0.075 |
| Ishigami (r=5, k=8, λ=0.3) | Bayesian | 0.95 | 0.654 | 1.13 | 0.47 |
| Ishigami | Frequentist | 0.95 | 0.627 | 1.05 | 0.47 |
| Friedman (r=5, k=8, λ=0.5) | Bayesian | 0.95 | 0.995 | 8.62 | 2.05 |
| Friedman | Frequentist | 0.95 | 0.990 | 7.62 | 2.05 |

Times (mean): fit ~0.01–0.93s; cov setup ~0.002–0.46s; pred-SE ~0.001–0.08s (n_test=200).

**Findings (PRELIMINARY):**

1. When TT rank is **grossly insufficient**, RMSE dominates and coverage collapses — not a covariance bug; do not claim coverage for under-ranked fits.
2. With richer rank, Ishigami coverage improves toward nominal (single-seed check ~0.98 at rank 6 / \(k=10\) / \(\lambda=0.1\)).
3. Near-null / heavy smoothing can show mild undercoverage (Marra–Wood-type behaviour).
4. Bayesian vs frequentist: similar coverage; Bayesian slightly wider on average.

Coverage is **not ESTABLISHED**. Larger factorial study is NEXT (not executed here).

## 11. Limitations

- Conditional on \(\hat r,\hat\lambda\) only.
- No smoothing-parameter uncertainty (cGCV / REML delta correction).
- No rank-selection uncertainty.
- Pointwise intervals only (not simultaneous bands).
- Approximate Bayesian covariance / local linearization — not an exact posterior.
- Gaussian first; GLM gates deferred.
- Dense factorization: limited by `edf_max_npar`.

## 12. Future: smoothing-parameter uncertainty

Possible delta route with \(\rho=\log\lambda\):

\[
V_C\approx V_B+J_\rho V_\rho J_\rho^\top.
\]

\(V_\rho\) must match the smoothing criterion. REML/LAML literature exists; **do not** transplant formulas blindly to cGCV. Paper 2 / future inference.

## 13. Future: rank-selection uncertainty

Not in v1. Candidate: bootstrap entire pipeline (resample → rank select → λ → fit → \(f(x)\)).

## References (conceptual)

- Wahba; Nychka — Bayesian confidence intervals for smoothers.
- Marra & Wood (2012) — coverage of GAM intervals; Bayesian vs frequentist covariance.
- Wood — mgcv `vcov` / `unconditional` machinery (adapted here to TT gauge).
