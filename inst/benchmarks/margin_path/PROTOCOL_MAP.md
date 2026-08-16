# Protocol — Margin Activity Path

**Name:** **Margin Activity Path** (do **not** abbreviate as "MAP"; that collides with *maximum a posteriori*).

**Status:** IDEA → example + exported helper `tt_margin_activity_path()`  
**Date:** 2026-08-17  
**Product lock:** KEEP `cGCV` (unchanged). Margin Activity Path is a **margin screening** layer; final smoothing still uses `cGCV` on the selected subset.  
**Code:** `R/margin_activity_path.R` · demo `inst/benchmarks/margin_path/run_map_example.R` · vignette `vignettes/margin-activity-path.Rmd`  
**Not claimed for Paper-1 SC** as a core product result yet.

---

## Research question

> In moderate / high \(d\), can a cheap **isotropic-\(\lambda\) activity path** rank margins so that a nested **CV** cut recovers the active set and improves holdout risk relative to full-\(d\) TT + `cGCV`?

**Primary object:** selection of a margin subset \(\hat{\mathcal A}\subset\{1,\ldots,d\}\).  
**Secondary:** predictive MSE after refit.  
**Not** the first question: group-lasso on TT cores, joint TT-gGCV, or replacing `cGCV`.

---

## 1. Setup

Data \((y_i,x_i)\), \(x_i\in\mathbb{R}^d\), \(i=1,\ldots,n\).  
TT–P-spline of rank \(r\) and marginal resolution \(k\), with margin-wise roughness weights \(\boldsymbol\lambda\).

**Screening fits** use a common isotropic path
\[
\lambda^{(\ell)}\in\Lambda,\qquad \ell=1,\ldots,L
\]
(e.g. \(\Lambda=\{10^{2},10^{1.5},\ldots,10^{-1}\}\)) and fit
\[
\hat f^{(\ell)}
=
\mathrm{TT}\bigl(y,X;\,r,k,\lambda^{(\ell)}\bigr)
\]
on all \(d\) columns (fixed \(r,k\); product `ttps`).

**Final fits** use anisotropic / conditional GCV:
\[
\lambda=\texttt{"cGCV"}
\]
on the chosen column subset (same \(r,k\) unless noted).

---

## 2. Partial effect and activity

Fix a fitted surface \(\hat f\) and a reference point \(x_{-j}\) (default: **coordinate-wise median** of the training covariates).

**1D partial curve** on a grid \(g_1,\ldots,g_G\) in the central range of \(x_j\) (e.g. 2.5%–97.5% quantiles for plots; for scores a fixed interior grid such as \(\{0.05,\ldots,0.95\}\) on \([0,1]\) covariates is acceptable if all margins share the same scale):
\[
\hat f_j(t)
=
\hat f\bigl(t;\,x_{-j}\bigr).
\]

**Partial range (activity at one \(\lambda\)):**
\[
A_j(\lambda)
=
\max_t \hat f_j(t)
-
\min_t \hat f_j(t).
\]

**Margin score (path aggregate):**
\[
S_j
=
\frac1L\sum_{\ell=1}^L A_j\bigl(\lambda^{(\ell)}\bigr).
\]

**Interpretation:** \(S_j\) large ⇒ the fitted surface varies materially when margin \(j\) moves, across a range of isotropic smoothness levels. Null margins stay near \(A_j\approx 0\).

**Invariance note:** \(A_j\) is a functional of \(\hat f\), not of TT gauges \(\|G_j\|_F\). Prefer this over core-norm penalties for screening.

---

## 3. Ranking and nested models

Order margins by decreasing score:
\[
S_{(1)}\ge S_{(2)}\ge\cdots\ge S_{(d)},
\]
with corresponding index permutation \(\pi\).

For size \(m\in\{0,1,\ldots,d\}\) define the nested active set
\[
\mathcal A_m
=
\begin{cases}
\emptyset & m=0,\\
\{\pi(1),\ldots,\pi(m)\} & m\ge 1.
\end{cases}
\]

Model \(M_m\): TT–P-spline (or intercept if \(m=0\); univariate smoother if \(m=1\), because product `ttps` requires \(d\ge 2\)) on columns \(\mathcal A_m\), with \(\lambda=\texttt{cGCV}\) inside each training fold.

---

## 4. Formal selection rule (primary)

### 4.1 Cross-validation

Partition \(\{1,\ldots,n\}\) into \(K\) folds (default \(K=5\)). For each \(m\),
\[
\mathrm{CV}(m)
=
\frac1K\sum_{k=1}^K
\frac1{|F_k|}
\sum_{i\in F_k}
\bigl(y_i-\hat f^{(-k)}_{\mathcal A_m}(x_i)\bigr)^2,
\]
and let \(\widehat{\mathrm{SE}}(m)\) be the usual fold-wise standard error of the mean.

\[
m_{\min}
=
\arg\min_{m\in\{0,\ldots,d\}} \mathrm{CV}(m).
\]

### 4.2 One-standard-error rule (default)

\[
\hat m
=
\min\Bigl\{
m:\ 
\mathrm{CV}(m)
\le
\mathrm{CV}(m_{\min})+\widehat{\mathrm{SE}}(m_{\min})
\Bigr\}.
\]

**Selected set:** \(\hat{\mathcal A}=\mathcal A_{\hat m}\).

**Rationale:** among nested models near the CV minimum, prefer the **smallest** margin set (Occam / stability), analogous to `cv.glmnet` / `lambda.1se`.

### 4.3 Plain min-CV (optional)

\[
\hat m=m_{\min}.
\]
Use when the path is smooth and under-selection risk dominates.

---

## 5. Diagnostic (not primary)

**Largest-gap rule** on sorted scores:
\[
m_{\mathrm{gap}}
=
\arg\max_{j=1,\ldots,d}
\bigl(S_{(j)}-S_{(j+1)}\bigr)
\quad\text{with }S_{(d+1)}:=0.
\]
Useful for plots and for cases with a clear elbow. **Do not** report as the formal selector when gaps are ambiguous; prefer §4.

Heuristic thresholds such as \(S_j\ge c\cdot\max_j S_j\) (e.g. \(c=0.3\)) are **exploratory only**.

---

## 6. Final estimator and baselines

1. **Margin Activity Path + CV:** refit `ttps` on \(\hat{\mathcal A}\) with \(\lambda=\texttt{cGCV}\) (full training sample).  
2. **Full TT:** `ttps` on all \(d\) margins with `cGCV`.  
3. **Linear lasso (optional baseline):** `glmnet` / `cv.glmnet` on raw \(X\) (\(\alpha=1\)); selection \(\hat\beta_j\neq 0\). Compares *linear* selection+prediction, not a TT method.

Report train / holdout MSE; if a known active set \(\mathcal A^\star\) exists (simulations), also TPR / FPR of \(\hat{\mathcal A}\) (and of lasso).

---

## 7. Algorithm (summary)

```text
Input: y, X (n x d), path Lambda, rank r, k, K folds, rule in {1SE, minCV}
1. For each lambda in Lambda:
     fit isotropic TT on all d margins
     for j = 1..d: A[j,lambda] <- partial_range(fit, j)
2. S[j] <- mean_lambda A[j,lambda];  pi <- order(S, decreasing)
3. For m = 0..d:
     CV[m], SE[m] <- K-fold risk of TT+cGCV on columns pi[1:m]
4. m_hat <- 1-SE(CV) or argmin(CV)
5. Refit TT+cGCV on pi[1:m_hat]; compare to full-d TT+cGCV (+ optional lasso)
Output: A_hat, scores S, CV path, holdout table, figures
```

---

## 8. Implementation map

| Piece | Location |
|-------|----------|
| Driver / demo | `inst/benchmarks/margin_path/run_map_example.R` |
| Package API | `tt_margin_activity_path()` in `R/margin_activity_path.R` |
| Vignette | `vignettes/margin-activity-path.Rmd` |
| This protocol | `inst/benchmarks/margin_path/PROTOCOL_MAP.md` |
| Tables | `results/{activity_path,margin_scores,cv_path,holdout_comparison,lasso_coefs}.csv` |
| Figures | `results/figures/01_*.png` … `05_map_summary.png`, `04b_lasso_coefs.png` |

**Defaults in demo:** \(n=500\), \(d=10\), \(n_{\mathrm{signal}}=5\), \(r=2\), \(k=5\), \(K=5\), `use_one_se=TRUE`.

**Edge cases:** \(m=0\) intercept; \(m=1\) `smooth.spline` GCV (TT API needs \(d\ge 2\)).

---

## 9. Epistemic status

| Claim | Status |
|-------|--------|
| Partial-range path separates clear additive signal vs null margins in the demo DGP | PRELIMINARY (single seed) |
| CV / 1-SE on nested top-\(m\) is the formal cut | DECISION for this protocol |
| Gap / \(c\cdot\max S\) as formal selectors | REJECTED for reporting (diagnostic only) |
| Margin Activity Path replaces `cGCV` or enters Paper-1 SC as product claim | REJECTED / out of scope |
| Group-lasso on \(\|G_j\|_F\) | IDEA parked (gauge non-invariance); prefer functional \(A_j\) |
| Multi-seed MC, correlated designs, interactions-only nulls | NEXT — **started:** `EVIDENCE_REPORT.md` + `run_evidence_mc.R` (B=20 × 2 DGPs) |

See **`EVIDENCE_REPORT.md`** for the current evidence synthesis (2026-08-17).

---

## 10. Open questions

1. Sensitivity of \(S_j\) to reference \(x_{-j}\) (median vs mean vs random profiles).  
2. Path \(\Lambda\): how coarse can it be?  
3. Correlated nulls / weak signal: does ranking stay consistent?  
4. After Margin Activity Path, should rank \(r\) be re-chosen by `tt_rank_select` on \(\hat{\mathcal A}\)?  
5. Fair nonparametric baseline: group-lasso on marginal B-spline bases (not linear `glmnet`).

---

## 11. Relation to other tracks

| Track | Relation |
|-------|----------|
| KEEP `cGCV` | Unchanged; Margin Activity Path only subsets columns |
| Global TT-gGCV / Sobol | Orthogonal; do not reopen product λ |
| Linear lasso in demo | External baseline only |
| TT group-lasso on cores | Not this protocol |

---

## Citation sketch (methods paragraph)

> We screen margins by a *margin activity path*: for an isotropic \(\lambda\)-grid, score each covariate by the mean range of its 1D partial effect with other covariates fixed at their training medians. Margins are ranked by this score. The subset size is chosen by \(K\)-fold CV over nested top-\(m\) TT–P-spline models with conditional GCV, using a one-standard-error rule. The selected margins are then refit with TT + cGCV.
