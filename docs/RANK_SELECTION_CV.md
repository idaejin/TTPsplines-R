# TT rank selection (CV + 1-SE) — method note for Paper 1

**Status:** ESTABLISHED in package API (2026-08-10).  
**Provenance:** Factual (implemented); simulation audit PRELIMINARY.

## DECISION — fitting vs selection stay separate

`ttpspline(..., rank = r)` means **exactly** that rank. No automatic rank
inside the fitter. Users who want data-driven \(r\) call:

```r
sel <- tt_rank_select(y, X, ranks = 1:5, rule = "1se", ...)
fit <- tt_rank_refit(sel)
```

## Complexity layers

\[
r = \text{structural capacity},\quad
\lambda = \text{directional smoothness},\quad
\mathrm{EDF} = \text{effective fitted flexibility}.
\]

\[
r \neq \lambda \neq \mathrm{EDF}.
\]

## Procedures

| Symbol | Definition | Role |
|--------|------------|------|
| \(r_{\mathrm{oracle}}\) | \(\arg\min_r \mathrm{RMSE}(\hat f_r, f_{\mathrm{true}})\) | simulation audit only |
| \(r_{\min}\) | \(\arg\min_r \overline{\mathrm{CV}}(r)\) | minimum-CV |
| \(r_{1\mathrm{SE}}\) | smallest \(r\) with \(\overline{\mathrm{CV}}(r)\le \overline{\mathrm{CV}}(r_{\min})+\mathrm{SE}(r_{\min})\) | **recommended** default |

Architecture when `lambda = "cGCV"`: outer CV over \(r\); inner cGCV on
**training folds only** (no validation leakage).

## Wording (avoid overclaim)

Prefer: “The 1-SE rule selects the smallest candidate rank whose
cross-validated error lies within one standard error of the minimum.”

Do **not** call this a formal test of the true tensor rank.

## Out of scope (still)

LRT / bootstrap rank tests, adaptive TT-SVD, non-uniform rank-chain search,
Paper-2 cFS.
