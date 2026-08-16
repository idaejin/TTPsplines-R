# Gate: profiled null-space TT (Gaussian)

**Status:** PRELIMINARY prototype gate (not Paper-1 claim).  
**Date:** 2026-08-17

## Claim under test

`null_space = "profiled"` implements the **joint profiled** Gaussian criterion

\[
\min_{\beta_0,\mathcal G}
\bigl\|y-X_0\beta_0-\mu_{\mathrm{TT}}(\mathcal G)\bigr\|^2
+
J_\lambda\bigl(\mu_{\mathrm{TT}}(\mathcal G)\bigr)
\]

via exact profiling \(\widehat\beta_0(\mathcal G)\) and ALS on \(Q_0 y\), \(Q_0 Z_k\).

## Estimators compared

| Mode | Meaning |
|------|---------|
| `joint` | Single TT on full \(\Theta\) (default) |
| `profiled` | \(Q_0\)-orthogonalized / profiled NSP–TT (Gaussian ALS) |

## Minimal scenarios (\(d=2\), \(q=2\), fixed \(\lambda\))

1. **Null only** — affine truth; \(\lambda\) large  
2. **Null + rough** — affine + interaction wiggle  
3. **Rough only** — mean-zero interaction  
4. **Correlated rough** — rough term with a linear trend component

## Pass criteria (mechanical)

- Profiled: \(X_0^\top \widehat\mu_{\mathrm{TT},\perp}\approx 0\) with \(\mu_{\perp}=Q_0\mu_{\mathrm{TT}}\)
- Profiled: in-sample identity \(\widehat y = P_0 y + Q_0\mu_{\mathrm{TT}}\)
- Null-only + large \(\lambda\): \(\|\mu_{\perp}\|\) small; good affine recovery
- No claim yet on cGCV / Poisson / global GCV

## Run

```bash
Rscript inst/benchmarks/null_space/run_gate_profiled.R
```

## Interpretation

> **PROFILED Gaussian — prototype ALS with \(Q_0\) residualization; gate below. Not yet NSP default; not GLM.**  
> Former `sequential` / `separate` offset split was removed (not equivalent to joint NSP).

TODO-SC-NULL remains open until profiled is validated more broadly.
