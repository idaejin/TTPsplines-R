# Protocol: fixed-λ ALS C++ path (global \(P_k^{\mathrm{full}}\))

**Package:** `TTPsplines` (`ttpsplines-pkg`)  
**Updated:** 2026-08-17

## Architecture (DECISION)

```text
outer λ (R: Sobol / BO / global-GCV / cGCV)
  +
fixed-λ ALS (C++, global P_k^full, parity vs R)
```

Do **not** port outer cGCV / global-GCV search to C++ first.

## P0 — Documentation — **PASS** (2026-08-17)

- README / `?tt_control` / `?ttps` / scalability vignette: ALS/PIRLS sweeps stay in **R**.
- `backend = "Rcpp"` warns; kernels only (Gram / \(P^{\mathrm{full}}\)).
- `tt_cgcv_fit_cpp` / `tt_glm_pirls_cgcv_cpp` marked **legacy own-margin**.
- Commit: `b27e12a`.

## P1 — Single core update — **PASS** (2026-08-17)

Gate: `tests/testthat/test-als-core-update-global.R` — **352 PASS** (R vs Rcpp).

### Scope

| In | Out |
|----|-----|
| Gaussian | GLM / PIRLS |
| Numeric fixed \(\boldsymbol\lambda\) | cGCV |
| `null_space = "joint"` | profiled NSP |
| One margin \(k\) | Full sweep / multi-sweep fitter |
| Global \(P_k^{\mathrm{full}}\) | Own-margin surrogate |

### API (internal)

- C++: `tt_als_core_update_global_cpp(...)`
- R: `tt_als_core_update_global(..., backend = "R"|"Rcpp")`

### Gate quantities (R vs C++)

\[
S=Z_k^\top W Z_k,\quad
b=Z_k^\top W y_c,\quad
P_k^{\mathrm{full}},\quad
g_k^{\mathrm{new}},\quad
\mathrm{RSS},\quad
J_\lambda,\quad
L_{\mathrm{pen}}=Q
\]

Relative objective gate on well-conditioned problems:

\[
\frac{|L_{\mathrm{Cpp}}-L_R|}{1+|L_R|}<10^{-8}
\]

(also matrix / \(g\) absolute tolerances in `tests/testthat/test-als-core-update-global.R`).

### Cases

- \(d \in \{2,3,5\}\), rank \(\in \{1,2,3\}\)
- isotropic and anisotropic \(\lambda\)
- edge and interior cores
- extreme \(\lambda\) (near-singular stress)
- observation weights

### Next (not P1)

- **P2** one full fixed-λ sweep (same core order as R)
- **P3** multi-sweep fitter + warm start / convergence
- **P4** outer cGCV calling C++ fixed-λ (optional)

## Benchmark note

Before investing in P2–P3, profile interfaces / Gram / \(P^{\mathrm{full}}\) / solve /
R↔C++ overhead on \(d\in\{3,5,10\}\), \(n\in\{500,5000,50000\}\), \(r\in\{2,5,10\}\).
