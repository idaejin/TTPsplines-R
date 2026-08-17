# Fixed-λ ALS velocity (R vs Rcpp) — 2026-08-17

**Machine:** Darwin arm64 · **R:** 4.6.1 · **Commit context:** after P3 (`b97e8a8`)  
**Script:** `run_fixed_lambda_velocity.R`

```bash
Rscript inst/benchmarks/global_gcv/run_fixed_lambda_velocity.R
TT_GGCV_VEL_MODE=QUICK Rscript inst/benchmarks/global_gcv/run_fixed_lambda_velocity.R
TT_GGCV_VEL_OUT=/tmp/tt_vel Rscript inst/benchmarks/global_gcv/run_fixed_lambda_velocity.R
```

## Fit (`tt_als_fit_fixed_global`, \(k=8\), \(\lambda=1\), `max_sweeps=20`)

| n | d | r | med R (s) | med Cpp (s) | speedup | sweeps |
|---|---|---|-----------|-------------|---------|--------|
| 500 | 3 | 2 | 0.057 | 0.030 | 1.90× | 20 |
| 500 | 5 | 2 | 0.135 | 0.077 | 1.75× | 20 |
| 5 000 | 3 | 2 | 0.307 | 0.205 | 1.50× | 20 |
| 5 000 | 5 | 2 | 0.801 | 0.551 | 1.45× | 20 |
| 5 000 | 3 | 5 | 2.68 | 2.45 | 1.10× | 20 |
| 50 000 | 3 | 2 | 1.86 | 1.31 | 1.42× | 13 |
| 5 000 | 10 | 2 | 1.87 | 1.04 | 1.79× | 13 |
| 5 000 | 5 | 5 | 7.48 | 6.87 | 1.09× | 20 |

Objective relative difference R vs Cpp: \(\sim 10^{-15}\).

## Micro (n=5 000, d=5, r=2)

| Component | med R (s) | med Cpp (s) | speedup |
|-----------|-----------|-------------|---------|
| one core (P1) | 0.008 | 0.008 | 1.00× |
| one sweep (P2) | 0.039 | 0.029 | 1.34× |

## Interpretation

- End-to-end fixed-λ fits: about **1.1–1.9×** faster in C++ on this machine.
- Larger relative gains when **R↔C++ overhead** matters (small \(n\), large \(d\)).
- When **rank is large**, Gram/solve dominate; R already uses C++ \(P_k^{\mathrm{full}}\) helpers → speedup ~1.1×.
- Useful for global-GCV (\(N_\lambda \times M_{\mathrm{GDF}}\) fits), not an order-of-magnitude win by itself.

See also `PROTOCOL_FIXED_LAMBDA_CPP.md`.
