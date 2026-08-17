# Rcpp / computational backends for TTPsplines

## Production model (canonical)

ALS / PIRLS **sweeps** live in R (`R/als.R`, `R/pirls.R`) under the classical
global discrete penalty on Theta (`P_k^full`). Compiled code accelerates
**kernels** called from that loop — there is no full C++ ALS fitter yet.

Requesting `backend = "Rcpp"` for ALS/PIRLS warns and keeps the R sweep
(`resolve_backend()` in `R/control.R`).

## Kernels (use these)

```
tt_als_core_update_global_cpp   P1: one-core fixed-λ update under P_k^full
tt_gram_rhs_cpp                 fused / blocked / kron / blas Gram+RHS (P2b)
                               OpenMP opt-in via n_threads (fused_blocked)
tt_gram_omp_available
tt_conditional_penalty_full_cpp
tt_global_penalty_value_cpp
tt_penalty_prepare_right_envs_cpp
tt_penalty_left_env_absorb_cpp
tt_penalty_from_envs_cpp
tt_conditional_penalty_full_env_cpp
```

Default conditional penalty path (2026-08-10):

```
tt_conditional_penalty_full(..., method = "auto")
  → tt_conditional_penalty_full_env_cpp (left/own/right bond environments)
```

Legacy unit-core reference (diagnostics / parity only): `method = "tt_cpp"`.

R fallbacks live in `R/penalties.R` / `R/tt_gram.R` if the shared library is
unavailable. Spectral cGCV workspace remains in R (`R/lambda.R`).

## Legacy full fitters (do not use for production)

| Symbol | Issue |
|--------|--------|
| `tt_cgcv_fit_cpp` | Own-margin \(P_k\) only — not global \(P_k^{\mathrm{full}}\) |
| `tt_glm_pirls_cgcv_cpp` | Same; unused by `ttps()` after backend resolution |

Planned next step for global-GCV: **P2** one fixed-λ sweep in C++ (parity vs R),
then **P3** multi-sweep fitter. Outer λ selection stays in R. See
`inst/benchmarks/global_gcv/PROTOCOL_FIXED_LAMBDA_CPP.md`.
