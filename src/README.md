# Rcpp / computational backends for TTPsplines
#
# Default conditional penalty path (2026-08-10):
#   tt_conditional_penalty_full(..., method = "auto")
#     → tt_conditional_penalty_full_env_cpp (left/own/right bond environments)
# Legacy unit-core reference (diagnostics / parity only):
#   method = "tt_cpp"
#
# Related exports / kernels:
#   tt_gram_rhs_cpp — fused / blocked / kron / blas Gram+RHS (P2b)
#                 OpenMP opt-in via n_threads (fused_blocked only; P2b.3)
#   tt_gram_omp_available
#   tt_penalty_prepare_right_envs_cpp
#   tt_penalty_left_env_absorb_cpp
#   tt_penalty_from_envs_cpp
#   tt_conditional_penalty_full_env_cpp
#
# R fallbacks live in R/penalties.R if the shared library is unavailable.
# Spectral cGCV workspace remains in R (lambda.R); small m_k, already O(m) per eval.
