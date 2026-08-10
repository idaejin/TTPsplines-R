# Rcpp / computational backends for TTPsplines
#
# Default conditional penalty path (2026-08-10):
#   tt_conditional_penalty_full(..., method = "auto")
#     → tt_conditional_penalty_full_env_cpp (left/own/right bond environments)
# Legacy unit-core reference (diagnostics / parity only):
#   method = "tt_cpp"
#
# Related exports:
#   tt_penalty_prepare_right_envs_cpp
#   tt_penalty_left_env_absorb_cpp
#   tt_penalty_from_envs_cpp
#   tt_conditional_penalty_full_env_cpp
#
# R fallbacks live in R/penalties.R if the shared library is unavailable.
# Spectral cGCV workspace remains in R (lambda.R); small m_k, already O(m) per eval.
