# Direct-likelihood optimizers (GD / LBFGS / Adam)

TTPsplines separates **structure-aware** estimation from **direct** minimization of the same penalized objective

\[
\mathcal L(G,\alpha;\lambda)
=
-\ell(y;\eta_G)
+\tfrac12\sum_k\lambda_k\,J_k(G_k).
\]

| `optimizer` | Philosophy | Role |
|-------------|------------|------|
| `ALS` / `PIRLS-ALS` | Conditional core solves (IRLS for GLM) | Production structure-aware path |
| `Damped-Newton-ALS` | Conditional Newton + Armijo on true \(Q_k\) | Experimental (see `OPTIMIZER_CONDITIONAL_GLM.md`) |
| `LBFGS-ALS` | Conditional L-BFGS on each \(Q_k\) | Experimental / benchmark |
| `GD` | First-order GD + Armijo on \(\mathcal L\) | Experimental / paper baseline |
| `LBFGS` | Quasi-Newton on \(\mathcal L\) | Default for Bernoulli (`auto`) |
| `Adam` | Adaptive first-order (Keras; stub) | Future |

`GD` reuses the **same analytical objective and gradient** as `LBFGS` (`R/optimizer_lbfgs.R`). It is **not** a PIRLS variant.

Controls: `gd_lr`, `gd_maxit`, `gd_tol`, `gd_linesearch`, `gd_step_factor`, `gd_step_min`, `gd_armijo_c`.

## Bernoulli gate result (2026-08-09, `--full`)

Script: `inst/benchmarks/benchmark_bernoulli_gd.R`  
Results: `inst/benchmarks/results/bernoulli_gd/`

**Decision (ranks 2 and 3, λ=1, 20 seeds):**

\[
\mathrm{GD}\approx\mathrm{LBFGS}\gg\mathrm{PIRLS\text{-}ALS}
\]

on test RMSE(\(\eta\)). Interpretation: the Bernoulli gap is driven by the **PIRLS/alternating trajectory**, not by L-BFGS being uniquely privileged. GD matches LBFGS quality but is ~40–50× slower.

| rank | method | med RMSE | med max\|η\| | med time |
|-----:|--------|---------:|-------------:|---------:|
| 2 | PIRLS-ALS | 0.933 | 5.52 | 0.7s |
| 2 | GD | 0.507 | 4.10 | 25s |
| 2 | LBFGS | 0.498 | 4.13 | 0.6s |
| 3 | PIRLS-ALS | 1.139 | 7.70 | 0.2s |
| 3 | GD | 0.693 | 5.02 | 40s |
| 3 | LBFGS | 0.731 | 5.37 | 0.9s |

At rank 4 / small λ the landscape is harder for all methods (including occasional LBFGS blowups); see `exp3_rank_lambda.csv`.

**Follow-up:** conditional solvers gate (`OPTIMIZER_CONDITIONAL_GLM.md`) — Scenario **E**: core-wise ALS (incl. LBFGS-ALS / Damped-Newton) still ≪ global LBFGS on Bernoulli.
