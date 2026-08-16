# Matched-budget gate report

**Mode:** SMOKE

**Date:** 2026-08-16 23:27

**Verdict:** GO

sobol_refine valley_rate=1.00 vs random=0.00 lhs=0.50 at N=32

Product lock: KEEP cGCV. Protocol: `PROTOCOL_MATCHED_BUDGET.md`.

## Aggregate (mean over design seeds)

        method  N      scenario in_valley_1pct   regret_rel n_als_total elapsed
1          lhs 32 smooth_smooth            0.5 0.0190263783         203  1.3995
2       random 32 smooth_smooth            0.0 0.0359909527         203  1.2325
3   sobol_pure 32 smooth_smooth            0.0 0.0454503692         203  1.3855
4 sobol_refine 32 smooth_smooth            1.0 0.0003402931        1217  5.1090
  n_valley n_seeds valley_rate
1        1       2         0.5
2        0       2         0.0
3        0       2         0.0
4        2       2         1.0

## Full summary

       scenario       method d  N seed_design seed_data    theta1    theta2
1 smooth_smooth         cgcv 2  0          NA         1 -1.958892 -2.925382
2 smooth_smooth       random 2 32           1         1 -1.596510 -2.063966
3 smooth_smooth          lhs 2 32           1         1 -1.735324 -2.546316
4 smooth_smooth   sobol_pure 2 32           1         1 -1.719913 -1.903761
5 smooth_smooth sobol_refine 2 32           1         1 -1.865340 -2.715552
6 smooth_smooth       random 2 32          11         1 -2.100250 -4.421411
7 smooth_smooth          lhs 2 32          11         1 -2.381167 -4.444595
8 smooth_smooth   sobol_pure 2 32          11         1 -1.758752 -2.338567
9 smooth_smooth sobol_refine 2 32          11         1 -1.865400 -2.715643
   q_search   q_final     q_ref in_valley_1pct   regret_rel    d_theta aniso_ok
1        NA 0.1766811 0.1766811           TRUE 0.0000000000 0.08519252       NA
2 0.1797577 0.1868330 0.1766811          FALSE 0.0574591341 1.01929549     TRUE
3 0.1724695 0.1779773 0.1766811           TRUE 0.0073366227 0.52524521     TRUE
4 0.1823404 0.1897310 0.1766811          FALSE 0.0738612219 1.13145413     TRUE
5 0.1717999 0.1767413 0.1766811           TRUE 0.0003410646 0.31471300     TRUE
6 0.1756358 0.1792470 0.1766811          FALSE 0.0145227713 1.42494136     TRUE
7 0.1789782 0.1821080 0.1766811          FALSE 0.0307161339 1.49403551     TRUE
8 0.1736982 0.1796916 0.1766811          FALSE 0.0170395164 0.70405565     TRUE
9 0.1717999 0.1767411 0.1766811           TRUE 0.0003395215 0.31460461     TRUE
  winner_source n_explore n_refine_starts n_refine_evals n_failed n_als_search
1          cgcv         0               0              0        0           NA
2  explore_best        32               0              0        0          192
3  explore_best        32               0              0        1          192
4  explore_best        32               0              0        1          192
5       refined        32               3            178        1         1260
6  explore_best        32               0              0        0          192
7  explore_best        32               0              0        0          192
8  explore_best        32               0              0        0          192
9       refined        32               3            160        0         1152
  n_als_total n_theta_miss elapsed
1          NA           NA   0.231
2         203           32   1.428
3         203           32   1.205
4         203           32   1.834
5        1271          210   6.062
6         203           32   1.037
7         203           32   1.594
8         203           32   0.937
9        1163          192   4.156
