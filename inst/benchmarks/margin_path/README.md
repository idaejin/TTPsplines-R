# Margin Activity Path

Formal screening of TT margins by partial-range activity + nested CV.

**Name:** Margin Activity Path (avoid acronym "MAP" = maximum a posteriori).

| Doc | Role |
|-----|------|
| [`PROTOCOL_MAP.md`](PROTOCOL_MAP.md) | **Canonical** method + selection rules |
| `tt_margin_activity_path()` | Exported package helper (`R/margin_activity_path.R`) |
| `vignette("margin-activity-path")` | User-facing walkthrough |
| [`run_map_example.R`](run_map_example.R) | Extended demo (+ optional linear lasso) |
| `results/` | CSV + `figures/` from the demo script |

```bash
Rscript /Users/daejin/Dropbox/IE/research/01_PROJECTS/ttpsplines-pkg/inst/benchmarks/margin_path/run_map_example.R
```

**Product:** KEEP `cGCV`. Screening layer only.
