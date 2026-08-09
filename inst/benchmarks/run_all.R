# Run all package benchmarks (development).
#
#   cd ttpsplines-pkg
#   Rscript inst/benchmarks/run_all.R
#
# Optional:
#   TTPSPLINES_BENCH_OUT=/tmp/ttps_bench Rscript inst/benchmarks/run_all.R
#   TTPSPLINES_BENCH_WHICH=gaussian,glam Rscript inst/benchmarks/run_all.R

bench_dir <- (function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
  }
  if (file.exists("inst/benchmarks/run_all.R")) {
    return(normalizePath("inst/benchmarks"))
  }
  getwd()
})()

pkg_root <- normalizePath(file.path(bench_dir, "../.."))
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  # sourced from inst/benchmarks → ../.. is package root
  pkg_root <- normalizePath(file.path(bench_dir, "..", ".."))
}
# bench_dir is .../inst/benchmarks → parent of inst is package root
pkg_root <- normalizePath(file.path(bench_dir, "..", ".."))
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  pkg_root <- normalizePath(file.path(bench_dir, ".."))
}

message("Package root: ", pkg_root)
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  stop("Need devtools::load_all for development benchmarks")
}

Sys.setenv(TTPSPLINES_BENCH_OUT = Sys.getenv(
  "TTPSPLINES_BENCH_OUT",
  unset = file.path(bench_dir, "results")
))
dir.create(Sys.getenv("TTPSPLINES_BENCH_OUT"), showWarnings = FALSE, recursive = TRUE)

which <- Sys.getenv("TTPSPLINES_BENCH_WHICH", unset = "all")
all_scripts <- c(
  gaussian = "benchmark_gaussian.R",
  poisson = "benchmark_poisson.R",
  bernoulli = "benchmark_bernoulli.R",
  rank = "benchmark_rank.R",
  glam = "benchmark_glam.R"
)
if (identical(which, "all")) {
  todo <- unname(all_scripts)
} else {
  keys <- trimws(strsplit(which, ",", fixed = TRUE)[[1]])
  todo <- unname(all_scripts[keys])
  todo <- todo[!is.na(todo)]
}

t_all <- proc.time()[["elapsed"]]
for (script in todo) {
  message("\n######## ", script, " ########")
  source(file.path(bench_dir, script), local = new.env())
}
message(sprintf("\nAll done in %.1fs. Results in %s",
                proc.time()[["elapsed"]] - t_all,
                Sys.getenv("TTPSPLINES_BENCH_OUT")))
