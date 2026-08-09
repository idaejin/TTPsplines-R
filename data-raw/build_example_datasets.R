# Build packaged example datasets (Ishigami, Sobol-g, Friedman).
#   Rscript data-raw/build_example_datasets.R

fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(fa)) {
  normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), ".."))
} else if (file.exists("DESCRIPTION")) {
  normalizePath(".")
} else if (file.exists("../DESCRIPTION")) {
  normalizePath("..")
} else {
  stop("Run from package root or via Rscript data-raw/build_example_datasets.R")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required")
}
devtools::load_all(root, quiet = TRUE)

ishigami <- as.data.frame(simulate_ishigami(n = 800, sigma = 0.15, seed = 1L))
attr(ishigami, "name") <- "ishigami"
attr(ishigami, "a") <- 7
attr(ishigami, "b") <- 0.1
attr(ishigami, "sigma") <- 0.15
attr(ishigami, "seed") <- 1L

sobol_g <- as.data.frame(simulate_sobol_g(
  n = 800, d = 4L, a = c(0, 0.5, 3, 9), sigma = 0.05, seed = 2L
))
attr(sobol_g, "name") <- "sobol_g"
attr(sobol_g, "a") <- c(0, 0.5, 3, 9)
attr(sobol_g, "sigma") <- 0.05
attr(sobol_g, "seed") <- 2L

friedman <- as.data.frame(simulate_friedman(n = 800, sigma = 1, seed = 3L))
attr(friedman, "name") <- "friedman"
attr(friedman, "sigma") <- 1
attr(friedman, "seed") <- 3L

dir.create(file.path(root, "data"), showWarnings = FALSE)
save(ishigami, file = file.path(root, "data", "ishigami.rda"), compress = "xz")
save(sobol_g, file = file.path(root, "data", "sobol_g.rda"), compress = "xz")
save(friedman, file = file.path(root, "data", "friedman.rda"), compress = "xz")
message("Wrote data/{ishigami,sobol_g,friedman}.rda under ", root)
