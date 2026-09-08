#!/usr/bin/env Rscript
# Build every vignette and report failures.
options(warn = 1)
root <- if (file.exists("DESCRIPTION")) {
  getwd()
} else if (file.exists("../DESCRIPTION")) {
  normalizePath("..")
} else {
  stop("Run from package root or scripts/")
}
setwd(root)

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Need devtools")
}
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Need rmarkdown")
}

devtools::load_all(".", quiet = TRUE)
outdir <- file.path(tempdir(), paste0("vig-", as.integer(Sys.time())))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE))
message("Found ", length(files), " vignettes; output -> ", outdir)

res <- data.frame(
  vignette = basename(files),
  ok = FALSE,
  error = NA_character_,
  time = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(files)) {
  f <- files[[i]]
  message("\n==== ", basename(f), " ====")
  t0 <- proc.time()[["elapsed"]]
  ok <- FALSE
  err <- NA_character_
  tryCatch(
    {
      rmarkdown::render(
        input = f,
        output_format = rmarkdown::html_vignette(),
        output_dir = outdir,
        quiet = TRUE,
        envir = new.env(parent = globalenv())
      )
      ok <- TRUE
    },
    error = function(e) {
      err <<- conditionMessage(e)
      message("ERROR: ", err)
    }
  )
  res$ok[i] <- ok
  res$error[i] <- err
  res$time[i] <- proc.time()[["elapsed"]] - t0
  message(if (ok) "OK" else "FAIL", sprintf(" (%.1fs)", res$time[i]))
}

print(res, right = FALSE)
out_csv <- file.path("vignettes", "_build_check.csv")
utils::write.csv(res, out_csv, row.names = FALSE)
message("Wrote ", out_csv)
if (!all(res$ok)) {
  quit(status = 1L)
}
