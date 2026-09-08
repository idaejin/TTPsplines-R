#!/usr/bin/env Rscript
# Render vignettes against the *installed* package (no load_all).
options(warn = 1)
root <- normalizePath(".")
stopifnot(file.exists(file.path(root, "DESCRIPTION")))

if ("package:TTPsplines" %in% search()) {
  detach("package:TTPsplines", unload = TRUE, character.only = TRUE)
}
try(unloadNamespace("TTPsplines"), silent = TRUE)

message("Installing package from ", root, " (no vignettes) ...")
utils::install.packages(
  root,
  repos = NULL,
  type = "source",
  quiet = TRUE,
  INSTALL_opts = "--no-byte-compile"
)

library(TTPsplines)
stopifnot("tt_ic" %in% getNamespaceExports("TTPsplines"))
message("tt_ic export OK")

files <- sort(list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE))
outdir <- file.path(tempdir(), paste0("vig-inst-", as.integer(Sys.time())))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("Rendering ", length(files), " vignettes -> ", outdir)

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
utils::write.csv(res, "vignettes/_build_check.csv", row.names = FALSE)
if (!all(res$ok)) {
  quit(status = 1L)
}
message("All vignettes OK against installed package")
