# Optional Adam / Keras backend (stub until ALS+cGCV+Rcpp+GLM are solid).

#' Is the optional Keras/TensorFlow stack available?
#'
#' Does **not** install anything. Returns `TRUE` only if `reticulate`, Python,
#' TensorFlow, and Keras can be loaded.
#'
#' @export
tt_has_keras <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  ok <- FALSE
  tryCatch({
    # Avoid initializing Python aggressively if not configured
    if (!reticulate::py_available(initialize = FALSE)) return(FALSE)
    tf <- reticulate::import("tensorflow", delay_load = TRUE)
    # Touch a lightweight attribute
    invisible(tf$`__name__`)
    ok <- TRUE
  }, error = function(e) {
    ok <<- FALSE
  })
  ok
}

#' Diagnostic status for the optional Adam/Keras backend.
#'
#' @return A named list (print-friendly) describing availability.
#' @export
tt_keras_status <- function() {
  has_ret <- requireNamespace("reticulate", quietly = TRUE)
  py_ok <- FALSE
  tf_ok <- FALSE
  keras_ok <- FALSE
  detail <- character()
  if (!has_ret) {
    detail <- "Install Suggests package 'reticulate' to enable optional Adam."
  } else {
    py_ok <- isTRUE(tryCatch(
      reticulate::py_available(initialize = FALSE),
      error = function(e) FALSE
    ))
    if (!py_ok) {
      detail <- paste(
        "Python not available via reticulate.",
        "Configure a Python env with tensorflow/keras; TTPsplines will not auto-install."
      )
    } else {
      tf_ok <- isTRUE(tryCatch({
        reticulate::import("tensorflow")
        TRUE
      }, error = function(e) FALSE))
      keras_ok <- isTRUE(tryCatch({
        reticulate::import("keras")
        TRUE
      }, error = function(e) {
        # TF 2.x often exposes keras as tensorflow.keras
        tryCatch({
          tf <- reticulate::import("tensorflow")
          invisible(tf$keras)
          TRUE
        }, error = function(e2) FALSE)
      }))
      if (!tf_ok || !keras_ok) {
        detail <- paste(
          "Python is available but TensorFlow/Keras were not found.",
          "Install tensorflow in the active Python env, then retry optimizer='Adam'."
        )
      } else {
        detail <- "Keras/TensorFlow appear available for optional Adam backend."
      }
    }
  }
  out <- list(
    reticulate = has_ret,
    python = py_ok,
    tensorflow = tf_ok,
    keras = keras_ok,
    ready = isTRUE(tf_ok && keras_ok),
    detail = detail
  )
  class(out) <- c("tt_keras_status", "list")
  out
}

#' @export
print.tt_keras_status <- function(x, ...) {
  cat("TTPsplines optional Keras/Adam status\n")
  cat(sprintf("  reticulate:  %s\n", x$reticulate))
  cat(sprintf("  python:      %s\n", x$python))
  cat(sprintf("  tensorflow:  %s\n", x$tensorflow))
  cat(sprintf("  keras:       %s\n", x$keras))
  cat(sprintf("  ready:       %s\n", x$ready))
  cat(sprintf("  note:        %s\n", x$detail))
  invisible(x)
}

#' Adam/Keras fit — not implemented yet (clear blocker message).
#' @keywords internal
tt_adam_fit <- function(y, basis, ranks, lambda_spec, control,
                        penalty_order = 2, init_cores = NULL,
                        family = NULL, offset = NULL, weights = NULL) {
  st <- tt_keras_status()
  msg <- paste0(
    "optimizer = \"Adam\" requires the optional Keras/TensorFlow backend ",
    "and is not implemented yet.\n",
    "Use optimizer = \"ALS\" (default) or \"LBFGS\".\n",
    "Keras status: ", st$detail, "\n",
    "Call tt_keras_status() / tt_has_keras() for diagnostics. ",
    "TTPsplines will not auto-install Python packages."
  )
  stop(msg, call. = FALSE)
}
