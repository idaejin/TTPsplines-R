#' Control parameters for [ttpspline()].
#'
#' Separates algorithmic knobs from the model API. Optimizer-specific fields
#' for L-BFGS / Adam are reserved; Adam/Keras is optional and not required
#' for the default ALS + cGCV workflows.
#'
#' @param max_sweeps Maximum ALS sweeps (Gaussian / inner GLM).
#' @param pirls_maxit Maximum outer PIRLS iterations (GLM).
#' @param tol Relative deviance / RSS change for stopping.
#' @param tol_lambda Relative log-λ change for cGCV outer convergence.
#' @param lambda_start Initial λ used when `lambda = "cGCV"`.
#' @param lambda_bounds Length-2 bounds for automatic λ search on log scale.
#' @param lambda_tol Alias of `tol_lambda`.
#' @param lambda_update Reserved (`"auto"`); future releases may extend.
#' @param backend `"auto"`, `"R"`, `"Rcpp"`, or `"keras"` (Adam only).
#' @param sparse `"auto"`, `TRUE`, or `FALSE` (v0: dense bases; hybrid reserved).
#' @param use_spectral_gcv Use spectral cache inside Brent cGCV when feasible.
#' @param outer_maxit Outer alternation iters (LBFGS/Adam + cGCV).
#' @param outer_tol Outer relative change tolerance.
#' @param lbfgs_maxit Max L-BFGS iterations.
#' @param adam_lr,adam_epochs,adam_batch_size,adam_patience Adam/Keras knobs
#'   (optional backend; `adam_batch_size = NULL` means full-batch).
#' @param trace Print iteration progress.
#' @param damping Deviance step control for Bernoulli PIRLS.
#' @param seed RNG seed for TT core initialization.
#' @param init_sd SD for random TT initialization.
#' @param als_sweeps_per_pirls Inner ALS sweeps per PIRLS iteration.
#' @param compute_edf Compute joint linearized EDF after fit (`TRUE`/`FALSE`).
#'   Skipped automatically when packed TT size exceeds `edf_max_npar`.
#' @param edf_max_npar Maximum packed TT parameters for joint EDF (memory guard).
#' @return A list of class `"tt_control"`.
#' @export
tt_control <- function(max_sweeps = 20,
                       pirls_maxit = 50,
                       tol = 1e-6,
                       tol_lambda = 1e-4,
                       lambda_start = 1,
                       lambda_bounds = c(1e-4, 1e4),
                       lambda_tol = NULL,
                       lambda_update = "auto",
                       backend = c("auto", "R", "Rcpp", "keras"),
                       sparse = c("auto", TRUE, FALSE),
                       use_spectral_gcv = TRUE,
                       outer_maxit = 20,
                       outer_tol = 1e-5,
                       lbfgs_maxit = 500,
                       adam_lr = 1e-3,
                       adam_epochs = 1000,
                       adam_batch_size = NULL,
                       adam_patience = 30,
                       trace = FALSE,
                       damping = TRUE,
                       seed = 1,
                       init_sd = 0.15,
                       als_sweeps_per_pirls = 4,
                       compute_edf = TRUE,
                       edf_max_npar = 2500L) {
  backend <- match.arg(backend)
  if (is.character(sparse) && length(sparse) == 1L) {
    sparse <- match.arg(sparse, c("auto", "TRUE", "FALSE", "true", "false"))
    if (sparse %in% c("TRUE", "true")) sparse <- TRUE
    else if (sparse %in% c("FALSE", "false")) sparse <- FALSE
    # else "auto"
  } else if (!is.logical(sparse) && !identical(sparse, "auto")) {
    sparse <- "auto"
  }
  if (!is.null(lambda_tol)) tol_lambda <- lambda_tol
  structure(
    list(
      max_sweeps = as.integer(max_sweeps),
      pirls_maxit = as.integer(pirls_maxit),
      tol = as.numeric(tol),
      tol_lambda = as.numeric(tol_lambda),
      lambda_start = as.numeric(lambda_start),
      lambda_bounds = as.numeric(lambda_bounds),
      lambda_update = as.character(lambda_update),
      backend = backend,
      sparse = sparse,
      use_spectral_gcv = isTRUE(use_spectral_gcv),
      outer_maxit = as.integer(outer_maxit),
      outer_tol = as.numeric(outer_tol),
      lbfgs_maxit = as.integer(lbfgs_maxit),
      adam_lr = as.numeric(adam_lr),
      adam_epochs = as.integer(adam_epochs),
      adam_batch_size = if (is.null(adam_batch_size)) NULL else as.integer(adam_batch_size),
      adam_patience = as.integer(adam_patience),
      trace = isTRUE(trace),
      damping = isTRUE(damping),
      seed = as.integer(seed),
      init_sd = as.numeric(init_sd),
      als_sweeps_per_pirls = as.integer(als_sweeps_per_pirls),
      compute_edf = isTRUE(compute_edf),
      edf_max_npar = as.integer(edf_max_npar)
    ),
    class = "tt_control"
  )
}

#' Resolve computational backend given optimizer preference.
#' @keywords internal
resolve_backend <- function(control, optimizer = "ALS") {
  be <- control$backend
  if (identical(be, "auto")) {
    if (identical(optimizer, "Adam")) {
      if (isTRUE(tt_has_keras())) return("keras")
      return("keras") # caller will error with install hint
    }
    if (.ttpsplines_has_rcpp()) return("Rcpp")
    return("R")
  }
  if (identical(be, "Rcpp") && !.ttpsplines_has_rcpp()) {
    warning("Rcpp backend unavailable; falling back to R.", call. = FALSE)
    return("R")
  }
  if (identical(be, "keras") && !identical(optimizer, "Adam")) {
    warning("backend='keras' is only used with optimizer='Adam'; ignoring.",
            call. = FALSE)
    if (.ttpsplines_has_rcpp()) return("Rcpp")
    return("R")
  }
  be
}

.ttpsplines_has_rcpp <- function() {
  exists("tt_fit_d_cpp", mode = "function") ||
    exists("tt_cgcv_fit_cpp", mode = "function")
}
