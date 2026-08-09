#' Control parameters for [ttpspline()].
#'
#' @param max_sweeps Maximum ALS sweeps (Gaussian / inner GLM).
#' @param pirls_maxit Maximum outer PIRLS iterations (GLM).
#' @param tol Relative deviance / RSS change for stopping.
#' @param tol_lambda Relative log-λ change for cGCV outer convergence.
#' @param lambda_bounds Length-2 bounds for automatic λ search.
#' @param backend `"auto"`, `"R"`, or `"Rcpp"`.
#' @param sparse `"auto"`, `TRUE`, or `FALSE` (v0 uses dense bases; flag reserved).
#' @param trace Print iteration progress.
#' @param damping Deviance step control for Bernoulli PIRLS (reserved / soft).
#' @param seed RNG seed for TT core initialization.
#' @param init_sd SD for random TT initialization.
#' @param als_sweeps_per_pirls Inner ALS sweeps per PIRLS iteration.
#' @return A list of class `"tt_control"`.
#' @export
tt_control <- function(max_sweeps = 20,
                       pirls_maxit = 50,
                       tol = 1e-6,
                       tol_lambda = 1e-4,
                       lambda_bounds = c(1e-4, 1e4),
                       backend = c("auto", "R", "Rcpp"),
                       sparse = c("auto", TRUE, FALSE),
                       trace = FALSE,
                       damping = TRUE,
                       seed = 1,
                       init_sd = 0.15,
                       als_sweeps_per_pirls = 4) {
  backend <- match.arg(backend)
  if (is.character(sparse)) sparse <- match.arg(sparse)
  structure(
    list(
      max_sweeps = as.integer(max_sweeps),
      pirls_maxit = as.integer(pirls_maxit),
      tol = as.numeric(tol),
      tol_lambda = as.numeric(tol_lambda),
      lambda_bounds = as.numeric(lambda_bounds),
      backend = backend,
      sparse = sparse,
      trace = isTRUE(trace),
      damping = isTRUE(damping),
      seed = as.integer(seed),
      init_sd = as.numeric(init_sd),
      als_sweeps_per_pirls = as.integer(als_sweeps_per_pirls)
    ),
    class = "tt_control"
  )
}

resolve_backend <- function(control) {
  be <- control$backend
  if (identical(be, "auto")) {
    if (.ttpsplines_has_rcpp()) return("Rcpp")
    return("R")
  }
  if (identical(be, "Rcpp") && !.ttpsplines_has_rcpp()) {
    warning("Rcpp backend unavailable; falling back to R.", call. = FALSE)
    return("R")
  }
  be
}

.ttpsplines_has_rcpp <- function() {
  exists("tt_fit_d_cpp", mode = "function") ||
    exists("tt_cgcv_fit_cpp", mode = "function")
}
