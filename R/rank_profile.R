#' Fit TT P-splines over a grid of ranks (in-sample rank diagnostic).
#'
#' Sweeps ranks on the **full** training sample and reports deviance /
#' complexity. This is **not** cross-validated predictive selection; for
#' that use [tt_rank_select()].
#'
#' @inheritParams ttpspline
#' @param ranks Numeric vector of scalar ranks to try.
#' @param ... Passed to [ttpspline()] (overrides `rank`).
#' @return Object of class `"tt_rank_profile"`.
#' @seealso [tt_rank_select()]
#' @export
tt_rank_profile <- function(y,
                            X,
                            ranks = c(1, 2, 3, 4, 6),
                            family = stats::gaussian(),
                            lambda = "cGCV",
                            k = 10,
                            control = tt_control(),
                            ...) {
  rows <- lapply(ranks, function(r) {
    t0 <- proc.time()[["elapsed"]]
    fit <- ttpspline(
      y, X, family = family, rank = r, k = k, lambda = lambda,
      control = control, ...
    )
    data.frame(
      rank = r,
      npar_tt = fit$npar_tt,
      npar_dense = fit$npar_dense,
      compression = fit$compression_ratio,
      deviance = fit$deviance,
      time_s = proc.time()[["elapsed"]] - t0,
      converged = fit$converged,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  structure(
    list(table = tab, family = normalize_family(family), lambda = lambda, k = k),
    class = "tt_rank_profile"
  )
}

#' @export
print.tt_rank_profile <- function(x, ...) {
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' @export
plot.tt_rank_profile <- function(x, ...) {
  tab <- x$table
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  on.exit(par(op), add = TRUE)
  plot(tab$rank, tab$deviance, type = "b", pch = 19,
       xlab = "rank r", ylab = "deviance / RSS", main = "Fit quality", ...)
  plot(tab$rank, tab$npar_tt, type = "b", pch = 19,
       xlab = "rank r", ylab = "TT parameters", main = "Complexity", ...)
  invisible(x)
}
