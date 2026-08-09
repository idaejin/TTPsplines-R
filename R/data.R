#' Ishigami test sample (d = 3)
#'
#' Scattered design on \eqn{[-\pi,\pi]^3} with Gaussian noise around the
#' Ishigami surface (see [f_ishigami()] / [simulate_ishigami()]).
#'
#' @format A data frame with 800 rows and 5 columns:
#' \describe{
#'   \item{x1,x2,x3}{Inputs on \eqn{[-\pi,\pi]}.}
#'   \item{y}{Noisy response.}
#'   \item{f}{True Ishigami value (noiseless).}
#' }
#' @seealso [simulate_ishigami()], [f_ishigami()]
#' @examples
#' data(ishigami)
#' X <- as.matrix(ishigami[, c("x1", "x2", "x3")])
#' fit <- ttpspline(ishigami$y, X, rank = 2, k = 8, lambda = 1,
#'                  control = tt_control(max_sweeps = 8, compute_edf = FALSE))
"ishigami"

#' Sobol g-function sample (d = 4)
#'
#' Scattered design on \eqn{[0,1]^4} with
#' \eqn{a=(0,0.5,3,9)} (see [f_sobol_g()] / [simulate_sobol_g()]).
#'
#' @format A data frame with 800 rows and 6 columns:
#' \describe{
#'   \item{x1,x2,x3,x4}{Inputs on \eqn{[0,1]}.}
#'   \item{y}{Noisy response.}
#'   \item{f}{True g-function value.}
#' }
#' @seealso [simulate_sobol_g()], [f_sobol_g()]
#' @examples
#' data(sobol_g)
#' X <- as.matrix(sobol_g[, paste0("x", 1:4)])
#' fit <- ttpspline(sobol_g$y, X, rank = 2, k = 6, lambda = 1,
#'                  control = tt_control(max_sweeps = 6, compute_edf = FALSE))
"sobol_g"

#' Friedman #1 test sample (d = 5)
#'
#' Scattered design on \eqn{[0,1]^5} (see [f_friedman()] / [simulate_friedman()]).
#'
#' @format A data frame with 800 rows and 7 columns:
#' \describe{
#'   \item{x1,...,x5}{Inputs on \eqn{[0,1]}.}
#'   \item{y}{Noisy response.}
#'   \item{f}{True Friedman surface.}
#' }
#' @seealso [simulate_friedman()], [f_friedman()]
#' @examples
#' data(friedman)
#' X <- as.matrix(friedman[, paste0("x", 1:5)])
#' fit <- ttpspline(friedman$y, X, rank = 3, k = 6, lambda = 1,
#'                  control = tt_control(max_sweeps = 8, compute_edf = FALSE))
"friedman"

#' Currie–Durbán–Eilers Poisson age × year array
#'
#' Simulated mortality-style counts on a regular age–period grid with known
#' smooth log-rate and exposures. Intended for [glam_fit_poisson()] demos
#' (Currie–Durbán–Eilers GLAM array methods).
#'
#' @format A list with components:
#' \describe{
#'   \item{Y}{41 × 31 integer count array.}
#'   \item{exposure}{Matching exposure array.}
#'   \item{age, year}{Grid axes.}
#'   \item{eta, mu}{True log-rate surface and mean counts.}
#' }
#' @seealso [simulate_glam_poisson()], [glam_fit_poisson()], [glam_grid_bases()]
#' @references
#' Currie, I. D., Durban, M. and Eilers, P. H. C. (2006).
#' Generalized linear array models with applications to multidimensional
#' smoothing. *JRSS-B*.
#' @examples
#' data(glam_poisson)
#' bb <- glam_grid_bases(list(age = glam_poisson$age, year = glam_poisson$year), k = 8)
#' fit <- glam_fit_poisson(glam_poisson$Y, bb$B, lambda = c(10, 1),
#'                         offset = log(glam_poisson$exposure))
#' fit$n_pirls
"glam_poisson"
