#' @keywords internal
"_PACKAGE"

#' @useDynLib TTPsplines, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats gaussian poisson binomial dbinom dpois dnorm
#' @importFrom stats optimize plogis qlogis rnorm runif rpois rbinom
#' @importFrom stats median quantile sd var cor
#' @importFrom graphics plot par abline image legend title matplot
#' @importFrom utils head modifyList
NULL
