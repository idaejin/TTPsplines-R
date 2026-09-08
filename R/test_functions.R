# Classic UQ / nonparametric regression test surfaces for TT-P-spline demos.
# These are *response surfaces*, not Sobol-index estimators.

#' Coerce x to an n x d matrix (length-d vectors become 1 x d).
#' @keywords internal
#' @noRd
.as_input_matrix <- function(x, d = NULL, name = "function", min_d = NULL) {
  if (is.null(dim(x)) && is.numeric(x)) {
    if (!is.null(d) && length(x) == d) {
      X <- matrix(x, nrow = 1L)
    } else if (!is.null(min_d) && length(x) >= min_d) {
      X <- matrix(x, nrow = 1L)
    } else {
      X <- matrix(x, ncol = 1L)
    }
  } else {
    X <- as.matrix(x)
  }
  storage.mode(X) <- "double"
  if (!is.null(d) && ncol(X) != d) {
    stop(name, " requires exactly ", d, " inputs.", call. = FALSE)
  }
  if (!is.null(min_d) && ncol(X) < min_d) {
    stop(name, " requires at least ", min_d, " inputs.", call. = FALSE)
  }
  X
}

#' Scale columns of X to a target interval (affine).
#' @keywords internal
#' @noRd
.scale_to_interval <- function(X, from = c(0, 1), to = c(0, 1)) {
  X <- as.matrix(X)
  a0 <- from[1]; b0 <- from[2]
  a1 <- to[1]; b1 <- to[2]
  u <- (X - a0) / (b0 - a0)
  a1 + u * (b1 - a1)
}

#' Ishigami test function (d = 3).
#'
#' Standard form on \eqn{x_j \in [-\pi,\pi]}:
#' \deqn{f(x)=\sin x_1 + a\sin^2 x_2 + b\, x_3^4\sin x_1}
#' with defaults \eqn{a=7}, \eqn{b=0.1}.
#'
#' @param x Numeric matrix / data frame with 3 columns (or vector of length 3).
#' @param a,b Ishigami parameters.
#' @return Numeric vector of function values.
#' @references Ishigami, T. and Homma, T. (1990).
#' @export
#' @examples
#' X <- matrix(runif(30, -pi, pi), 10, 3)
#' f_ishigami(X)
f_ishigami <- function(x, a = 7, b = 0.1) {
  X <- .as_input_matrix(x, d = 3L, name = "Ishigami")
  x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
  sin(x1) + a * sin(x2)^2 + b * x3^4 * sin(x1)
}

#' Sobol g-function (Saltelli).
#'
#' On \eqn{x_j\in[0,1]^d}:
#' \deqn{g(x)=\prod_{j=1}^d\frac{|4x_j-2|+a_j}{1+a_j}.}
#' Smaller \eqn{a_j} ⇒ stronger main effect for coordinate \eqn{j}.
#'
#' @param x Numeric matrix / data frame (`n × d`).
#' @param a Nonnegative coefficients (length `d`, or recycled). Classic
#'   Saltelli example uses `a = c(0, 0.5, 3, 9, 99, 99, ...)`.
#' @return Numeric vector of function values.
#' @references Saltelli, A. et al. (2000). *Sensitivity Analysis in Practice*.
#' @export
#' @examples
#' X <- matrix(runif(40), 10, 4)
#' f_sobol_g(X, a = c(0, 0.5, 3, 9))
f_sobol_g <- function(x, a = c(0, 0.5, 3, 9)) {
  a <- as.numeric(a)
  if (is.null(dim(x)) && is.numeric(x)) {
    X <- .as_input_matrix(x, d = length(a), name = "Sobol g-function")
  } else {
    X <- .as_input_matrix(x, name = "Sobol g-function")
  }
  d <- ncol(X)
  a <- rep(a, length.out = d)
  if (any(a < 0)) stop("Sobol g-function requires a_j >= 0.", call. = FALSE)
  g <- rep(1, nrow(X))
  for (j in seq_len(d)) {
    g <- g * (abs(4 * X[, j] - 2) + a[j]) / (1 + a[j])
  }
  g
}

#' Friedman #1 test function (d = 5).
#'
#' On \eqn{x_j\in[0,1]^5}:
#' \deqn{f(x)=10\sin(\pi x_1 x_2)+20(x_3-1/2)^2+10 x_4+5 x_5.}
#'
#' @param x Numeric matrix / data frame with at least 5 columns (extra ignored).
#' @return Numeric vector of function values.
#' @references Friedman, J. H. (1991). Multivariate adaptive regression splines.
#' @export
#' @examples
#' X <- matrix(runif(50), 10, 5)
#' f_friedman(X)
f_friedman <- function(x) {
  X <- .as_input_matrix(x, min_d = 5L, name = "Friedman #1")
  10 * sin(pi * X[, 1] * X[, 2]) +
    20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] +
    5 * X[, 5]
}

#' Build a scattered regression sample from a test surface.
#'
#' @param fun Function of a numeric matrix `X`.
#' @param n Sample size.
#' @param d Dimension (columns of `X`).
#' @param domain Length-2 range for each margin, or `d × 2` matrix of bounds.
#' @param sigma Gaussian noise SD (`0` ⇒ noiseless).
#' @param seed Optional RNG seed.
#' @param family One of `"gaussian"`, `"poisson"`, `"binomial"`. For Poisson /
#'   binomial, `fun(X)` is treated as a latent surface then mapped to mean /
#'   probability (see details in [simulate_ishigami()]).
#' @param ... Passed to `fun`.
#' @return A list with `X`, `y`, `f` (latent/true mean on the link / Gaussian
#'   scale as appropriate), `family`, and metadata.
#' @keywords internal
#' @noRd
.simulate_test_surface <- function(fun, n, d, domain = c(0, 1), sigma = 0.1,
                                   seed = NULL, family = "gaussian", ...) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  n <- as.integer(n); d <- as.integer(d)
  if (is.matrix(domain) || is.data.frame(domain)) {
    dom <- as.matrix(domain)
    if (nrow(dom) != d || ncol(dom) != 2L) {
      stop("domain matrix must be d x 2.", call. = FALSE)
    }
  } else {
    domain <- as.numeric(domain)
    if (length(domain) != 2L) stop("domain must be length 2 or d x 2.", call. = FALSE)
    dom <- matrix(domain, nrow = d, ncol = 2L, byrow = TRUE)
  }
  U <- matrix(stats::runif(n * d), n, d)
  X <- matrix(NA_real_, n, d)
  for (j in seq_len(d)) {
    X[, j] <- dom[j, 1] + U[, j] * (dom[j, 2] - dom[j, 1])
  }
  colnames(X) <- paste0("x", seq_len(d))
  f <- as.numeric(fun(X, ...))
  family <- match.arg(family, c("gaussian", "poisson", "binomial"))
  if (identical(family, "gaussian")) {
    y <- f + stats::rnorm(n, 0, sigma)
  } else if (identical(family, "poisson")) {
    # centre and shift so means are moderate
    eta <- f - mean(f) + log(3)
    mu <- exp(eta)
    y <- stats::rpois(n, mu)
    f <- eta
  } else {
    eta <- 1.2 * (f - mean(f)) / stats::sd(f)
    p <- stats::plogis(eta)
    y <- stats::rbinom(n, 1L, p)
    f <- eta
  }
  list(
    X = X, y = y, f = f, family = family,
    n = n, d = d, domain = dom, sigma = sigma, seed = seed
  )
}

#' Simulate an Ishigami sample for TT-P-spline demos.
#'
#' Inputs are drawn uniformly on \eqn{[-\pi,\pi]^3}.
#'
#' @param n Sample size.
#' @param sigma Gaussian noise SD (ignored unless `family = "gaussian"`).
#' @param seed RNG seed.
#' @param a,b Ishigami parameters.
#' @param family `"gaussian"`, `"poisson"`, or `"binomial"`.
#' @return List with `X`, `y`, `f`, and metadata; also class `"ttps_sim"`.
#' @export
#' @examples
#' dat <- simulate_ishigami(n = 400, sigma = 0.15, seed = 1)
#' fit <- ttps(dat$y, dat$X, rank = 2, k = 8, lambda = 1,
#'                  control = tt_control(max_sweeps = 8, compute_edf = FALSE))
simulate_ishigami <- function(n = 800, sigma = 0.15, seed = 1,
                              a = 7, b = 0.1, family = "gaussian") {
  out <- .simulate_test_surface(
    f_ishigami, n = n, d = 3L, domain = c(-pi, pi),
    sigma = sigma, seed = seed, family = family, a = a, b = b
  )
  out$name <- "ishigami"
  out$a <- a
  out$b <- b
  class(out) <- c("ttps_sim", "list")
  out
}

#' Simulate a Sobol g-function sample.
#'
#' @param n Sample size.
#' @param d Dimension (default 4).
#' @param a Coefficients (recycled to length `d`).
#' @param sigma Gaussian noise SD.
#' @param seed RNG seed.
#' @param family `"gaussian"`, `"poisson"`, or `"binomial"`.
#' @return A `"ttps_sim"` list (`X`, `y`, `f`, …).
#' @export
#' @examples
#' dat <- simulate_sobol_g(n = 500, d = 4, seed = 2)
#' fit <- ttps(dat$y, dat$X, rank = 2, k = 6, lambda = 1,
#'                  control = tt_control(max_sweeps = 6, compute_edf = FALSE))
simulate_sobol_g <- function(n = 800, d = 4,
                             a = c(0, 0.5, 3, 9, 99, 99),
                             sigma = 0.05, seed = 2,
                             family = "gaussian") {
  d <- as.integer(d)
  a <- rep(as.numeric(a), length.out = d)
  out <- .simulate_test_surface(
    f_sobol_g, n = n, d = d, domain = c(0, 1),
    sigma = sigma, seed = seed, family = family, a = a
  )
  out$name <- "sobol_g"
  out$a <- a
  class(out) <- c("ttps_sim", "list")
  out
}

#' Simulate a Friedman #1 sample (d = 5).
#'
#' @param n Sample size.
#' @param sigma Gaussian noise SD.
#' @param seed RNG seed.
#' @param family `"gaussian"`, `"poisson"`, or `"binomial"`.
#' @return A `"ttps_sim"` list (`X`, `y`, `f`, …).
#' @export
#' @examples
#' dat <- simulate_friedman(n = 600, sigma = 1, seed = 3)
#' fit <- ttps(dat$y, dat$X, rank = 3, k = 6, lambda = 1,
#'                  control = tt_control(max_sweeps = 8, compute_edf = FALSE))
simulate_friedman <- function(n = 800, sigma = 1, seed = 3,
                              family = "gaussian") {
  out <- .simulate_test_surface(
    f_friedman, n = n, d = 5L, domain = c(0, 1),
    sigma = sigma, seed = seed, family = family
  )
  out$name <- "friedman"
  class(out) <- c("ttps_sim", "list")
  out
}

#' Dette--Pepelyshev 8D test function.
#'
#' The function is defined on \eqn{[0,1]^8} and has non-negligible effects in
#' all eight coordinates as well as strong low-order interactions.
#' Reference: Dette and Pepelyshev (2010),
#' \url{https://www.sfu.ca/~ssurjano/detpep108d.html}.
#'
#' @param x Numeric matrix / data frame with at least 8 columns (extra ignored).
#'   Inputs must lie in \eqn{[0,1]^8}.
#' @return Numeric vector of function values.
#' @references Dette, H. and Pepelyshev, A. (2010). Generalized Latin hypercube
#'   design for computer experiments. *Technometrics* **52**, 421--429.
#' @export
#' @examples
#' X <- matrix(runif(80), 10, 8)
#' f_dette(X)

# Internal: vectorised Dette-Pepelyshev (all 8 terms, standard form).
#' @keywords internal
#' @noRd
.f_dette_internal <- function(x) {
  X <- as.matrix(x)
  x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
  x4 <- X[, 4]; x5 <- X[, 5]; x6 <- X[, 6]
  x7 <- X[, 7]; x8 <- X[, 8]
  4 * (x1 - 2 + 8 * x2 - 8 * x2^2)^2 +
    (3 - 4 * x2)^2 +
    16 * sqrt(x3 + 1) * (2 * x3 - 1)^2 +
    (x4 - x3)^2 +
    (x5 - 1)^2 +
    (x6 - 1)^2 +
    (x7 - 1)^2 +
    (x8 - 1)^2
}

#' @export
f_dette <- function(x) {
  X <- .as_input_matrix(x, min_d = 8L, name = "Dette-Pepelyshev 8D")
  .f_dette_internal(X)
}

#' Simulate a Dette--Pepelyshev 8D point-process / Poisson sample.
#'
#' Inputs are drawn uniformly on \eqn{[0,1]^8}.  For `d < 8`, only the first
#' `d` inputs vary; the rest are fixed at 0.5 (matching the PDF experiment).
#'
#' @param n_events Approximate number of Poisson events (controls `rate_scale`).
#' @param d Number of active dimensions (1--8; remaining fixed at 0.5).
#' @param n_bins Number of grid bins per axis.
#' @param rate_scale Scalar that controls the overall event rate.
#' @param rate_offset Log-intensity offset (shifts overall level).
#' @param holdout_frac Fraction of occupied cells held out for testing.
#' @param seed RNG seed.
#' @return A list with components:
#'   \describe{
#'     \item{`X_train`, `X_test`}{Cell centres (\eqn{n \times d}) for
#'       training and test occupied cells.}
#'     \item{`y_train`, `y_test`}{Event counts.}
#'     \item{`eta_true`}{True log-intensity function (closure).}
#'     \item{`d`, `n_bins`, `n_events`, `seed`}{Metadata.}
#'   }
#' @export
simulate_dette <- function(n_events    = 500,
                           d           = 4L,
                           n_bins      = 12L,
                           rate_scale  = 3,
                           rate_offset = 0,
                           holdout_frac = 0.10,
                           seed        = 0L) {
  d <- as.integer(d)
  stopifnot(d >= 1L, d <= 8L)

  ## Calibrate: sample reference points to find raw range (deterministic offset)
  set.seed(as.integer(seed) + 99L)
  X_ref <- matrix(stats::runif(5000L * d), 5000L, d)
  X8_ref <- matrix(0.5, 5000L, 8L)
  X8_ref[, seq_len(d)] <- X_ref
  raw_ref <- .f_dette_internal(X8_ref)
  raw_min <- min(raw_ref); raw_max <- max(raw_ref)

  ## True log-intensity on [0,1]^d (remaining dims fixed at 0.5)
  ## eta(x) = rate_offset + rate_scale * (f(x) - f_min) / (f_max - f_min)
  ## So eta ranges from rate_offset to rate_offset + rate_scale.
  ## lambda_max = exp(rate_offset + rate_scale) -- but we express per-cell.
  ## Per-cell mean = cell_vol * lambda(x) where cell_vol = (1/n_bins)^d.
  ## To get ~n_events total: n_events = cell_vol * sum lambda(x_i) over grid.
  ## We calibrate rate_offset so that the total expected count = n_events.
  eta_true_fn <- function(X) {
    X <- as.matrix(X)
    X8 <- matrix(0.5, nrow(X), 8L)
    X8[, seq_len(d)] <- X[, seq_len(d)]
    raw <- .f_dette_internal(X8)
    rate_offset + rate_scale * (raw - raw_min) / (raw_max - raw_min)
  }

  ## Calibrate rate_offset to hit n_events total on the grid
  ## E[N] = cell_vol * sum_{all cells} exp(eta(centre_i))
  cell_vol <- (1 / n_bins)^d
  ## Approximate sum using reference sample
  eta_ref  <- eta_true_fn(X_ref)
  ## E[N] = cell_vol * n_bins^d * mean(exp(eta)) = mean(exp(eta))
  ## (since cell_vol * n_bins^d = 1 and mean over uniform approx integral)
  mean_lambda <- mean(exp(eta_ref))
  ## Adjust rate_offset so mean_lambda = n_events / 1 (unit box, intensity per unit vol)
  ## target: mean(exp(eta + delta)) = n_events  => delta = log(n_events) - log(mean_lambda)
  delta <- log(n_events) - log(mean_lambda)
  rate_offset_adj <- rate_offset + delta

  eta_true_fn_adj <- function(X) {
    eta_true_fn(X) + delta
  }

  ## Inhomogeneous Poisson process by thinning on [0,1]^d
  ## Expected events = integral of lambda(x) dx = mean(exp(eta)) over [0,1]^d
  ## By calibration above, mean(exp(eta_adj)) = n_events.
  ## lambda_max = exp(rate_offset_adj + rate_scale)  [eta is capped at this]
  ## We draw a homogeneous PP with rate lambda_max, then thin.
  ## E[N_homogeneous] = lambda_max; keep each with prob lambda(x)/lambda_max.
  ## E[N_thinned] = integral lambda(x) dx = n_events.
  set.seed(as.integer(seed))
  lambda_max_est <- exp(rate_offset_adj + rate_scale)
  ## Draw number of homogeneous proposals from Poisson(lambda_max)
  n_propose <- stats::rpois(1L, lambda_max_est)
  n_propose <- max(n_propose, ceiling(lambda_max_est * 0.5))
  X_prop <- matrix(stats::runif(n_propose * d), n_propose, d)
  eta_prop <- eta_true_fn_adj(X_prop)
  lambda_prop <- exp(eta_prop)
  accept <- stats::runif(n_propose) < lambda_prop / lambda_max_est
  X_pts <- X_prop[accept, , drop = FALSE]
  cat(sprintf("Thinning: %d events accepted (approx target %d)\n",
              nrow(X_pts), n_events))

  ## Bin onto regular grid
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  centres <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  cell_idx <- matrix(0L, nrow(X_pts), d)
  for (j in seq_len(d)) {
    cell_idx[, j] <- findInterval(X_pts[, j], breaks, rightmost.closed = TRUE)
    cell_idx[, j] <- pmax(1L, pmin(n_bins, cell_idx[, j]))
  }
  ## Convert multi-index to scalar key and tally
  strides <- cumprod(c(1L, rep(n_bins, d - 1L)))
  keys    <- as.integer(cell_idx %*% strides - sum(strides) + 1L)
  tbl     <- tabulate(keys, nbins = prod(rep(n_bins, d)))
  occupied <- which(tbl > 0L)
  y_occ   <- tbl[occupied]

  ## Reconstruct centres for occupied cells
  multi_idx <- function(linear_idx) {
    idx <- matrix(0L, length(linear_idx), d)
    rem <- linear_idx - 1L
    for (j in seq_len(d)) {
      idx[, j] <- (rem %% n_bins) + 1L
      rem <- rem %/% n_bins
    }
    idx
  }
  cidx  <- multi_idx(occupied)
  X_occ <- matrix(centres[cidx], ncol = d)

  cat(sprintf("Occupied cells: %d / %d (%.4f%%)\n",
              length(occupied),
              n_bins^d,
              100 * length(occupied) / n_bins^d))

  ## Train / test split
  n_occ  <- length(occupied)
  n_test <- max(1L, round(holdout_frac * n_occ))
  test_idx  <- sort(sample.int(n_occ, n_test))
  train_idx <- setdiff(seq_len(n_occ), test_idx)

  list(
    X_train      = X_occ[train_idx, , drop = FALSE],
    y_train      = y_occ[train_idx],
    X_test       = X_occ[test_idx,  , drop = FALSE],
    y_test       = y_occ[test_idx],
    centres      = centres,
    breaks       = breaks,
    cell_X_occ   = X_occ,
    cell_y_occ   = y_occ,
    train_idx    = train_idx,
    test_idx     = test_idx,
    eta_true_fn  = eta_true_fn_adj,   # calibrated: integral = n_events
    d            = d,
    n_bins       = n_bins,
    n_events     = nrow(X_pts),
    seed         = seed,
    rate_offset_adj = rate_offset_adj,
    delta           = delta
  )
}

#' @export
print.ttps_sim <- function(x, ...) {
  cat(sprintf("TTPsplines simulation: %s\n", x$name %||% "custom"))
  cat(sprintf("  n=%d, d=%d, family=%s", x$n, x$d, x$family))
  if (identical(x$family, "gaussian")) {
    cat(sprintf(", sigma=%.3g", x$sigma))
  }
  cat("\n")
  invisible(x)
}

#' Convert a simulation list to a data frame.
#'
#' @param x A `"ttps_sim"` object.
#' @param row.names Optional row names.
#' @param optional Unused.
#' @param ... Unused.
#' @return A `data.frame` with columns `x1..xd`, `y`, `f`.
#' @export
as.data.frame.ttps_sim <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- data.frame(x$X, check.names = FALSE)
  names(df) <- colnames(x$X) %||% paste0("x", seq_len(ncol(x$X)))
  df$y <- as.numeric(x$y)
  df$f <- as.numeric(x$f)
  if (!is.null(row.names)) rownames(df) <- row.names
  df
}
