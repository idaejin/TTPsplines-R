% Please edit documentation in R/families.R
\name{poisson_intensity_ll}
\alias{poisson_intensity_ll}
\title{Held-out log-likelihood of a binned Poisson intensity}
\usage{
poisson_intensity_ll(y_te, eta_te, eta_U, n_cells, delta = 1, p_te = 1)
}
\arguments{
\item{y_te, eta_te}{Held-out counts and log-intensity at those cells.}

\item{eta_U}{Log-intensity on the domain sample.}

\item{n_cells}{Number of domain cells \eqn{M}.}

\item{delta}{Cell volume \eqn{\Delta}. Use \code{1/M} on the unit cube, and
\code{1} when \eqn{\exp(\eta)} is already a cell mean.}

\item{p_te}{Weight on the integral, usually \code{1 - p_tr}.}
}
\value{
The log-likelihood
\eqn{\sum y^{\mathrm{te}}\eta - p_{\mathrm{te}}\sum_j\mu_j},
up to a term that does not depend on \eqn{\eta}.
}
\description{
Cell means are \eqn{\mu_j = \Delta\exp(\eta(x_j))}. The integrated intensity
is \eqn{M\Delta} times the mean of \eqn{\exp(\eta)} on \code{eta_U}: every
cell, or a uniform sample of cells. That sample is not the empty-cell
quadrature used inside a Berman--Turner fit.
}
