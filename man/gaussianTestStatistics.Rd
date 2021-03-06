% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussianTestStatistics.R
\name{gaussianTestStatistics}
\alias{gaussianTestStatistics}
\title{Simulate Gaussian test statistics}
\usage{
gaussianTestStatistics(
  m,
  B,
  pi0 = 1,
  SNR = 0,
  dep = c("equi", "Toeplitz"),
  param = 0
)
}
\arguments{
\item{m}{Number of hypotheses}

\item{B}{Number of simulations}

\item{pi0}{Proportion of true null hypotheses}

\item{SNR}{Signal to noise ratio. Either a numeric value (a measure of
distance between H0 and H1) or a vector of length \code{m*(1-pi0)}}

\item{dep}{A character value, the type of dependency between test statistics. Can be one of "equi" for equi-correlation, or "Toeplitz". Defaults to "equi".}

\item{param}{A numeric value defaulting to \code{0}. If \code{dep=="equi"}, \code{param} is the level of equi-correlation between pairs of variables. If \code{dep=="Toeplitz"}, the first row of the Toeplitz matrix will be \code{(1:m)^(param)}.}
}
\value{
A list with elements \describe{
\item{x}{A vector of length \eqn{m} test statistics}
\item{X0}{An \eqn{m x B} matrix of test statistics under the null
hypothesis}
\item{H}{A vector of length \eqn{m}, the status of each
hypothesis: 0 for true null hypothesis, and 1 for true alternative
hypothesis} }
}
\description{
Simulate Gaussian test statistics
}
\details{
B
}
\examples{

m <- 123
B <- 100

# independent statistics under the full null
sim <- gaussianTestStatistics(m, B)

# equi-correlated statistics under the full null
sim <- gaussianTestStatistics(m, B, dep = "equi", param = 0.2)

# equi-correlated statistics with some signal
sim <- gaussianTestStatistics(m, B, pi0 = 0.8, SNR = 1, dep = "equi", param = 0.2)

## show test statistics
stat <- sim$x
pch <- 20
colStat <- 1+sim$H
plot(stat, col=colStat, main="Test statistics", pch=pch)
legend("topleft", c("H0", "H1"), pch=pch, col=1:2)

# Toeplitz statistics with some signal
sim <- gaussianTestStatistics(m, B, pi0 = 0.8, SNR = 1, dep = "Toeplitz", param = -0.5)

}
\author{
Gilles Blanchard, Pierre Neuvial and Etienne Roquain
}
