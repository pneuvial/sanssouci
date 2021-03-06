% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confidenceEnvelope.R
\name{confCurveFromFam}
\alias{confCurveFromFam}
\title{Confidence bounds on the true/false positives among most significant items}
\usage{
confCurveFromFam(
  p.values,
  refFamily,
  param,
  K = length(p.values),
  what = c("TP", "FDP")
)
}
\arguments{
\item{p.values}{A vector containing the p-values for all m hypotheses, sorted
increasingly}

\item{refFamily}{A character value, the reference family to be used. Should
be either "Simes" (or equivalenlty, "Linear"), "Beta", or "Oracle".}

\item{param}{A numeric value or vector of parameters for the reference family.}

\item{K}{For JER control over \code{1:K}, ie joint control of all
k-FWER, k <= K. Defaults to m.}

\item{what}{A character vector, the names of the post hoc bounds to be
computed, among:
\itemize{
\item FP: Upper bound on the number of false positives in the 'x' most significant items
\item TP: Lower bound on the number of true positives in the 'x' most significant items
\item FDP: Upper bound on the proportion of false positives in the 'x' most significant items
\item TP: Lower bound on the proportion of true positives in the 'x' most significant items.
}

Defaults to \code{c("TP", "FDP")}}
}
\value{
A \code{data.frame} with \code{m} rows and 5 columns:
\itemize{
\item x: Number of most significant items selected
\item family: Matches input argument \code{refFamily}
\item param: Matches argument \code{param}
\item procedure: Label for the procedure, typically of the form 'refFamily(param)'
\item bound: Value of the post hoc bound
\item stat: Type of post hoc bound, as specified by argument \code{bound}
}
}
\description{
Confidence bounds on the true/false positives among most significant items
}
\details{
\code{param} should be a numeric value unless \code{refFamily == "Oracle"}. In the latter case, `param`` should be a boolean vector of
length m indicating whether each null hypothesis is true or false.
}
\examples{

# Generate Gaussian data and perform multiple tests
sim <- gaussianSamples(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
rwt <- rowWelchTests(sim$X, sim$categ, alternative = "greater")

# calculate, print, and plot confidence bound
cb <- confCurveFromFam(rwt$p.value, refFamily = "Simes", param = 0.1)
head(cb)
plotConfCurve(cb, xmax = 200) 

}
\author{
Gilles Blanchard, Pierre Neuvial and Etienne Roquain
}
