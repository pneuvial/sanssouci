% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SansSouciStruct-class.R
\name{fit.SansSouciStruct}
\alias{fit.SansSouciStruct}
\title{Fit SansSouciStruct object}
\usage{
\method{fit}{SansSouciStruct}(
  object,
  alpha,
  p.values,
  family = c("DKWM", "HB", "trivial", "Simes", "Oracle"),
  flavor = c("tree", "partition"),
  ...
)
}
\arguments{
\item{object}{An object of class \code{SansSouciStruct}}

\item{alpha}{Target risk (JER) level}

\item{p.values}{A vector of length \code{nHyp(object)},}

\item{family}{A character value describing how the number of true nulls in a set is estimated. Can be either:
\itemize{
\item "DKWM": estimation by the Dvoretzky-Kiefer-Wolfowitz-Massart inequality (related to the Storey estimator of the proportion of true nulls), valid for independent p-values
\item "HB": estimation by the Holm-Bonferroni method, always valid
\item "trivial": dummy estimation as the the size of the set
\item "Simes": estimation via the Simes inequality, valid for positively-dependent (PRDS) p-values
\item "Oracle": true number of true null hypotheses Truth" must be available in \code{object$input$truth}
}}

\item{flavor}{A character value which can be
\itemize{
\item "tree" (default value): the reference family is the entire tree structure
\item "partition": the reference family is the partition corresponding to the leaves of the tree
}}

\item{...}{Not used}
}
\value{
A 'fitted' object of class 'SansSouciStruct'. It is a list of three elements
\itemize{
\item input: see \link{SansSouciStruct}
\item param: the input parameters, given as a list
\item output: a list of two elements
\itemize{
\item p.values: the input argument 'p.values'
\item ZL: the output of the "zeta function" associated to the input parameter 'family', see e.g. \link{zeta.DKWM}
}
}
}
\description{
Fit SansSouciStruct object
}
\details{
In the particular case where \code{family=="Simes"} or \code{family=="Oracle"}, the return value is actually of class \code{SansSouci} and not \code{SansSouciStruct}
}
\examples{
s <- 100
q <- 7
m <- s*2^q
obj <- SansSouciDyadic(m, leaf_size = s, direction = "top-down")

mu <- gen.mu.leaves(m = m, K1 = 8, d = 0.9, grouped = TRUE, 
  setting = "const", barmu = 3, leaf_list = obj$input$leaves)
pvalues <- gen.p.values(m = m, mu = mu, rho = 0)

alpha <- 0.05
S1 <- which(mu != 0)

res_DKWM <- fit(obj, alpha, pvalues, "DKWM")
predict(res_DKWM, S = S1, what = "FP")

res_Simes <- fit(obj, alpha, pvalues, "Simes")
predict(res_Simes, S = S1, what = "FP")

}
\references{
Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.

Dvoretzky, A., Kiefer, J., and Wolfowitz, J. (1956). Asymptotic minimax character of the sample distribution function and of the classical multinomial estimator. The Annals of Mathematical Statistics, pages 642-669.

Holm, S. A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6 (1979), pp. 65-70.

Massart, P. (1990). The tight constant in the Dvoretzky-Kiefer-Wolfowitz inequality. The Annals of Probability, pages 1269-1283.
}
\seealso{
zeta.DKWM zeta.HB, zeta.tricial
}
