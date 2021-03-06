% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posthoc-bounds.R
\name{curveMaxFP}
\alias{curveMaxFP}
\title{Upper bound for the number of false discoveries among most significant items}
\usage{
curveMaxFP(p.values, thr, flavor = c("BNR2016", "Mein2006", "BNR2014"))
}
\arguments{
\item{p.values}{A vector containing all \eqn{m} p-values, sorted non-decreasingly}

\item{thr}{A vector of \eqn{K} JER-controlling thresholds, sorted non-decreasingly}

\item{flavor}{The algorithm to compute the bound 'BNR2014' and
'BNR2016' give identical results. Both should be slightly better
than 'Mein2006' (example situation?). Flavor 'BNR2016' has a
linear time complexity, hence it is much faster than 'Mein2006'
and much much faster than 'BNR2014'.}
}
\value{
A vector of size \eqn{m} giving an joint upper confidence bound on
the number of false discoveries among the \eqn{k} most significant items
for all \eqn{k \in \{1\ldots m\}}.
}
\description{
Upper bound for the number of false discoveries among most significant items
}
\author{
Gilles Blanchard, Pierre Neuvial and Etienne Roquain
}
