% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-normal.R
\name{ci.asymptotic}
\alias{ci.asymptotic}
\title{Confidence intervals using Normal approximation}
\usage{
ci.asymptotic(x, w, level, trans = c("none", "log", "cube.root", "skew"))
}
\arguments{
\item{x}{a vector of counts}

\item{w}{a vector of weights}

\item{level}{the level of confidence}

\item{trans}{transformation to apply}
}
\value{
a vector with the lower and upper bound of the confidence
interval.The estimate of the directly standardised rate and the
of confidence are returned as attributes to this vector
}
\description{
Confidence intervals for directly standardized rate
estimate. A common approach for constructing confidence intervals
around an MLE is to use a normal approximation of the MLE or transformed
MLE.
}
\details{
Four different transformations are implented -
following Ng et al (2008).
\itemize{
\item none (no transformation)
\item log
\item cubic root
\item Edgeworth correction for skewness
}
}
\references{
Ng, Filardo, & Zheng (2008). 'Confidence interval estimating
procedures for standardized incidence rates.'
\emph{Computational Statistics and Data Analysis} \strong{52} 3501--3516.
\doi{doi:10.1016/j.csda.2007.11.004}
}
