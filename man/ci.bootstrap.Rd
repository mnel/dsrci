% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-normal.R
\name{ci.bootstrap}
\alias{ci.bootstrap}
\alias{ci.abc}
\title{Approximate Bootstrap Confidence Interval (ABC)}
\usage{
ci.bootstrap(x, w, level, ...)
}
\arguments{
\item{x}{a vector of counts}

\item{w}{a vector of weights}

\item{level}{the level of confidence}

\item{...}{currently ignored}
}
\value{
a vector with the lower and upper bound of the
confidence interval.The estimate of the directly standardised
rate and the level of confidence are returned as attributes
this vector
}
\description{
Confidence intervals for directly standardized rates
using the approximate bootstrap method  derived by Swift (1995)
}
\references{
Swift, MB (1995). 'Simple confidence intervals for
standardized rates based on the approximate bootstrap method',
\emph{Statistics in Medicine}, \strong{14}, 1875--1888.
\doi{doi:10.1002/sim.4780141704}.
}
