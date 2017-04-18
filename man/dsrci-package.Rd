\name{dsrci-package}
\alias{dsrci-package}
\alias{dsrci}
\docType{package}
\title{\packageTitle{dsrci}}
\description{\packageDescription{dsrci}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{dsrci}
  \packageIndices{dsrci}
}
\author{
  \packageAuthor{dsrci}
  
  Maintainer: \packageMaintainer{dsrci}
}
\references{
Fay & Feuer (1997). 'Confidence intervals for directly standardized rates: a method based on the gamma distribution.' \emph{Statistics in Medicine}. \strong{16}: 791-801.

Fay & Kim (2017). 'Confidence intervals for directly standardized rates using mid-p gamma intervals.' \emph{Biometrical Journal}. \doi{doi:10.1002/bimj.201600111}.

Ng, Filardo, & Zheng (2008). 'Confidence interval estimating procedures for standardized incidence rates.' \emph{Computational Statistics and Data Analysis} \strong{52}: 3501-3516.

Tiwari, Clegg, & Zou (2006). 'Efficient interval estimation for age-adjusted cancer rates.' \emph{Statistical Methods in Medical Research} \strong{15}: 547-569.

Swift, MB (1995). 'Simple confidence intervals for standardized rates based on the approximate bootstrap method', \emph{Statistics in Medicine}, \strong{14}, 1875—1888.

Dobson, AJ, Kuulasmaa, K, Eberle, E and Scherer, J (1991)
'Confidence intervals for weighted sums of Poisson parameters',
\emph{Statistics in Medicine}, \strong{10}: 457—462.
}
\keyword{package}
\seealso{
\code{\link[asht]{wspoissonTest}}, \code{\link[dsrci]{dsr}}
}
\examples{
\dontrun{
## using the example from asht::wspoissonTest
## birth data on Down's syndrome from Michigan, 1950-1964
## see Table II  of Fay and Feuer (1997)
## xfive= counts for mothers who have had 5 or more children
## nfive and ntotal are number of live births 
xfive<-c(0,8,63,112,262,295)
nfive<-c(327,30666,123419,149919,104088,34392)
ntotal<-c(319933,931318,786511,488235,237863,61313)
## use mult =10^5 to give rates per 100,000
## gamma method of Fay and Feuer (1997) is default
dsr(xfive, nfive, ntotal, ci.method = "gamma", type = "max", mult = 1e5)
## Dobson et al (1991)
dsr(xfive, nfive, ntotal, ci.method = "moments", type = "dobson", mult = 1e5)
  }
}
