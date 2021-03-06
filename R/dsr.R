#' @title Calculate directly standardized rate
#' @name dsr
#' @description Calculate a directly standardised rate with 
#' confidence interval based on a specified method.
#' @details Five groupds of methods can be specified using 
#' `'ci.method'`, with variations on each depending on the method. 
#'  five groups are:
#' 
#' \tabular{ll}{
#' \code{'asymptotic'} \tab Using the normal approximation of the 
#' MLE distribution (or transformed MLE distribution). 
#' See \code{\link{ci.asymptotic}} for more details and  
#'  the currently implemented transformations. \cr
#' \code{'moments'} \tab Moment matching based on Dobson et al (1991). 
#' A variety methods for constructing the confidence interval on the 
#' unweighted sum of `x`. - See \code{\link{ci.moments}} for more details. \cr
#' \code{'gamma'} \tab Based on the gamma distribution (Fay & Feuer 1997). 
#' See \code{\link{ci.gamma}} for more details and modifications implemented. \cr
#' \code{'beta'} \tab Based on the beta distribution as proposed by 
#' Tiwari et al (2006). See \code{\link{ci.beta}} for details and 
#' the modifications implemented. \cr
#' \code{'bootstrap'} \tab Appromiximate Bootstrap Confidence 
#'  proposed by Swift (1995). \cr
#' }
#' @return `dsr` returns an object of class "`dsr`" . The function `confint` is used
#' to return a confidence interval using a specified method
#' 
#' An object of class "`dsr`" is a list containing the following components:
#' \tabular{ll}{
#' \code{estimate} \tab the estimate of the directly standardised rate \cr
#' \code{lower} \tab lower bound of the confidence interval \cr
#' \code{upper} \tab upper bound of the confidence interval \cr
#' \code{level} \tab the level of confidence \cr
#' \code{ci.method} \tab method used to calcalate the confidence interval \cr
#' \code{method.arg} \tab additional argument passed to ci method function \cr
#' \code{call} \tab the matched call \cr
#' \code{mult} \tab The multiplicative factor 
#' to scale the final estimate \cr
#' \code{strata} \tab number of strata or summands \cr
#' }
#
#' @references 
#' Dobson, AJ, Kuulasmaa, K, Eberle, E and Scherer, J (1991) 
#' 'Confidence intervals for weighted sums of Poisson parameters', 
#' *Statistics in Medicine*, **10**: 457--462.
#' \doi{doi:10.1002/sim.4780100317}
#' 
#' Swift, MB (1995). 'Simple confidence intervals for 
#' standardized rates based on the approximate bootstrap method', 
#' *Statistics in Medicine*, **14**, 1875--1888.
#' \doi{doi:10.1002/sim.4780141704}. 
#' 
#' Fay & Feuer (1997). 'Confidence intervals for directly 
#' standardized rates: a method based on the gamma distribution.
#' Statistics in Medicine*. **16**: 791--801.
#' \url{https://doi.org10.1002/(SICI)1097-0258(19970415)16:7<791::AID-SIM500>3.0.CO;2-\%23}
#' 
#' Tiwari, Clegg, & Zou (2006). 'Efficient interval estimation 
#' for age-adjusted cancer rates.' 
#' *Statistical Methods in Medical Research* **15**: 547--569. 
#' \doi{doi:10.1177/0962280206070621}
#' 
#' Ng, Filardo, & Zheng (2008). 'Confidence interval estimating 
#' procedures for standardized incidence rates.' 
#' *Computational Statistics and Data Analysis* **52** 3501-3516.
#' \doi{doi:10.1016/j.csda.2007.11.004}
#' 
#' @param x a vector of strata-specific counts 
#' @param n a vector of strata-specific time bases for counts
#' @param w a vector of strata-specific weights (or standard populations)
#' @param mult a factor to multiply the estimate to give rates per `mult`
#' @param ci.method  method used to calculate the confidence interval. See details
#' @param level the confidence level required
#' @param ... Further arguments passed to the confidence interval function.
#' @export
dsr <- function(x, n, w, 
                ci.method = c("asymptotic", "moments", "gamma", 
                              "beta", "bootstrap"),
                level = 0.95, mult = 1000, ...){
  # formula / stats::model.frame method from lm
  cl <- match.call()
  # Check input
   if((is.null(x)) | !is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }
  if((is.null(n)) | !is.numeric(n)) {
    stop("'n' must be a numeric vector")
  }
  if(is.null(w) | !is.numeric(w)){
    stop("'w' must be a numeric vector")
  }
  # valid level
  if(!(is.numeric(level))||((level <0.5) |(level>=1))){
    stop("'level' must be a numeric value in [0.5,1)")
  }
  # check for missing values
  errorNA(x)
  errorNA(n)
  errorNA(w)
  # lengths
  if(length(unique(lengths(list(x,n,w))))!=1){
    stop("'x','n' and 'w' must have the same lengths")
  }
  
  # method argument
  ci.method <- match.arg(ci.method)
  # recalculate weights
  w = w /(n * sum(w))
  # 
  
  ci <- switch(ci.method,
    "asymptotic" = ci.asymptotic(x, w, level, ...),
    "moments" = ci.moments(x, w, level, ...),
    "gamma" = ci.gamma(x, w, level, ...),
    "beta" = ci.beta(x, w, level, ...),
    "bootstrap" = ci.bootstrap(x, w, level, ...))
  
  data.frame(
    estimate = attr(ci, "estimate") * mult, 
    lower = ci[1] * mult, 
    upper = ci[2] * mult,
    level = level,
    ci.method = ci.method,
    method.arg = attr(ci, "method.arg"),
    mult = mult,
    strata = length(x))  

}

