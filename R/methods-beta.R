#' @title Confidence intervals based on the beta distribution
#' @aliases beta-method
#' @description Confidence intervals for directly standardized rates 
#' based on the beta distribution.
#' @details Uses the methods proposed by Tiwari et al (2006) and further explored by
#' Ng et al (2008). 
#' @return a vector with the lower and upper bound of the confidence interval.
#' The estimate of the directly standardised rate and the level of confidence are 
#' returned as attributes to this vector.
#' @param x a vector of stratum-specific counts of events
#' @param w a vector of stratum-specific weights
#' @param level confidence level for the returned confidence interval
#' @param type type of modification for the beta confidence interval
#' @references 
#' Ng, Filardo, & Zheng (2008). 
#' 'Confidence interval estimating procedures for 
#' standardized incidence rates.' *Computational Statistics and Data Analysis* **52**: 3501-3516.
#' 
#' Tiwari, Clegg, & Zou (2006). 
#' 'Efficient interval estimation for age-adjusted cancer rates.'
#'  *Statistical Methods in Medical Research* **15**: 547-569.
#' 
#' @importFrom stats qbeta
#' @export
ci.beta <- function(x, w, level, 
                     type = c("tcz", "cc")){
  type = match.arg(type)

  y = sum(x*w)
  v = sum(x*w^2)
  f.a = function(y,v){y*(y*(1-y)/v - 1)}
  f.b = function(y,v){(1-y)*(y*(1-y)/v - 1)}
  a = f.a(y,v)
  b = f.b(y,v)
  yst = y + mean(w)
  vst = v + mean(w^2)
  ci <- switch(type,
    "tcz" = stats::qbeta(alpha(level), f.a(y,v), f.b(y,v)),
    "cc"  = stats::qbeta(alpha(level), f.a(yst,vst), f.b(yst,vst)))
  attr(ci, "estimate") <- y
  attr(ci, "level") <- level
  attr(ci, "method.arg") <- type
  ci
}
