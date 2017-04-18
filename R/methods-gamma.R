#' @title Confidence intervals using Gamma methods  
#' @aliases gamma-method
#' @seealso \code{\link[asht]{wspoissonTest}}
#' @description Confidence intervals for directly standardized 
#' rates based on the gamma distribution.
#' @details Calls `asht::wspoissonTest` with the appropriate 
#' type modification for the gamma confidence interval (`wmtype`) and
#' `'mid.p'` value. Note that the `'mid.p'` method uses `wmtype = 'max'`.
#' @param x a vector of stratum-specific counts of events
#' @param w a vector of stratum-specific weights
#' @param level confidence level for the returned confidence interval
#' @param type type of modification for the gamma confidence interval
#' @param ... passed to `asht::wspoissonTest` when type = `'midp'`
#' @references 
#' Fay  MP (2017). '`asht`: Applied Statistical Hypothesis Tests'. 
#' *`R` package version 0.9.*  \CRANpkg{ahst}
#' 
#' Fay & Feuer (1997). 'Confidence intervals for directly 
#' standardized rates: a method based on the gamma distribution.' 
#' *Statistics in Medicine*. **16**: 791--801.
#' \url{https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7<791::AID-SIM500>3.0.CO;2-\%23} 
#' 
#' Fay & Kim (2017). 'Confidence intervals for directly 
#' standardized rates using mid-p gamma intervals.' 
#' *Biometrical Journal*. \doi{doi:10.1002/bimj.201600111}.
#' 
#' Ng, Filardo, & Zheng (2008). 'Confidence interval estimating procedures for 
#' standardized incidence rates.' *Computational Statistics and Data Analysis* 
#' **52** 3501--3516. \doi{doi:10.1016/j.csda.2007.11.004}
#'  
#' Tiwari, Clegg, & Zou (2006). 'Efficient interval estimation for age-adjusted 
#' cancer rates.'  *Statistical Methods in Medical Research* **15**: 547--569.
#' \doi{doi:10.1177/0962280206070621}
#' 
#' @importFrom asht wspoissonTest
#' @export
ci.gamma <- function(x, w, level, 
                     type = c("max", "midp", "mean", "minmaxavg", "tcz"),
                     ... ){
  type = match.arg(type)
  if(type == 'midp'){
    ci <- asht::wspoissonTest(x, w, conf.level = level, midp = TRUE, ...)[['conf.int']]
  } else {
    ci <- asht::wspoissonTest(x, w, conf.level = level, wmtype = type)[['conf.int']]
  }
  attr(ci, "estimate") <- sum(x*w)
  attr(ci, "level") <- level
  attr(ci, "method.arg") <- type
  ci
}
