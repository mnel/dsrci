#' @title Confidence intervals for Directly Standardised Rates using moment match methods
#' @aliases ci.dobson
#' 
#' @description Confidence intervals for directly standardized rates based on
#' the approximation proposed by Dobson, Kuulasmaa, Eberle and Scherer (1991).
#' In addition to the method proposed by Dobson et al, the various approximations
#' implemented by Ng, Filardo & Zheng (2008) are implemented.
#' @return a vector with the lower and upper bound of the confidence interval.
#' The estimate of the directly standardised rate and the level of confidence are 
#' returned as attributes to this vector 
#' @param x a vector of stratum-specific counts of events
#' @param w a vector of stratum-specific weights
#' @param level confidence level for the returned confidence interval
#' @param type type of approximation for the poisson confidence interval of the unweighted sum
#' @references 
#' Dobson, AJ, Kuulasmaa, K, Eberle, E and Scherer, J (1991) 
#' 'Confidence intervals for weighted sums of Poisson parameters', 
#' *Statistics in Medicine*, **10**: 457—462.
#' 
#' Ng, Filardo, & Zheng (2008). 'Confidence interval estimating procedures for standardized incidence rates.' 
#' *Computational Statistics and Data Analysis* **52**: 3501-3516.
#' @importFrom exactci poisson.exact
#' @export
ci.moments <- function(x, w, level, type = 
 c("dobson", "boise.monson", "normal", "wilson.hilferty",
   "byar", "midp", "approx.midp", "simple.midp")){
  type = match.arg(type)
  y <- sum(x*w)
  v <- sum(w^2*x)
  X <- sum(x)
  X01 <- X+c(0,1)
  z <- stats::qnorm(alpha(level))
  dobson <- function(CI, Y = y,S = sqrt(v/X), Xt = X){Y +S*(CI - Xt)}
  ci <- dobson(
    switch(type,
    "dobson" = stats::qgamma(alpha(level),X01),
    "boise.monson" = exp(log(X) + z/sqrt(X)),
    "normal" = X + z^2/2 + z * sqrt(X + z^2/4),
    "wilson.hilferty" = X + (2*z^2 + 1)/6 + (0.5 + z *sqrt(X + (z^2+2)/18-0.5)),
    "byar" = X01*(1 - 1/(9*X01) - z/(3*sqrt(X01)))^3,
    "midp" = exactci::poisson.exact(X, midp = TRUE)[['conf.int']],
    "approx.midp" = (X+0.5)*(1 - 1/(9*sqrt(X+0.5))  + z/(3*sqrt(X+0.5))),
    "simple.midp" = (sqrt(X+0.5) + z/2)^2)
    )
  attr(ci, "estimate") <- y
  attr(ci, "level") <- level
  ci
  }