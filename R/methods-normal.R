#' @title Confidence intervals using Normal approximation
#' @name ci.asymptotic
#' 
#' @description Confidence intervals for directly standardized rate
#' estimate. A common approach for constructing confidence intervals 
#' around an MLE is to use a normal approximation of the MLE or transformed
#' MLE. 
#' @details Four different transformations are implented - following Ng et al (2008).
#' - none (no transformation) 
#' - log
#' - cubic root
#' - Edgeworth correction for skewness
#' @return a vector with the lower and upper bound of the confidence interval.
#' The estimate of the directly standardised rate and the level of confidence are 
#' returned as attributes to this vector 
#' @references 
#' Ng, Filardo, & Zheng (2008). 'Confidence interval estimating procedures for standardized incidence rates.' *
#' Computational Statistics and Data Analysis* **52** 3501-3516.
#' @param x a vector of counts
#' @param w a vector of weights
#' @param level the level of confidence
#' @param trans transformation to apply
#' @importFrom stats qnorm qlnorm
#' @export
ci.asymptotic <- function(x, w, level, 
    trans = c('none','log','cube.root','skew')){
  # match args
  trans <- match.arg(trans)
  # level to alpha
  alph <- alpha(level)
  z <- stats::qnorm(alph)
  # recalculate estimate
  y <- sum(x*w)
  v <- sum(x*w^2)
  # apply appropriate transformation
  ci <- switch(trans,
    "none" = y + z*sqrt(v),
    "log" = stats::qlnorm(alph, log(y), sqrt(v/y^2)),
    "cube.root" =  (y^(1/3) + z*sqrt(v/(9*(y^(4/3)))))^3,
    "skew" =   y + (z - edgeworth.skew(x,w,rev(z)))*sqrt(v))
  attr(ci, "estimate") <- y
  attr(ci, "level") <- level
  attr(ci, "method.arg") <- trans
  ci
}

#' @title Approximate Bootstrap Confidence Interval (ABC)
#' @aliases ci.abc
#' 
#' @description Confidence intervals for directly standardized rates 
#' using the approximate bootstrap method  derived by Swift (1995)
#' @return a vector with the lower and upper bound of the confidence interval.
#' The estimate of the directly standardised rate and the level of confidence are 
#' returned as attributes to this vector 
#' @references 
#' Swift, MB (1995). 'Simple confidence intervals for standardized rates based on 
#' the approximate bootstrap method', *Statistics in Medicine*, **14**, 1875—1888.
#' @param x a vector of counts
#' @param w a vector of weights
#' @param level the level of confidence
#' @param ... currently ignored
#' @export
ci.bootstrap <- function(x, w, level, ...){
  z <- stats::qnorm(alpha(level))
  y <- sum(w*x)
  v <- sum(w^2*x)
  z0 <- a <- sum((w^3)*x)/(6*v^(3/2))
  numerator <- (z0 + z)
  ci <- y+sqrt(v)*numerator/ ((1-a*numerator)^2)
  attr(ci, "estimate") <- y
  attr(ci, "level") <- level
  attr(ci, "method.arg") <- NA_character_
  ci
}

