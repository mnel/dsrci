# Some utility functions
alpha <- function(conf.level){
  c(1 - conf.level, 1 + conf.level) / 2
}

errorNA <- function(var){
  deparse(substitute(var))
  not.ok <- anyNA(var)
  if(not.ok) stop("missing values in '",deparse(substitute(var)),"'") else (invisible(NULL))
}

# Extract relevant infomration from 
# htest object
to.list <- function(htest){
  list(
    estimate = unname(htest[['estimate']]), 
    lower =  unname(htest[['conf.int']][1]),
    upper = unname(htest[['conf.int']][2]))
  
}

# g
`edgeworth.skew` <- function(x,w,t){
  v <- sum(x*w^2)
  b1 <- (sum(x*w^3) / (v^(3/2)))^2
  b2 <- sum((x*w^4)*(1+3*x))/(v^2)
  sqrt(b1)/6 + ((b2-3)/8 +15*b1/72) - 
    (sqrt(b1)/6)*t^2 - ((b2-3)/24 + 4*b1/36)*t^3 -
    (b1/72)*t^5

}