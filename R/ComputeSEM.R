

ComputeSEM <- function(x){
  #' Compute standard error of the mean for an input vector
  #' 
  #' @param x Input vector (any NA will be ignored).
  #' 
  #' @return Standard error of the mean (SEM) of \code{x}. 
  #' @examples 
  #'   ## random vector:
  #'   ComputeSEM(rnorm(100))
  #'   
  #'   ## with missing entries:
  #'   x <- runif(1000)*5
  #'   x[sample(1000, 50)] <- NA
  #'   ComputeSEM(x)
  #' @export
  
  return(stats::sd(x, na.rm=TRUE) / sqrt(length(x)-sum(is.na(x))))
}
