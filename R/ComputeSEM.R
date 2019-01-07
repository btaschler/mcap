

ComputeSEM <- function(x){
  #' Compute standard error of the mean for an input vector
  #' 
  #' @param x Input vector (any NA will be ignored).
  #' 
  #' @return Standard error of the mean (SEM) of \code{x}. 
  #' @export
  
  return(stats::sd(x, na.rm=TRUE) / sqrt(length(x)-sum(is.na(x))))
}
