
Rand1 <- function(n){
  #' Return a random sign vector
  #'
  #' Create a random sign vector of length \code{n} (iid ~ +1/-1).
  #' 
  #' @param n Number of elements
  #' 
  #' @return Sign vector of length \code{n}.
  #' @export
  
  return((round(stats::runif(n)) * 2) - 1)
}