
Precision2PartialCor <- function(xx){
  #' Compute partial correlation matrix from inverse covariance matrix
  #'
  #' This function takes an inverse covariance (or precision) matrix as input 
  #' and computes the partial correlation matrix.
  #' 
  #' @param xx The inverse covariance matrix.
  #' 
  #' @return \item{pcor}{ The partial correlation matrix}. 
  #' @examples 
  #'   A <- matrix(rnorm(100),10)
  #'   Precision2PartialCor(xx = A%*%t(A))
  #' @export
  
  x <- as.matrix(xx)
  p <- dim(xx)[1]
  d <- diag(xx)
  
  diag_prec <- diag(diag(xx)^(-.5))
  pcor <- diag(2,p) - diag_prec %*% xx %*% diag_prec
  #note: diag(2,p) is only needed to turn diag entries of pcor from -1 into +1
  
  return(pcor)
}