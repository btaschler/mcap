
Rank2Normality <- function (xx){
  #' Transform data matrix to normality based on rank
  #'
  #' This function performs a rank to normality transformation of the columns
  #' of an input matrix X.
  #' 
  #' @seealso \code{\link[stats]{qnorm}}
  #' 
  #' @param xx The data matrix.
  #' 
  #' @return tt Data matrix transformed to normality.
  #' @export
  
  # number of samples
  n <- nrow(xx)
  
  # rank transformation for each column
  rr <- apply(xx, 2, rank)
  
  # inverse Normal cdf
  tt <- stats::qnorm(rr/(n+1))
  
  return(tt)
}