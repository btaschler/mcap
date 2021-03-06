
CentrePerGroup <- function(xx, true_labels = numeric(nrow(xx))){
  #' Centre a data matrix to mean zero (per group)
  #'
  #' Mean centering of an input matrix. When true labels are provided, the 
  #' groups/clusters are centred individually. 
  #' 
  #' @seealso \code{\link{colMeans}}, \code{\link{scale}}
  #' 
  #' @param xx The data matrix (n x p).
  #' @param true_labels Vector of true cluster assignments. Default: 0 (all 
  #'                    elements are treated as belonging to one group).
  #' 
  #' @return \item{xx_centred}{ Mean centred data matrix.}
  #' @examples 
  #'   ## 10x5 matrix, 1 group:
  #'   CentrePerGroup(xx=matrix(rnorm(50),10,5))
  #'   
  #'   ## 2 groups:
  #'   CentrePerGroup(xx=matrix(rnorm(50),10,5), true_labels=round(runif(10)))
  #'   
  #'   ## 10 groups:
  #'   CentrePerGroup(xx=matrix(rnorm(5000),100,50), true_labels=round(runif(100)*10))
  #' @export
  
  # preliminaries ------------------------------------------------------------
  n <- nrow(xx)
  p <- ncol(xx)
  xx_centred <- matrix(numeric(), n, p)
  
  ## get group labels
  labels_names <- unique(true_labels)
  k <- length(labels_names)
  
  ## mean centre each group separately
  for(j in seq(k)){
    curr_k_idx <- which(true_labels == labels_names[j])
    curr_k_data <- xx[curr_k_idx,]
    xx_centred[curr_k_idx,] <- curr_k_data - 
      t(matrix(rep(colMeans(curr_k_data), nrow(curr_k_data)), 
               ncol(curr_k_data), nrow(curr_k_data)))        
  }    
  
  return(xx_centred)
}