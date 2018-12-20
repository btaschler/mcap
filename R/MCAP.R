
MCAP <- function(xx, k, projection = 'PCA',
                 true_labels = NULL, parallel = FALSE, verbose = FALSE){
  #' Model based clustering via adaptive (linear) projections
  #'
  #' Model based clustering using full variance Gaussian mixtures in a lower
  #' dimensional projected space obtained via adaptive (linear) projections.
  #' Projection variants include PCA-based and random projection. 
  #' The wrapper performs multiple restarts in case true labels are provided 
  #' (in order to optimise cluster assignments w.r.t. the Rand index).
  #' 
  #' @param xx The data matrix (n x p).
  #' @param k The number of clusters.
  #' @param true_labels Vector of true cluster assignments (when provided, it is 
  #'                    used to compute the Rand index). 
  #' @param parallel Logical, when true: perform line search over projection 
  #'                 dimension in parallel. 
  #' @param verbose Logical, when true: print some progress information. 
  #' 
  #' @return @param model_fit Model fit (GMM output of `mixglasso()`).
  #' @return @param aRI Adjusted Rand index (when `true_labels` are provided).
  #' @export
  
  ## input checks
  if(!is.null(true_labels)){ 
    stopifnot(length(true_labels) == nrow(xx)) 
    stopifnot(length(unique(true_labels)) == k)
  }
  
  
  
  return(NA)
}