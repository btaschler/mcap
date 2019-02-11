
KMwrapper <- function(xx, k, true_labels = NULL, verbose = FALSE){
  #' Wrapper to do K-means clustering
  #'
  #' Standard K-means clustering. 
  #' The wrapper performs multiple restarts in case true labels are provided 
  #' (in order to optimise cluster assignments w.r.t. the Rand index).
  #' 
  #' @seealso \code{\link[stats]{kmeans}}, 
  #'          \code{\link[mclust]{adjustedRandIndex}}
  #' 
  #' @param xx The data matrix (n x p).
  #' @param k The number of clusters.
  #' @param true_labels Vector of true cluster assignments (when provided, it is 
  #'                    used to compute the Rand index). 
  #' @param verbose Logical, when true: print progress information. 
  #' 
  #' @return \item{model_fit}{ Model fit (output of \code{\link[stats]{kmeans}}).}
  #'         \item{aRI}{ Adjusted Rand index (when \code{true_labels} is provided).}
  #' @examples 
  #'   KMwrapper(xx = matrix(rnorm(500),50,10), k = 2)
  #' @export
  
  ## input checks
  if(!is.null(true_labels)){ 
    stopifnot(length(true_labels) == nrow(xx)) 
    stopifnot(length(unique(true_labels)) == k)
  }
  
  ## perform K-means clustering
  if(length(true_labels)>0){
    n_repeats <- 5
    aRI_arr <- numeric(n_repeats)
    mod_fit_arr <- vector('list', n_repeats)
    
    for(i in seq(n_repeats)){
      mod_fit_arr[[i]] <- stats::kmeans(xx, centers = k, iter.max = 50, nstart = 20)
      aRI_arr[i] <- mclust::adjustedRandIndex(true_labels, mod_fit_arr[[i]]$cluster)  
    }
    max_idx <- which(aRI_arr == max(aRI_arr))[1]
    mod_fit <- mod_fit_arr[[max_idx]]
    aRI <- mean(aRI_arr, na.rm = TRUE)  
    
  }else{  #no true labels provided
    mod_fit <- stats::kmeans(xx, centers = k, iter.max = 50, nstart = 100)
    aRI <- NA
  }
  return(list('model_fit' = mod_fit, 'aRI' = aRI)) 
}