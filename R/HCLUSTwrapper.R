
HCLUSTwrapper <- function(xx, k, method = 'ward.D', 
                          true_labels = NULL, verbose = FALSE){
  #' Wrapper to do hierarchical clustering
  #'
  #' Hierarchical clustering using Euclidean distances and (by default) 
  #' the Ward criterion. 
  #' The wrapper can perform an optimisation over clustering methods in case 
  #' true labels are provided (i.e. optimising cluster assignments w.r.t. the Rand index).
  #' 
  #' @author Bernd Taschler: \email{bernd.taschler@dzne.de}
  #' @seealso \code{\link[stats]{hclust}}, 
  #'          \code{\link[mclust]{adjustedRandIndex}}
  #' 
  #' @param xx The data matrix (n x p).
  #' @param k The number of clusters.
  #' @param method Clustering method (see \code{\link[stats]{hclust}}). 
  #'               If \code{method = "all"} and \code{true_labels} is provided: 
  #'               optimise clustering w.r.t. Rand index  and return method 
  #'               with highest adjusted RI. Default: \code{method = "Ward.D"}.
  #' @param true_labels Vector of true cluster assignments (when provided, it is 
  #'                    used to compute the Rand index). 
  #' @param verbose Logical, when true: print progress information. 
  #' 
  #' @return \item{model_fit}{ Model fit (output of \code{\link[stats]{hclust}}).}
  #'         \item{aRI}{ Adjusted Rand index (when \code{true_labels} is provided).}
  #' @examples 
  #'   HCLUSTwrapper(xx = matrix(rnorm(500),50,10), k = 2)
  #' @export
  
  ## input checks
  if(!is.null(true_labels)){ 
    stopifnot(length(true_labels) == nrow(xx)) 
    stopifnot(length(unique(true_labels)) == k)
    stopifnot(method %in% c('ward.D', 'single', 'complete', 'average', 
                            'mcquitty','median', 'centroid', 'all'))
  }else{
    stopifnot(method %in% c('ward.D', 'single', 'complete', 'average', 
                            'mcquitty','median', 'centroid'))
  }
  
  ## perform hierarchical clustering
  if(length(true_labels) > 0){ 
    if(method == 'all'){
      method_arr <- c('ward.D', 'single', 'complete', 'average',
                      'mcquitty', 'median', 'centroid')
    }else{ method_arr <- method }
    
    n_repeats <- length(method_arr)
    aRI_arr <- numeric(n_repeats)
    mod_fit_arr <- vector('list', n_repeats)    
    for(i in seq(n_repeats)){
      mod_fit_arr[[i]] <- stats::hclust(stats::dist(xx, method = 'euclidean'), 
                                        method = method_arr[i])
      aRI_arr[i] <- mclust::adjustedRandIndex(true_labels, 
                                              stats::cutree(mod_fit_arr[[i]], 
                                                            k=k)) 
    }
    max_idx <- which(aRI_arr == max(aRI_arr, na.rm = TRUE))[1]
    mod_fit <- mod_fit_arr[[max_idx]]
    aRI <- aRI_arr[max_idx]  
    
  }else{  #no true labels provided
    mod_fit <- stats::hclust(stats::dist(xx, method = 'euclidean'), 
                             method = method)
    aRI <- NA
  }
  
  return(list('model_fit' = mod_fit, 'aRI' = aRI)) 
}