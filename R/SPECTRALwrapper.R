
SPECTRALwrapper <- function(xx, k, true_labels = NULL, verbose = FALSE){
  #' Wrapper to do spectral clustering
  #'
  #' Spectral clustering using an RBF kernel. 
  #' The wrapper performs multiple restarts in case true labels are provided 
  #' (in order to optimise cluster assignments w.r.t. the Rand index).
  #' 
  #' @author Bernd Taschler \email{bernd.taschler@dzne.de}
  #' @seealso \code{\link{MCAPfit}}, 
  #'          \code{\link{GMMwrapper}}, 
  #'          \code{\link[kernlab]{specc}}, 
  #'          \code{\link[mclust]{adjustedRandIndex}}
  #' 
  #' @param xx The data matrix (n x p).
  #' @param k The number of clusters.
  #' @param true_labels Vector of true cluster assignments (when provided, it is 
  #'                    used to compute the Rand index). 
  #' @param verbose Logical, when true: print progress information. 
  #' 
  #' @return \item{model_fit}{ Model fit (output of \code{\link{mixglasso}})}.
  #'         \item{aRI}{ Adjusted Rand index (when \code{true_labels} is provided)}.
  #' @export
  
  ## input checks
  if(!is.null(true_labels)){ 
    stopifnot(length(true_labels) == nrow(xx)) 
    stopifnot(length(unique(true_labels)) == k)
  }

  ## preliminaries
  p_max <- 2e4
  if(ncol(xx) > p_max){
    warning(' Spectral clustering: data matrix too large - aborting ...')
    return(list('model_fit' = NA, 'aRI' = NA))
  }  
  
  if(length(true_labels)>0){
    n_repeats <- 5
    aRI_arr <- numeric(n_repeats)
    mod_fit_arr <- vector('list', n_repeats)
    
    init_max <- 100
    for(i in seq(n_repeats)){
      init_count <- 0
      while (init_count < init_max) {
        mod_fit_arr[[i]] <- tryCatch(kernlab::specc(xx, centers = k,  
                                                    kernel = 'rbfdot', 
                                                    iterations = 1000),
                                     error = identity)
        if(!methods::is(mod_fit_arr[[i]], 'error')){ break }
        init_count <- init_count + 1
      }
      if(init_count >= init_max){ aRI_arr[i] <- NA
      }else{ aRI_arr[i] <- mclust::adjustedRandIndex(true_labels, 
                                                     mod_fit_arr[[i]]@.Data) }
    }
    max_idx <- which(aRI_arr == max(aRI_arr, na.rm = TRUE))[1]
    if(length(max_idx) < 1){
      warning('! Spectral clustering: all initialisations returned an error !')
      mod_fit <- NA
      aRI <- NA
    }else{
      mod_fit <- mod_fit_arr[[max_idx]]
      aRI <- mean(aRI_arr, na.rm = TRUE)  
    }
    
  }else{  #no true labels provided
    init_count <- 0
    while (init_count < init_max) {
      mod_fit <- tryCatch(kernlab::specc(xx, centers = k,  kernel = 'rbfdot', 
                                         iterations = 1000),
                          error = identity)
      if (!methods::is(mod_fit, 'error')) {break}
      mod_fit <- NA
      init_count <- init_count + 1
    }
    if(is.na(mod_fit)){ 
      warning('! Spectral clustering: all initialisations returned an error !') 
    }
    aRI <- NA
  }
  
  return(list('model_fit' = mod_fit, 'aRI' = aRI)) 
}