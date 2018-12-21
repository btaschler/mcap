
GMMwrapper <- function(xx, k, criterion = 'mmdl', 
                       true_labels = NULL, verbose = FALSE){
  #' Wrapper to do clustering with Gaussian mixture models
  #'
  #' Full covariance Gaussian mixture modelling based on the nethet package. 
  #' The wrapper performs multiple restarts in case true labels are provided 
  #' (in order to optimise cluster assignments w.r.t. the Rand index).
  #' 
  #' @author Bernd Taschler: \email{bernd.taschler@dzne.de}
  #' @author Sach Mukherjee: \email{sach.mukherjee@dzne.de}
  #' @seealso \code{\link{runMCAP}}
  #' @seealso \code{\link{OptDimClusterStability}}
  #' 
  #' @param xx The data matrix (n x p).
  #' @param k The number of clusters.
  #' @param criterion Optimisation criterion (\code{bic} or \code{mmdl}). 
  #'                  Default: \code{mmdl}.
  #' @param true_labels Vector of true cluster assignments (when provided, it is 
  #'                    used to compute the Rand index). 
  #' @param verbose Logical, when true: print progress information. 
  #' 
  #' @return \item{model_fit}{ Model fit (output of \code{\link[nethet]{mixglasso}}).}
  #' @return \item{aRI}{ Adjusted Rand index (when \code{true_labels} is provided).}
  #' @export
  
  ## input checks
  stopifnot(criterion %in% c('bic', 'mmdl'))
  if(!is.null(true_labels)){ 
    stopifnot(length(true_labels) == nrow(xx)) 
    stopifnot(length(unique(true_labels)) == k)
  }
  
  ## preliminaries
  max_p <- 1000
  if(ncol(xx) > max_p){
    warning(' GMM-mixglasso: data matrix too large - aborting ...')
    return(list('model_fit' = NA, 'bic' = NA, 'mmdl' = NA, 'aRI' = NA))
  }
  n_repeats <- 50
  mod_fit_arr <- vector('list',n_repeats)
  mod_results <- dplyr::tibble('bic' = numeric(n_repeats),
                               'mmdl' = numeric(n_repeats),
                               'aRI' = numeric(n_repeats))
  
  ## compute mixture models based on different initialisation variants
  init_max <- 10
  for(i in seq(n_repeats)){
    if(verbose & i%%10 == 0){ cat('=') }
    
    init_count <- 0
    if(i==1){
      while(init_count < init_max){
        curr_mod <- tryCatch(nethet::mixglasso(xx, n.comp = k,
                                    lambda = 0, 
                                    init = 'kmeans', nstart.kmeans = 50,
                                    iter.max.kmeans = 10, term = 1e-3,
                                    modelname.hc = 'VVV'), 
                          error=identity)
        if(!methods::is(curr_mod, 'error')){ break }
        init_count <- init_count+1
      }
      
    }else if(i==2){
      while(init_count < init_max){
        curr_mod <- tryCatch(nethet::mixglasso(xx, n.comp = k,
                                    lambda = 0, 
                                    init = 'hc', term = 1e-3,
                                    modelname.hc = 'VVV'), 
                          error=identity)
        if(!methods::is(curr_mod, 'error')){ break }
        init_count <- init_count+1
      }
      
    }else if(i==3){
      while(init_count < init_max){
        curr_mod <- tryCatch(nethet::mixglasso(xx, n.comp = k,
                                    lambda = 0, 
                                    init = 'kmeans.hc', nstart.kmeans = 50,
                                    iter.max.kmeans = 10, term = 1e-3,
                                    modelname.hc = 'VVV'), 
                          error=identity)
        if(!methods::is(curr_mod, 'error')){ break }
        init_count <- init_count+1
      }
      
    }else{
      while(init_count < init_max){
        curr_mod <- tryCatch(nethet::mixglasso(xx, n.comp = k,
                                    lambda = 0, 
                                    init = 'r.means', term = 1e-3,
                                    modelname.hc = 'VVV'), 
                          error=identity)
        if (!methods::is(curr_mod, 'error')) {break}
        init_count <- init_count+1
      }
    }
    
    if(is.null(curr_mod)){ 
      mod_fit_arr[[i]] <- NA 
    }else{ 
      mod_fit_arr[[i]] <- curr_mod 
    }
    
    ## store model output for current initialisation
    if(init_count < init_max){ #not all sub-initialisations failed to converge
      mod_results$bic[i] <- curr_mod$bic
      mod_results$mmdl[i] <- curr_mod$mmdl
      
      if(length(true_labels)>0){
        mod_results$aRI[i] <- mclust::adjustedRandIndex(curr_mod$comp, true_labels)
      }
    }
  }
  if(verbose){ cat(' done \n') }
  
  ## get optimal model fit
  if(sum(is.na(mod_results$bic))==n_repeats){ #all repeats over all initialisations fail
    if(verbose){ warning('! GMM: all initialisations returned an error !') }
    return(list('model_fit' = NA, 'bic' = NA, 'mmdl' = NA, 'aRI' = NA))
    
  }else if(criterion == 'mmdl'){  #take model fit that minimises MMDL
    min_idx <- which(mod_results$mmdl == min(mod_results$mmdl, na.rm=TRUE))[1]
    return(list('model_fit' = mod_fit_arr[[min_idx]], 
                'bic' = mod_results$bic[min_idx], 
                'mmdl' = mod_results$mmdl[min_idx],
                'aRI' = mod_results$aRI[min_idx])) 
  }else if(criterion == 'bic'){  #take model fit that minimises BIC
    min_idx <- which(mod_results$bic == min(mod_results$bic, na.rm=TRUE))[1]
    return(list('model_fit' = mod_fit_arr[[min_idx]], 
                'bic' = mod_results$bic[min_idx], 
                'mmdl' = mod_results$mmdl[min_idx],
                'aRI' = mod_results$aRI[min_idx])) 
  }
}