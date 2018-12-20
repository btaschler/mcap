
OptDimClusterStability <- function(xx, k, method = 'PCA', 
                                   n_grid = 5, 
                                   q_max = min(ncol(xx), sqrt(10*nrow(xx)/k)),
                                   true_labels = NULL, 
                                   parallel = FALSE, verbose = FALSE){
  #' Determine optimal projection dimension (PCA or random projections) based
  #' on cluster stability
  #'
  #' Find the optimal projection dimension for PCA or random projections 
  #' based on cluster stability by performing a line search over the target 
  #' dimension q.
  #' 
  #' @param xx The data matrix (n x p).
  #' @param k The number of clusters.
  #' @param method Projection method ('PCA' or random projections: 'gaussian',
  #'               'achlioptas' or 'li'). Default: `"PCA"`.
  #' @param n_grid Number of values to be used in the line search for optimal
  #'               projection dimension. Default: 5.
  #' @param q_max Maximum target dimension to be used in line search. (Note: 
  #'              the smallest target dimension is always `k`, the maximum may
  #'              not exceed the total dimensionality p). Default: sqrt(10n/k).
  #' @param true_labels Vector of true cluster assignments (if provided, it is 
  #'                    used to compute the Rand index and q_star).
  #' @param parallel Logical, if true: perform line search over q in parallel. 
  #' @param verbose Logical, if true: print progress information. 
  #' 
  #' @return @param q_opt Optimal target dimension (maximises cluster stability).
  #' @return @param stab_score Stability measure for q_opt. 
  #' @return @param q_star Optimal ("oracle") target dimension (maximises adj. 
  #'                       Rand index). Only available if true labels have been 
  #'                       provided. 
  #' @export
  
  ## preliminaries
  n <- nrow(xx)
  p <- ncol(xx)

  ## input checks
  stopifnot(q_max <= p)
  stopifnot(method %in% c('PCA', 'gaussian', 'achlioptas', 'li'))
  
  ## hyperparameter for line search over target dimension
  c_arr <- c(k^3 / n, seq(0.5, max(k, q_max^2*k/n), length = n_grid-1)) 
  
  ## set up parallel or serial computation
  if(parallel){
    num_cores <- min(length(c_arr), detectCores()-2)
    cl <- makeCluster(num_cores)
    registerDoParallel(cores = num_cores)
    if(verbose){ 
      cat('\n run line search over target dimension in parallel on ', 
          getDoParWorkers(), 'cores ... \n')
    }
  }else{
    registerDoSEQ()
    if(verbose){ 
      cat('\n run line search over target dimension sequentially ... \n')
    }
  }
  
  ## perform line search over target dimension q
  results_tbl <- foreach(i=1:length(c_arr), 
                         .export=c('GMMwrapper',
                                   'GramPCA',
                                   'RandProject',
                                   'ClusterStability'),
                         .packages=c('mclust',
                                     'nethet',
                                     'pcaMethods',
                                     'RandPro',
                                     'tidyverse'),
                         .combine = rbind, .multicombine=FALSE,
                         .init = tibble('q' = integer(), 
                                        'aRI' = numeric(),
                                        'stability' = numeric())) %dopar% 
  { 
    if(parallel){ setMKLthreads(1) }  #prevent multi-threading (!)
    
    ## set current target dimension
    curr_q <- max(k, min(p, floor(sqrt(c_arr[i] * n / k))))
    
    ## compute cluster stability for given target dimension
    if(method == 'PCA'){
      curr_stability <- mean(ClusterStability(xx = GramPCA(xx, curr_q)$zz, 
                                              k = k, B = 10, 
                                              frac_subsample = 0.75), na.rm=TRUE)
      if(!is.null(true_labels)){
        curr_aRI <- GMMwrapper(xx = GramPCA(xx, npc = curr_q)$zz, 
                               k = k, true_labels = true_labels)$aRI
      }else{ curr_aRI <- NA }
     
    }else{ #random projections
      curr_stability <- mean(ClusterStability(xx = RandProject(xx, q = curr_q, 
                                                               method = method), 
                                              k = k, B = 10, 
                                              frac_subsample = 0.75), na.rm=TRUE)
      if(!is.null(true_labels)){
        curr_aRI <- GMMwrapper(xx = RandProject(xx, q = curr_q, method = method), 
                               k = k, true_labels = true_labels)$aRI
      }else{ curr_aRI <- NA }
    }
    
    ## combine results
    tbl_out <- tibble('q' = curr_q, 
                      'aRI' = curr_aRI,
                      'stability' = curr_stability)
    return(tbl_out)
  }
  ## shut down parallel computing
  if(parallel){ stopCluster(cl) }
  closeAllConnections()
  
    
  ## find optimal target dimension
  tmp <- (results_tbl %>% filter(., stability == max(stability, na.rm = TRUE)))[1,]
  q_opt <- tmp$q
  stab_max <- tmp$stability
  
  ## find oracle target dimension
  if(!is.null(true_labels)){
    if(sum(is.na(results_tbl$aRI)) == nrow(results_tbl)){
      q_oracle <- NA 
    }else{
      q_oracle <- (results_tbl %>% filter(., aRI == max(aRI, na.rm = TRUE)))$q[1]
    }
  }else{ q_oracle <- NA }

  return(list('q_opt' = q_opt, 'stab_score' = stab_max, 'q_oracle' = q_oracle)) 
}