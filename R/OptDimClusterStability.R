
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
  #' @author Bernd Taschler \email{bernd.taschler@dzne.de}
  #' @author Sach Mukherjee \email{sach.mukherjee@dzne.de}
  #' @seealso \code{\link{GMMwrapper}}
  #' @seealso \code{\link{MCAPfit}}
  #' 
  #' @param xx The data matrix (n x p).
  #' @param k The number of clusters.
  #' @param method Projection method (\code{"PCA"} or random projections: 
  #'               \code{"gaussian"}, \code{"achlioptas"} or \code{"li"}). 
  #'               Default: \code{"PCA"}.
  #' @param n_grid Number of values to be used in the line search for optimal
  #'               projection dimension. Default: 5.
  #' @param q_max Maximum target dimension to be used in line search. (Note: 
  #'              the smallest target dimension is always \code{k}, 
  #'              the maximum may not exceed the total dimensionality p). 
  #'              Default: \code{sqrt(10n/k)}.
  #' @param true_labels Vector of true cluster assignments (if provided, it is 
  #'                    used to compute the Rand index and \code{q_star}).
  #' @param parallel Logical, if true: perform line search over \code{q} in parallel. 
  #' @param verbose Logical, if true: print progress information. 
  #' 
  #' @return \item{q_opt}{ Optimal target dimension (maximises cluster stability).}
  #' @return \item{stab_score}{ Stability measure for \code{q_opt}.}
  #' @return \item{q_star}{ Optimal ("oracle") target dimension (maximises adj. 
  #'                        Rand index). Only available if true labels have been 
  #'                        provided.}
  #' @importFrom magrittr %>%
  #' @importFrom foreach %dopar%
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
    num_cores <- min(length(c_arr), parallel::detectCores()-2)
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cores = num_cores)
    if(verbose){ 
      cat('\n run line search over target dimension in parallel on ', 
          foreach::getDoParWorkers(), 'cores ... \n')
    }
  }else{
    foreach::registerDoSEQ()
    if(verbose){ 
      cat('\n run line search over target dimension sequentially ... \n')
    }
  }
  
  ## perform line search over target dimension q
  i <- 0  #initialise to prevent warning in code check that i is a global variable
  results_tbl <- foreach::foreach(i=1:length(c_arr), 
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
                                  .init = dplyr::tibble('q' = integer(), 
                                                            'aRI' = numeric(),
                                                            'stability' = numeric())) %dopar% 
  { 
    if(parallel){ RevoUtilsMath::setMKLthreads(1) }  #prevent multi-threading (!)
    
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
    tbl_out <- dplyr::tibble('q' = curr_q, 
                                 'aRI' = curr_aRI,
                                 'stability' = curr_stability)
    return(tbl_out)
  }
  ## shut down parallel computing
  if(parallel){ parallel::stopCluster(cl) }
  closeAllConnections()
  
    
  ## find optimal target dimension
  tmp <- (results_tbl %>% dplyr::filter(results_tbl$stability, 
            results_tbl$stability == max(results_tbl$stability, na.rm = TRUE)))[1,]
  q_opt <- tmp$q
  stab_max <- tmp$stability
  
  ## find oracle target dimension
  if(!is.null(true_labels)){
    if(sum(is.na(results_tbl$aRI)) == nrow(results_tbl)){
      q_oracle <- NA 
    }else{
      q_oracle <- (results_tbl %>% dplyr::filter(results_tbl$aRI, 
                     results_tbl$aRI == max(results_tbl$aRI, na.rm = TRUE)))$q[1]
    }
  }else{ q_oracle <- NA }

  return(list('q_opt' = q_opt, 'stab_score' = stab_max, 'q_oracle' = q_oracle)) 
}