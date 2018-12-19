
OptDimClusterStability <- function(xx, k, method, true_labels = NULL, 
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
  #' @param method Projection method (PCA or random projections).
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
  c.arr <- c(k^3 / n, seq(0.5, max(k, 10), length = 5)) 
  #c.arr <- c(seq(2,30,2),35,40)                #for fine grid search only
  stab.max <- -Inf
  aRI.max <- -Inf
  oracle_dim <- 0
  curr.aRI <- 0
  
  ## run stability computation in parallel over c.arr
  if(parallel){ 
    n_cores <- min(detectCores()-4, length(c.arr))
    cl <- makeCluster(n_cores)
    registerDoParallel(cores=n_cores)
    par_results <- foreach(i=1:length(c.arr), 
                           .export=c('MyGMM','MyGMM3','GramPCA','RandPro',
                                     'ComputeClusterStability','MyKM','MyMCLUST'),
                           .packages=c('mclust','nethet','pcaMethods','RandPro'),
                           .combine = rbind, .multicombine=TRUE, 
                           .init = vector(mode='list', 3)) %dopar% 
      { 
        c <- c.arr[i]
        curr_dim <- max(k, min(p, floor(sqrt(c * n / k))))
        #curr_dim <- c.arr[i]                           #for fine grid search only
        if(method == 'PCA'){
         curr.stability <- mean(ComputeClusterStability(GramPCA(xx, curr_dim)$zz, 
                                                        k, true.labels), na.rm=TRUE)
         curr.aRI <- MyGMM3(GramPCA(xx, npc = curr_dim)$zz, 
                            k = k, true.labels = true.labels)$aRI
        }else if(method == 'RP'){
         curr.stability <- mean(ComputeClusterStability(RandPro(xx, q = curr_dim, 
                                                                method = 'Gaussian'), 
                                                        k, true.labels), na.rm=TRUE)
         curr.aRI <- MyGMM3(RandPro(xx, q = curr_dim, method = 'Gaussian'), 
                            k = k, true.labels = true.labels)$aRI
        }else if(method %in% c('achlioptas', 'li')){
         curr.stability <- mean(ComputeClusterStability(RandPro(xx, q = curr_dim, 
                                                                method = method), 
                                                        k, true.labels), na.rm=TRUE)
         curr.aRI <- MyGMM3(RandPro(xx, q = curr_dim, method = method), 
                            k = k, true.labels = true.labels)$aRI
        }else{message(' ! unknown method !'); return(NA)}
        
        return(list(curr.stability, curr_dim, curr.aRI))
      }
    stopCluster(cl)
    closeAllConnections()
    
    ## combine results
    stab_arr <- unlist(par_results[,1])
    dim_arr <- unlist(par_results[,2])
    ari_arr <- unlist(par_results[,3])
    
    #cat(stab_arr)   #diagnosis
    
    ## find optima in combined results
    stab_max_idx <- which(stab_arr == max(stab_arr, na.rm=TRUE))[1]
    opt_dim <- dim_arr[stab_max_idx]
    
    ari_max_idx <- which(ari_arr == max(ari_arr, na.rm=TRUE))[1]
    if(length(ari_max_idx)>0){ 
      oracle_dim <- dim_arr[ari_max_idx] 
    }else{ 
      oracle_dim <- NA 
    }
    
    ## run stability computation in series  
  }else{ 
    for(c in c.arr){
      curr_dim <- max(k, min(p, floor(sqrt(c * n / k))))
      if(verbose){ message(sprintf(' compute cluster stability for q=%i', curr_dim)) }
      if(method == 'PCA'){
        curr.stability <- mean(ComputeClusterStability(GramPCA(xx, curr_dim)$zz, 
                                                       k, true.labels), na.rm=TRUE)
        curr.aRI <- MyGMM3(GramPCA(xx, npc = curr_dim)$zz, 
                           k = k, true.labels = true.labels)$aRI
      }else if(method == 'RP'){
        curr.stability <- mean(ComputeClusterStability(RandPro(xx, q = curr_dim, 
                                                               method = 'Gaussian'), 
                                                       k, true.labels), na.rm=TRUE)
        curr.aRI <- MyGMM3(RandPro(xx, q = curr_dim, method = 'Gaussian'), 
                           k = k, true.labels = true.labels)$aRI
      }else if(method %in% c('achlioptas', 'li')){
        curr.stability <- mean(ComputeClusterStability(RandPro(xx, q = curr_dim, 
                                                               method = method), 
                                                       k, true.labels), na.rm=TRUE)
        curr.aRI <- MyGMM3(RandPro(xx, q = curr_dim, method = method), 
                           k = k, true.labels = true.labels)$aRI
      }else{message(' ! unknown method !'); return(NA)}
      
      #print(curr.stability)
      if(curr.stability > stab.max){
        stab.max <- curr.stability
        opt_dim <- curr_dim
      }
      #print(curr.aRI)
      if(length(curr.aRI)>0 && !is.na(curr.aRI) && curr.aRI > aRI.max){
        aRI.max <- curr.aRI
        oracle_dim <- curr_dim
      }
    }
  }
  return(list('opt_q'=opt_dim, 'stab_score' = stab.max, 'q_star'=oracle_dim)) 
}