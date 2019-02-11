
OptDimClusterStability_multimodal <- function(xa, xb, k, method = 'PCA', 
                                   n_grid = 5, 
                                   qa_max = min(ncol(xa), sqrt(10*nrow(xa)/k)),
                                   qb_max = min(ncol(xb), sqrt(10*nrow(xb)/k)),
                                   true_labels = NULL, 
                                   parallel = FALSE, verbose = FALSE){

  
  ## preliminaries
  n <- nrow(xa)
  pa <- ncol(xa)
  pb <- ncol(xb)
  p <- pa + pb
  
  ## input checks
  stopifnot(qa_max <= p)
  stopifnot(method %in% c('PCA', 'gaussian', 'achlioptas', 'li'))
  
  ## hyperparameter for line search over target dimension
  ca_arr <- c(k^3 / n, seq(sqrt(0.5), sqrt(max(k, qa_max^2*k/n)), length = n_grid-1)^2)
  cb_arr <- c(k^3 / n, seq(sqrt(0.5), sqrt(max(k, qa_max^2*k/n)), length = n_grid-1)^2)
  
  ## set up parallel or serial computation
  if(parallel){
    num_cores <- min(length(ca_arr)+length(cb_arr), parallel::detectCores()-2)
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
  ij <- 0  #initialise to prevent warning in code check that i is a global variable
  results_tbl <- foreach::foreach(ij=1:(length(ca_arr)*length(cb_arr)), 
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
    
    i <- ceiling(ij / length(cb_arr))
    j <- (ij %% length(cb_arr)) + 1
    
    curr_qa <- max(k, min(pa, floor(sqrt(ca_arr[i] * n / k))))
    curr_qb <- max(k, min(pb, floor(sqrt(cb_arr[j] * n / k))))
    
    
    ## compute cluster stability for given target dimension
    if(method == 'PCA'){
      curr_stability <- mean(ClusterStability(xx = cbind(GramPCA(xa, curr_qa)$zz, 
                                                         GramPCA(xb, curr_qb)$zz),
                                              k = k, B = 10, 
                                              frac_subsample = 0.75), na.rm=TRUE)
      if(!is.null(true_labels)){
        curr_aRI <- GMMwrapper(xx = cbind(GramPCA(xa, curr_qa)$zz, 
                                          GramPCA(xb, curr_qb)$zz),
                               k = k, true_labels = true_labels)$aRI
      }else{ curr_aRI <- NA }
      
    }else{ #random projections
      curr_stability <- mean(ClusterStability(xx = cbind(RandProject(xa, q = curr_qa, 
                                                               method = method),
                                                         RandProject(xb, q = curr_qb, 
                                                                     method = method)),
                                              k = k, B = 10, 
                                              frac_subsample = 0.75), na.rm=TRUE)
      if(!is.null(true_labels)){
        curr_aRI <- GMMwrapper(xx = cbind(RandProject(xa, q = curr_qa, 
                                                      method = method),
                                          RandProject(xb, q = curr_qb, 
                                                      method = method)),
                               k = k, true_labels = true_labels)$aRI
      }else{ curr_aRI <- NA }
    }
    
    ## combine results
    tbl_out <- dplyr::tibble('qa' = curr_qa, 
                             'qb' = curr_qb,
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
  qa_opt <- tmp$qa
  qb_opt <- tmp$qb
  stab_max <- tmp$stability
  
  ## find oracle target dimension
  if(!is.null(true_labels)){
    if(sum(is.na(results_tbl$aRI)) == nrow(results_tbl)){
      q_oracle <- NA 
    }else{
      tmp <- (results_tbl %>% dplyr::filter(results_tbl$aRI,
                                            results_tbl$aRI == max(results_tbl$aRI, na.rm = TRUE)))[1,]
      qa_oracle <- tmp$qa
      qb_oracle <- tmp$qb
      aRI_max <- tmp$aRI
    }
  }else{ q_oracle <- NA }
  
  return(list('qa_opt' = qa_opt, 'qb_opt' = qb_opt, 'stab_score' = stab_max, 
              'qa_oracle' = qa_oracle, 'qb_oracle' = qb_oracle, 
              'aRI_oracle' = aRI_max, 'full_tbl' = results_tbl)) 
  }