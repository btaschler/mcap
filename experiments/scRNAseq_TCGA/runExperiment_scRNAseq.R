
runExperiment_scRNAseq <- function(fid_data, methods_arr, dir_out = NULL, 
                                   doCentre = FALSE, parallel = TRUE, 
                                   verbose = FALSE){
  # RunMethodsParallel(fid_data, methods_arr, doCentre = FALSE,
  #                    parallel=TRUE, verbose=FALSE)
  # Purpose: run clustering methods on input data sets in parallel
  # Note: parallelisation is over input data sets. Running this switches the 
  #       optimisation of projection dimension (GGMM_a, RP_a) to serial. Thus,
  #       RunMethodsParallel() is only useful if the number of data sets is >10,
  #       otherwise, RunMethods() should be faster. 
  #
  # Args:
  #   fid_data: data files (.rds) containing a list with data sets 
  #             for a given scenario
  #   methods_arr: array of clustering method names to be used 
  #   dir_out: output directory
  #   doCentre: indicator whether to centre data befor clustering 
  #             [default = FALSE]
  #   parallel: logical, if true run in parallel
  #             [default: parallel = TRUE]
  #   verbose: logical, if true print process/computation information
  #            [default: verbose = FALSE]
  #
  # OUTPUT:
  #   results: adjusted Rand index for each method
  
  
  ## preliminaries ------------------------------------------------------------
  
  # ## perform clustering [in parallel]
  # if(parallel){
  #   no_cores <- min(length(fid_data), detectCores()-2)
  #   cl <- makeCluster(no_cores)
  #   registerDoParallel(cores=no_cores)
  #   cat('\n run in parallel on ', getDoParWorkers(), 'cores ... \n')
  #   verbose <- FALSE
  # }else{
  #   registerDoSEQ()
  #   cat('\n run sequentially ... \n')
  # }
  # par_results <- foreach(ds_par=1:length(fid_data), 
  #                        .export=c('OptDimClusterStability',
  #                                  'ClusterStability',
  #                                  'MCAPfit',
  #                                  'GMMwrapper',
  #                                  'GramPCA',
  #                                  'RandProject',
  #                                  'KMwrapper',
  #                                  'HCLUSTwrapper',
  #                                  'SPECTRALwrapper',
  #                                  'CentrePerGroup'
  #                        ),
  #                        .packages=c('mclust',
  #                                    'mvtnorm',
  #                                    'matrixcalc',
  #                                    'pcaMethods',
  #                                    'kernlab',
  #                                    'flexclust',
  #                                    'nethet',
  #                                    'RandPro')) %dopar% { 
     
                                    
   for(ds_par in seq_along(fid_data)){   #HACK because parallel does not work !!
     
   ### clustering - main body -------------------------------------------------
   ds <- fid_data[ds_par]
   results <- vector('list', 6)
   names(results) <- c('dataID', 'method', 'projDim', 'gridpoint', 'gridvalue', 'aRI')
   
   if(parallel){ setMKLthreads(1) }  #important to prevent multi-threading (!)
   
   if(verbose){ cat('\n load data set:', ds, ' ...\n') } 
   data <- readRDS(ds)
   n_grid <- length(data$xx)
   
   ## loop over grid points in the data set
   for(g in seq(n_grid)){ 
     if(verbose){ cat(' gridpoint', g, '/',n_grid, '...') } 
     curr_xx <- data$xx[[g]]
     curr_labels <- data$labels[[g]]
     curr_k <- length(unique(curr_labels))
     curr_n <- nrow(curr_xx)
     curr_p <- ncol(curr_xx)
     
     ## mean centre the data per group
     if(doCentre){ curr_xx <- CentrePerGroup(curr_xx, curr_labels) }
     
     if(any(c('GGMM_a', 'GGMM_o', 'mclust_a') %in% methods_arr)){
       if(verbose){ cat(' opt. dim. (PCA)') } 
       dummy1 <- OptDimClusterStability(xx = curr_xx, k = curr_k, 
                                        method='PCA', n_grid = 5,
                                        true_labels = curr_labels,
                                        parallel = TRUE)         #NOTE: switch off if outer parallel loop enabled!
       q_opt_pca <- dummy1$q_opt
       q_star_pca <- dummy1$q_oracle
     }
     if(any(c('RP_a', 'RP_o') %in% methods_arr)){
       if(verbose){ cat(' opt. dim. (RP-Gaussian)') } 
       dummy2 <- OptDimClusterStability(xx = curr_xx, k = curr_k, 
                                        method='gaussian', n_grid = 5,
                                        true_labels = curr_labels,
                                        parallel = TRUE)         #NOTE: switch off if outer parallel loop enabled!
       q_opt_rp <- dummy2$q_opt
       q_star_rp <- dummy2$q_oracle
     }
     
     ## loop over clustering methods
     for(m in methods_arr){ 
       curr_q <- ncol(curr_xx)   #default: no projection

       if(m=='GGMM'){
         if(verbose){ cat('| GGMM') } 
         curr_q <- curr_k
         curr_fit <- GMMwrapper(curr_xx, 
                            k = curr_k, true_labels = curr_labels)     
       }else if(m=='GGMM_k'){
         if(verbose){ cat('| MCAP-K') } 
         curr_q <- curr_k
         curr_fit <- GMMwrapper(GramPCA(curr_xx, npc = curr_q)$zz, 
                            k = curr_k, true_labels = curr_labels)  
       }else if(m=='GGMM_r10'){
         if(verbose){ cat('| MCAP-r10') } 
         curr_q <- max(curr_k, min(curr_p, floor(sqrt(10 * curr_n / curr_k))))
         curr_fit <- GMMwrapper(GramPCA(curr_xx, npc = curr_q)$zz, 
                            k = curr_k, true_labels = curr_labels)
       }else if(m=='GGMM_a'){
         if(verbose){ cat('| MCAP-adapt.') } 
         curr_q <- q_opt_pca
         curr_fit <- GMMwrapper(GramPCA(curr_xx, npc = curr_q)$zz, 
                            k = curr_k, true_labels = curr_labels)
       }else if(m=='GGMM_o'){
         if(verbose){ cat('| MCAP-oracle') } 
         curr_q <- q_star_pca
         curr_fit <- GMMwrapper(GramPCA(curr_xx, npc = curr_q)$zz, 
                            k = curr_k, true_labels = curr_labels)

       }else if(m=='RP_k'){
         if(verbose){ cat('| MCAP-RP-Gaussian-K') } 
         curr_q <- curr_k
         curr_fit <- GMMwrapper(RandProject(curr_xx, q = curr_q, method = 'gaussian'), 
                            k = curr_k, true_labels = curr_labels)
       }else if(m=='RP_a'){
         if(verbose){ cat('| MCAP-RP-Gaussian-adapt.') } 
         curr_q <- q_opt_rp
         curr_fit <- GMMwrapper(RandProject(curr_xx, q = curr_q, method = 'gaussian'), 
                            k = curr_k, true_labels = curr_labels)
       }else if(m=='RP_o'){
         if(verbose){ cat('| MCAP-RP-Gaussian-oracle') } 
         curr_q <- q_star_rp
         curr_fit <- GMMwrapper(RandProject(curr_xx, q = curr_q, method = 'gaussian'), 
                            k = curr_k, true_labels = curr_labels)
       }else if(m=='RP_sparse_a'){
         if(verbose){ cat('opt. dim. (RP-Achlioptas)') } 
         dummy <- OptDimClusterStability(curr_xx, k = curr_k, 
                                         true_labels = curr_labels, 
                                         method = 'achlioptas',
                                         parallel = TRUE)  #NOTE: switch off if outer parallel loop enabled!
         if(verbose){ cat('| MCAP-RP-Achlioptas-adapt.') } 
         curr_q <- dummy$q_opt
         curr_fit <- GMMwrapper(RandProject(curr_xx, q = curr_q, method = 'achlioptas'), 
                            k = curr_k, true_labels = curr_labels)
       }else if(m=='RP_verysparse_a'){
         if(verbose){ cat('opt. dim. (RP-Li)') } 
         dummy <- OptDimClusterStability(curr_xx, k = curr_k, 
                                         true_labels = curr_labels, 
                                         method = 'li',
                                         parallel = TRUE)  #NOTE: switch off if outer parallel loop enabled!
         if(verbose){ cat('| MCAP-RP-Li-adapt.') } 
         curr_q <- dummy$q_opt
         curr_fit <- GMMwrapper(RandProject(curr_xx, q = curr_q, method = 'li'), 
                            k = curr_k, true_labels = curr_labels)
         
       }else if(m=='KM'){
         if(verbose){ cat('| KM') } 
         curr_fit <- KMwrapper(curr_xx, k = curr_k, true_labels = curr_labels)
       }else if(m=='KMPP'){
         if(verbose){ cat('| KM++') } 
         curr_fit <- KMwrapperPP(curr_xx, k = curr_k, true_labels = curr_labels)
       }else if(m=='hclust'){
         if(verbose){ cat('| hclust') } 
         curr_fit <- HCLUSTwrapper(curr_xx, k = curr_k, true_labels = curr_labels)
       }else if(m=='specc'){
         if(verbose){ cat('| spectral') } 
         curr_fit <- SPECTRALwrapper(curr_xx, k = curr_k, true_labels = curr_labels)

       }else{
         message(sprintf('\n ! unknown method: %s !', m))
         error(' ...aborting...')
       }
       
       ## collect results
       results$dataID <- substr(ds, nchar(ds)-10, nchar(ds)-4)
       results$method <- m
       results$projDim <- curr_q
       results$gridpoint <- g
       results$gridvalue <- data$gridvalues[g]
       results$aRI <- unlist(curr_fit['aRI'])
       
       ## save results to file
       if(is.null(dir_out)){ dir_out <- getwd() }
       fid_out <- file.path(dir_out, 'scRNAseq_')
       if(doCentre){ fid_out <- paste0(fid_out, 'd0_') }
       fid_out <- paste0(fid_out, results$dataID[length(results$dataID)], '.rds')
       
       if(file.exists(fid_out)){ previous_res <- readRDS(fid_out) 
       }else{ previous_res <- NULL }
       curr_res <- rbind(previous_res, results)
       saveRDS(curr_res, file = fid_out)
     }
     if(verbose){ cat(' ... done! \n\n') }
   }
  }
  
  
  ## shut down parallel computing
  if(parallel){ stopCluster(cl) }
  closeAllConnections()
  return(NA)
}


