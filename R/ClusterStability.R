
ClusterStability <- function(xx, k, B = 10, frac_subsample = 0.75){
  #' Compute cluster stability of MCAP output
  #'
  #' Compute cluster stability w.r.t. pairs of data points over B random 
  #' subsamples of size m < n (where n is the total sample size). Default for 
  #' m=0.75n. Clustering is done using a full covariance Gaussian mixture model.
  #' Measure of cluster stability is based on the adjusted Rand index
  #' of shared points given two subsets of the data. 
  #' 
  #' @param xx The data matrix (n x p).
  #' @param k The number of clusters. 
  #' @param B The number of subsamples to be used (default: B=10). Note that if
  #'          B is too small (<10), the variance of the stability estimate increases.
  #' @param frac_subsample Fraction of total samples to be used in each subsample (default: 0.75).
  #' 
  #' @return @param stab Measure of cluster stability.
  #' @export
  
  ## preliminaries
  n_tot <- nrow(xx)                         #total number of samples
  n_sub <- floor(frac_subsample * n_tot)    #number of samples per subsample
  
  I <- matrix(0, B, n_sub)                  #rows: sample indices for each subsample
  A <- I                                    #rows: assignments for each subsample
  
  ## make matrices of indices and assignments
  for(b in seq(B)){ 
    idx <- sample(n_tot, n_sub)
    I[b,] <- idx
    xx_sub <- xx[idx,]
    
    ## perform clustering (full covariance Gaussian mixtures)
    clust_model <- MyGMM3(xx_sub, k)
    A[b,] <- clust_model$mod.fit$comp
    
    ## [LEGACY] emergency solution if none of the initialisations for MyGMM3 work
    #A[b,] <- MyMCLUST(xx_sub, k)$mod.fit$classification
    #A[b,] <- kmeans(xx_sub, k)$cluster
  }
  
  TT <- t(combn(1:B, 2))   #combinations of pairs of subsamples, i.e. (B choose 2)
  M <- nrow(TT)            #total number of pairs, i.e. B!/(2*(B-2)!)
  stab <- numeric(M) - 1   #vector of stability measures for each pair of assignments
  
  ## compute stability of assignments of points in the intersection
  for(m in seq(M)){
    i <- TT[m, 1] 
    j <- TT[m, 2]
    idx1 <- I[i,]
    idx2 <- I[j,]
    assignments1 <- A[i,] 
    assignments2 <- A[j,]
    
    shared_idx <- intersect(idx1, idx2)
    shared_idx1 <- which(idx1 %in% shared_idx)
    shared_idx2 <- which(idx2 %in% shared_idx)
    
    ## compute the Rand index of the two assignments on the shared points
    stab[m] <- mclust::adjustedRandIndex(assignments1[shared_idx1], 
                                         assignments2[shared_idx2])
  }
  return(stab)
}