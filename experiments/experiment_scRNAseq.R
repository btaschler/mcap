
# run_scRNAseq <- function(xx, k, projection = 'PCA',
#                                 true_labels = NULL, centering_per_group = FALSE, 
#                                 parallel = FALSE, verbose = FALSE, ...){
#   
  # RunMethodsParallel(fid.data, methods.arr, isTCGA = FALSE, doCentre = FALSE,
  #                    parallel=TRUE, verbose=FALSE)
  # Purpose: run clustering methods on input data sets in parallel
  # Note: parallelisation is over input data sets. Running this switches the 
  #       optimisation of projection dimension (GGMM_a, RP_a) to serial. Thus,
  #       RunMethodsParallel() is only useful if the number of data sets is >10,
  #       otherwise, RunMethods() should be faster. 
  #
  # Args:
  #   fid.data: data files (.rds) containing a list with data sets 
  #             for a given scenario
  #   methods.arr: array of clustering method names to be used 
  #   isTCGA: logical indicator whether input data is TCGA or not
  #           [default: isTCGA = FALSE]
  #   doCentre: indicator whether to centre data befor clustering 
  #             [default = FALSE]
  #   parallel: logical, if true run in parallel
  #             [default: parallel = TRUE]
  #   verbose: logical, if true print process/computation information
  #            [default: verbose = FALSE]
  #
  # OUTPUT:
  #   results: adjusted Rand index for each method
  
  

library(mclust, quietly=TRUE)
library(nethet)
library(MCMCpack)
library(mvtnorm)
library(matrixcalc)
library(pcaMethods)
library(kernlab)
library(flexclust)
library(RandPro)
library(ggplot2)
library(doParallel)
library(huge)
library(tidyverse)
library(mcap)


source('./experiments/runExperiment_scRNAseq.R')




### preliminaries
rm_group_means <- FALSE   #remove group means before clustering 
run_parallel <- FALSE      #run in parallel over multiple input files 
verbose <- TRUE          #print progress information
dir_out <- file.path(getwd(), 'experiments', 'results')    #output directory
if(!dir.exists(dir_out)){ dir.create(dir_out) }


## data
fid_path <- '//fileserver.dzne.de/taschlerb/DZNE/data/MCAP/generated_data/'
fid_list <- list.files(fid_path)[1:80]

## select files to run
fid_run <- file.path(fid_path, fid_list[c(1:5, 11:15, 21:25, 31:35)+5])



## select methods to use
methods_all <- c('mclust', 'mclust_k', 'mclust_r10', 'mclust_a',           
                 'GGMM', 'GGMM_k', 'GGMM_r10', 'GGMM_a', 'GGMM_o',
                 'RP_k', 'RP_rn', 'RP_r10', 'RP_a', 'RP_o',                 
                 'RP_sparse_a', 'RP_verysparse_a',
                 'KM', 'KMPP', 'hclust', 'specc', 'mixGLasso')                      

methods_run <- c('GGMM_a', 'RP_a', 'RP_sparse_a', 'RP_verysparse_a')




### run experiment
results <- runExperiment_scRNAseq(fid_data = fid_run,
                               methods_arr = methods_run,
                               dir_out = dir_out,
                               doCentre = rm_group_means, 
                               parallel = run_parallel, verbose = verbose)     
closeAllConnections()




# }

