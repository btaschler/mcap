
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


source('./experiments/runExperiment_TCGA.R')




### preliminaries
rm_group_means <- FALSE   #remove group means before clustering 
run_parallel <- FALSE     #run in parallel over multiple input files 
verbose <- TRUE           #print progress information
dir_out <- file.path(getwd(), 'experiments', 'results','TCGA')    #output directory
if(!dir.exists(dir_out)){ dir.create(dir_out) }


## data
fid_path <- '//fileserver.dzne.de/ag-mukherjee/MCAP/TCGA_data/'

fid_list <- paste0(fid_path,
                   c(paste0('fig1_n/tcga_gex_norm_n_',seq(10),'.rds'),                    #  0
                     paste0('fig2_p/tcga_gex_norm_p_',seq(10),'.rds'),                    # 10
                     paste0('fig3_d/tcga_gex_norm_d_',seq(10),'.rds'),                    # 20
                     paste0('fig4_k/tcga_gex_norm_k_',seq(10),'.rds')                     # 30
                   ))


## select files to run
fid_run <- fid_list[31:35]
fid_run = fid_list[c(1:3,11:13,31:33)]; rm_group_means = TRUE


## select methods to use
methods_all <- c('mclust', 'mclust_k', 'mclust_r10', 'mclust_a',           
                 'GGMM', 'GGMM_k', 'GGMM_r10', 'GGMM_a', 'GGMM_o',
                 'RP_k', 'RP_rn', 'RP_r10', 'RP_a', 'RP_o',                 
                 'RP_sparse_a', 'RP_verysparse_a',
                 'KM', 'KMPP', 'hclust', 'specc', 'mixGLasso')                      

methods_run <- c('RP_sparse_a', 'RP_verysparse_a')




### run experiment
results <- runExperiment_TCGA(fid_data = fid_run,
                              methods_arr = methods_run,
                              dir_out = dir_out,
                              doCentre = rm_group_means, 
                              parallel = run_parallel, verbose = verbose)     
closeAllConnections()




# }

