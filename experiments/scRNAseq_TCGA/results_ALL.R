

## script to combine and process MCAP output

library(tidyverse)
library(mcap)



### results from latest run (incl. sparse RP methods)

fid_path <- file.path(getwd(), 'experiments', 'scRNAseq_TCGA', 'results', 'SC')
fid_list <- file.path(fid_path, list.files(fid_path))


### process raw results
dd_n_d0_1_5 <- ProcessRawResults_scRNAseq(fid_list[1:5], param_setting = 'n', 
                                          d_setting = 'zero', dir_out = NULL)
#NOTE: need to fix parameter settings ! 
dd_n_d0_6_10 <- ProcessRawResults_scRNAseq(fid_list[6:10], param_setting = 'n', 
                                           d_setting = 'zero', dir_out = NULL)
dd_n_d0 <- dplyr::bind_rows(dd_n_d0_1_5, dd_n_d0_6_10)

dd_n_1_5 <- ProcessRawResults_scRNAseq(fid_list[36:40], param_setting = 'n', 
                                          d_setting = 'original', dir_out = NULL)
#NOTE: need to fix parameter settings ! 
dd_n_6_10 <- ProcessRawResults_scRNAseq(fid_list[41:45], param_setting = 'n', 
                                           d_setting = 'original', dir_out = NULL)
dd_n <- dplyr::bind_rows(dd_n_1_5, dd_n_6_10)



dd_p <- ProcessRawResults_scRNAseq(fid_list[46:55], param_setting = 'p', 
                                      d_setting = 'original', dir_out = NULL)

dd_p_d0 <- ProcessRawResults_scRNAseq(fid_list[11:20], param_setting = 'p', 
                                   d_setting = 'zero', dir_out = NULL)

dd_large_p <- ProcessRawResults_scRNAseq(fid_list[76:85], param_setting = 'sim_p', 
                                   d_setting = 'original', dir_out = NULL, sim = TRUE)

dd_large_p_d0 <- ProcessRawResults_scRNAseq(fid_list[31:35], param_setting = 'sim_p', 
                                      d_setting = 'zero', dir_out = NULL, sim = TRUE)



dd_d <- ProcessRawResults_scRNAseq(fid_list[56:65], param_setting = 'd', 
                                   d_setting = 'original', dir_out = NULL)

dd_k <- ProcessRawResults_scRNAseq(fid_list[66:75], param_setting = 'k', 
                                   d_setting = 'original', dir_out = NULL)

dd_k_d0 <- ProcessRawResults_scRNAseq(fid_list[21:30], param_setting = 'k', 
                                      d_setting = 'zero', dir_out = NULL)


dd_sim_d_k2 <- ProcessRawResults_scRNAseq(fid_list[86:95], param_setting = 'sim_d', 
                                          d_setting = 'sim_k2', dir_out = NULL, 
                                          sim = TRUE)

dd_sim_d_k4 <- ProcessRawResults_scRNAseq(fid_list[96:100], param_setting = 'sim_d', 
                                          d_setting = 'sim_k4', dir_out = NULL, 
                                          sim = TRUE)



fid_path <- file.path(getwd(), 'experiments', 'results', 'TCGA')
fid_list <- file.path(fid_path, list.files(fid_path))
dd_k_d0_tcga <- ProcessRawResults_TCGA(fid_list[1:3], param_setting = 'k',
                                       d_setting = 'zero', dir_out = NULL)

dd_n_d0_tcga <- ProcessRawResults_TCGA(fid_list[4:6], param_setting = 'n',
                                    d_setting = 'zero', dir_out = NULL)

dd_p_d0_tcga <- ProcessRawResults_TCGA(fid_list[7:9], param_setting = 'p',
                                    d_setting = 'zero', dir_out = NULL)

dd_d_tcga <- ProcessRawResults_TCGA(fid_list[10:14], param_setting = 'd',
                                    d_setting = 'original', dir_out = NULL)

dd_k_tcga <- ProcessRawResults_TCGA(fid_list[15:19], param_setting = 'k',
                                    d_setting = 'original', dir_out = NULL)

dd_n_tcga <- ProcessRawResults_TCGA(fid_list[20:24], param_setting = 'n',
                                    d_setting = 'original', dir_out = NULL)

dd_p_tcga <- ProcessRawResults_TCGA(fid_list[25:29], param_setting = 'p',
                                    d_setting = 'original', dir_out = NULL)



## combine all parameter scenarios
dd_sc_sparse_raw <- dplyr::bind_rows(dd_n, dd_n_d0, dd_p, dd_p_d0, dd_d, dd_k, dd_k_d0,
                                     dd_large_p, dd_large_p_d0,
                                     dd_sim_d_k2, dd_sim_d_k4)
dd_sc_sparse_raw <- dd_sc_sparse_raw %>% filter(method %in% c('RP_sparse_a',
                                                              'RP_verysparse_a'))

dd_tcga_sparse_raw <- dplyr::bind_rows(dd_n_tcga, dd_p_tcga, dd_d_tcga, dd_k_tcga,
                                       dd_n_d0_tcga, dd_p_d0_tcga, dd_k_d0_tcga)
                                       




### mixglasso and old MCAP results (incl. other methods)
fid_mgl_raw <- file.path(getwd(), 'experiments', 'scRNAseq_TCGA', 'old_results', 
                         'combined_results_new_ALL_raw_v2.rds')

fid_mgl_summary <- file.path(getwd(), 'experiments', 'scRNAseq_TCGA', 'old_results', 
                             'combined_results_new_ALL_summary_v2.rds')


dd_old_raw <- readRDS(fid_mgl_raw)

dd_mgl_raw <- dd_old_raw %>% filter(method %in% c('mgl_subsample_500', 'mixglasso',
                                                   'GGMM_k', 'GGMM_o',
                                                   'GGMM_a', 'RP_a',
                                                   'hclust', 'KM', 'mclust', 
                                                   'specc', 'kmeans'))

dd_mgl_raw <- dd_mgl_raw %>% mutate('d_setting' = setting)
dd_mgl_raw <- dd_mgl_raw %>% 
                mutate(d_setting = fct_recode(d_setting, 
                                              'zero' = 'fig1_n_d0', 
                                              'zero' = 'fig2_p_d0', 
                                              'zero' = 'fig4_k_d0',
                                              'zero' = 'vig1_p_d0',
                                              'original' = 'fig1_n', 
                                              'original' = 'fig2_p',
                                              'original' = 'fig3_d', 
                                              'original' = 'fig4_k',
                                              'original' = 'vary_n', 
                                              'original' = 'vary_p',
                                              'original' = 'vary_d', 
                                              'original' = 'vary_k',
                                              'original' = 'vig1_p', 
                                              'K = 2' = 'vig2_d_k2',
                                              'K = 4' = 'vig2_d_k4'))
                       

dd_mgl_raw <- dd_mgl_raw %>% 
                mutate(data = fct_recode(data, 
                                         'scRNAseq' = 'SC',
                                         'scRNAseq' = 'SC (d=0)',
                                         'sim' = 'Sim. (k=2)',
                                         'sim' = 'Sim. (k=4)'))
dd_mgl_raw <- dd_mgl_raw %>% 
                mutate(setting = fct_recode(setting,
                                            'n' = 'fig1_n',
                                            'n' = 'fig1_n_d0',
                                            'n' = 'vary_n',
                                            'p' = 'fig2_p',
                                            'p' = 'fig2_p_d0',
                                            'p' = 'vary_p',
                                            'sim_p' = 'vig1_p',
                                            'sim_p' = 'vig1_p_d0',
                                            'd' = 'fig3_d', 
                                            'd' = 'vary_d',
                                            'sim_d' = 'vig2_d_k2',
                                            'sim_d' = 'vig2_d_k4',
                                            'k' = 'fig4_k',
                                            'k' = 'fig4_k_d0',
                                            'k' = 'vary_k'))




## combine mixglasso, previous runs and new runs
dd_raw_all <- dplyr::bind_rows(dd_sc_sparse_raw, dd_tcga_sparse_raw, dd_mgl_raw)



## tidy up factors
dd_raw_all$data <- as.factor(dd_raw_all$data)
dd_raw_all$setting <- as.factor(dd_raw_all$setting)
dd_raw_all$d_setting <- as.factor(dd_raw_all$d_setting)
dd_raw_all$method <- as.factor(dd_raw_all$method)
dd_raw_all <- dd_raw_all %>% 
                mutate(method = fct_recode(method, 
                                           'MixGLasso' = 'mixglasso',
                                           'MGL-sub500' = 'mgl_subsample_500',
                                           'MCAP-PCA' = 'GGMM_a',
                                           'MCAP-PCA-K' = 'GGMM_k',
                                           'MCAP-PCA-oracle' = 'GGMM_o',
                                           'MCAP-RP-Gauss' = 'RP_a',
                                           'MCAP-RP-Li' = 'RP_sparse_a',
                                           'MCAP-RP-Achl' = 'RP_verysparse_a',
                                           'K-means' = 'KM',
                                           'spectral' = 'specc'))


## sort rows according to column values
dd_raw_all <- dd_raw_all %>% arrange(data, setting, d_setting, rep, param_val, method)

## tidy up aRI values
dd_raw_all$aRI[dd_raw_all$aRI<0] <- 0


dd_raw_all <- dd_raw_all %>% 
                mutate(d_setting = fct_recode(d_setting, 
                                              'mean centered' = 'zero',
                                              'K = 2' = 'sim_k2',
                                              'K = 4' = 'sim_k4'))


### result summaries (mean & SEM)
## compute mean and SEM across realisations of the same parameter set for each method
dd_summary_all <- dd_raw_all %>% 
                    group_by(data, setting, d_setting, param_val, method) %>% 
                    summarise_at(vars(aRI), 
                                 funs('aRI_mean' = mean(., na.rm=TRUE), 
                                      'aRI_SEM' = ComputeSEM)) 


dd_summary_all$setting <- recode(dd_summary_all$setting,
                                 'p' = 'p (log10-scale)', 
                                 'k' = 'K')

dd_summary_all$param_val[dd_summary_all$setting == 'p (log10-scale)'] <- 
  log10(dd_summary_all$param_val[dd_summary_all$setting == 'p (log10-scale)'])
dd_summary_all$param_val[dd_summary_all$setting == 'sim_p'] <- 
  log10(dd_summary_all$param_val[dd_summary_all$setting == 'sim_p'])



### output
fid_path_out <- file.path(getwd(), 'experiments', 'scRNAseq_TCGA', 'results_processed')
fid_out_raw <- file.path(fid_path_out, 'combined_results_ALL_raw_v3.rds')
saveRDS(dd_raw_all, file = fid_out_raw)

fid_out_summary <- file.path(fid_path_out, 'combined_results_ALL_summary_v3.rds')
saveRDS(dd_summary_all, file = fid_out_summary)







