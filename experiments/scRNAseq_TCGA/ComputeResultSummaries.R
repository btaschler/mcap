
ComputeResultSummaries <- function(fid_raw, dir_out){
  # ComputeResultSummaries(fid_raw)
  # Purpose: compute summaries (mean, SEM) across multiple realisations and 
  #          for different parameter scenarios
  #
  # Args:
  #   fid_raw: array of data files (.rds) containing processed clustering results
  #   dir_out: output directory
  #
  # OUTPUT:
  #   tibble (results summary) saved to file
  
  
  ### preliminaries
  res_tbl <- tibble('data' = character(), 
                    'setting' = character(), 
                    'd_setting' = character(),
                    'rep' = integer(), 
                    'param_val' = numeric(),
                    'method' = character(), 
                    'aRI' = numeric())
  
  ## combine multiple raw result files
  for(f in fid_raw){
    curr_res <- readRDS(f)
    res_tbl <- dplyr::bind_rows(res_tbl, curr_dd)
  }
  
  
  ## rename method names
  res_tbl$method <- recode(res_tbl$method,
                           'GGMM_a' = 'MCAP',
                           'GGMM_k' = 'MCAP-K', 
                           'RP_a'   = 'MCAP-RPgauss',
                           'RP_sparse_a' = 'MCAP-RPachl',
                           'RP_verysparse_a' = 'MCAP-RPli',
                           'mgl_subsample_500' = 'MGL-sub500',
                           'mixglasso' = 'MixGLasso',
                           'KM' = 'K-means',
                           'specc' = 'spectral')
  
  res_summary$d_setting <- recode(res_summary$d_setting, 
                                  'zero'   = 'd=0',
                                  'sim_k2' = 'K=2',
                                  'sim_k4' = 'K=4')
  
  
  ## reorder factors
  res_tbl$method <- factor(res_summary$method, 
                           levels = c('MCAP','MCAP-K','MCAP-RPgauss','MCAP-RPachl',
                                      'MCAP-RPli', 'MixGLasso','MGL-sub500',
                                      'K-means','hclust','spectral','mclust'))
  
  

  ### compute mean and SEM across realisations of the same parameter set for each method
  res_summary <- res_tbl %>% 
    group_by(data, setting, d_setting, param_val, method) %>% 
    summarise_at(vars(aRI), 
                 funs('aRI_mean' = mean(., na.rm=TRUE), 
                      'aRI_SEM' = ComputeSEM)) 
  

  ## save results tibble
  saveRDS(res_tbl, file = file.path(dir_out, paste0(data_id, '_SUMMARY.rds')))
}


