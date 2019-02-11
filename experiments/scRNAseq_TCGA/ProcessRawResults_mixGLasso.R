
ProcessRawResults_mixGLasso <- function(dir_out){
  # ProcessRawResults_mixGLasso(fid_res, setting)
  # Purpose: post-process result files from experiments
  #
  # Args:
  #   dir_out: output directory
  #
  # OUTPUT:
  #   tibble (raw results) saved to file
  
  
  ## preliminaries
  d_setting <- 'original'
  
  
  ## add mixGLasso results to data frame
  fid_mgl_sc <- './frank/results_mgl_SC.csv'
  fid_mgl_tcga <- './frank/results_mgl_TCGA.csv'
  
  dd_mgl <- as.tibble(rbind(read.csv(fid_mgl_sc), read.csv(fid_mgl_tcga)))
  colnames(dd_mgl) <- c('aRI','param_val','rep','method','data','setting')
  
  ## fix numerical issue in representation of parameter values (some trailing e-17 uncertainty)
  dd_mgl$param_val <- round(dd_mgl$param_val,6)
  
  ## map factor levels to new names
  dd_mgl <- dd_mgl %>% mutate(data = fct_recode(data, 'scRNAseq'='SingleCell'))
  #dd_mgl <- dd_mgl %>% mutate(setting = fct_recode(setting, 'fig1_n'='vary_n'))
  #dd_mgl <- dd_mgl %>% mutate(setting = fct_recode(setting, 'fig2_p'='vary_p'))
  #dd_mgl <- dd_mgl %>% mutate(setting = fct_recode(setting, 'fig3_d'='vary_d'))
  #dd_mgl <- dd_mgl %>% mutate(setting = fct_recode(setting, 'fig4_k'='vary_k'))
  
  ## fix representation of parameter value for vary_n scenarios (total n instead of n_k)
  dd_mgl$param_val[dd_mgl$data=='SC' & dd_mgl$setting=='vary_n'] <- 
    dd_mgl$param_val[dd_mgl$data=='SC' & dd_mgl$setting=='vary_n']*2
  dd_mgl$param_val[dd_mgl$data=='SC' & dd_mgl$setting=='vary_n' & dd_mgl$param_val==1728] <- 2178
  
  
  ## re-order columns
  dd_mgl <- dd_mgl[,c(5,6,3,2,4,1)]
  
  ## sort rows according to column values
  dd_mgl <- dd_mgl %>% arrange(data, setting, rep, param_val, method)
  
  ## check looks of results df
  head(dd_mgl,20)
  
  
  ## add my results
  dd_res <- tibble('data' = character(), 
                   'setting' = character(), 
                   'd_setting' = character(),
                   'param_val' = numeric(),
                   'rep' = integer(), 
                   'method' = character(), 
                   'aRI' = numeric())
  
  ## SC results
  for(f in fid_res_full[1:11]){
    if(nchar(f)==94){
      curr_setting <- substr(f,83,88)
    }else{
      curr_setting <- substr(f,83,91)
    }
    
    curr_res <- readRDS(f)
    n <- length(curr_res$aRI)
    
    curr_dd <- tibble('data' = rep('SC',n), 
                      'setting' = rep(curr_setting,n), 
                      'rep' = as.numeric(substr(curr_res$dataID,6,7)),
                      'param_val' = rep(as.numeric(curr_res$gridvalue), 
                                        length(curr_res$gridpoint)/length(curr_res$gridvalue)), 
                      'method' = as.character(curr_res$method),
                      'aRI' = as.numeric(curr_res$aRI))
    
    if(f %in% fid_res_full[1:2]){
      #note: this hack is needed to overwrite the incorrect gridvalues 
      #      (stored are n_k instead of total n) for vary_n scenarios
      curr_dd$param_val <- rep(rep(c(1,seq(2,14,2),17,21.78)*100, 
                                   each = length(unique(curr_res$method))), 
                               length(curr_res$gridpoint)/length(curr_res$gridvalue))
    }
    dd_res <- rbind(dd_res, curr_dd)
  }
  
  ## TCGA results
  for(f in fid_res_full[12:18]){
    if(nchar(f)==99){
      curr_setting <- substr(f,88,93)
    }else{
      curr_setting <- substr(f,88,96)
    }
    
    curr_res <- readRDS(f)
    n <- length(curr_res$aRI)
    
    curr_dd <- tibble('data' = rep('TCGA',n), 
                      'setting' = rep(curr_setting,n), 
                      'rep' = as.numeric(substr(curr_res$dataID,7,7)),
                      'param_val' = rep(as.numeric(curr_res$gridvalue), 
                                        length(curr_res$gridpoint)/length(curr_res$gridvalue)), 
                      'method' = as.character(curr_res$method),
                      'aRI' = as.numeric(curr_res$aRI))
    
    curr_dd$rep[curr_dd$rep==0] = 10   #hack to correct for rep=10 not correctly parsed
    
    dd_res <- dplyr::bind_rows(dd_res, curr_dd)
  }
  
  
  
  
  
  
  ## initialise results tibble
  res_tbl <- tibble('data' = character(), 
                    'setting' = character(), 
                    'd_setting' = character(),
                    'rep' = integer(), 
                    'param_val' = numeric(),
                    'method' = character(), 
                    'aRI' = numeric())
  
  
  ## read in raw results
  for(f in fid_res){
    curr_res <- readRDS(f)
    n <- length(curr_res$aRI)
    
    if(param_setting == 'n'){
      curr_dd <- tibble('data' = rep(data_id, n), 
                        'setting' = rep(param_setting, n), 
                        'd_setting' = rep(d_setting, n),
                        'rep' = as.numeric(substr(curr_res$dataID, 6,7)),
                        'param_val' = rep(c(100,200,400,600,800,1000,1200,1400,1700,2178), 
                                          each = length(unique(curr_res$method))), 
                        'method' = as.character(curr_res$method),
                        'aRI' = as.numeric(curr_res$aRI))
      
    }else{
      curr_dd <- tibble('data' = rep(data_id, n), 
                        'setting' = rep(param_setting, n), 
                        'd_setting' = rep(d_setting, n),
                        'rep' = as.numeric(substr(curr_res$dataID, 6,7)),
                        'param_val' = rep(as.numeric(curr_res$gridvalue), 
                                          length(curr_res$gridpoint)/length(curr_res$gridvalue)), 
                        'method' = as.character(curr_res$method),
                        'aRI' = as.numeric(curr_res$aRI))
    }
    
    res_tbl <- dplyr::bind_rows(res_tbl, curr_dd)
  }
  
  
  ## convert to factors
  res_tbl$data <- as.factor(res_tbl$data)
  res_tbl$setting <- as.factor(res_tbl$setting)
  res_tbl$d_setting <- as.factor(res_tbl$d_setting)
  res_tbl$method <- as.factor(res_tbl$method)
  
  
  ## save results tibble
  saveRDS(res_tbl, file = file.path(dir_out, 'results_mixglasso_raw.rds')))
}


