
ProcessRawResults_TCGA <- function(fid_res, param_setting, d_setting, 
                                   dir_out = NULL){
  # ProcessRawResults_TCGA(fid_res, setting)
  # Purpose: post-process result files from experiments
  #
  # Args:
  #   fid_res: array of data files (.rds) containing clustering results
  #   param_setting: scenario / parameter setting ('n', 'p', 'd', 'k')
  #   d_setting: separation setting for group means ('original', 'zero')
  #   dir_out: (optional) output directory [defaultL: NULL]
  #   sim: logical, indicate simulation data [default: FALSE]
  #
  # OUTPUT:
  #   tibble (raw results) saved to file
  
  
  ## preliminaries
  stopifnot(param_setting %in% c('n','p','d','k'))
  stopifnot(d_setting %in% c('original', 'zero'))
  
  data_id <- 'TCGA'
  
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
    curr_res <- as.tibble(apply(readRDS(f), 2, unlist))
    #curr_res <- curr_res[1:40,]    # hack for results from ID1006-ID1010 !
    n <- length(curr_res$aRI)
    
    if(param_setting == 'n'){
      curr_dd <- tibble('data' = rep(data_id,n), 
                        'setting' = rep(param_setting,n), 
                        'd_setting' = rep(d_setting, n),
                        'rep' = as.numeric(substr(curr_res$dataID,7,7)),
                        'param_val' = rep(c(149,298,447,596,745,894,1043,1192,1341,1490), 
                                          each = length(unique(curr_res$method))), 
                        'method' = as.character(curr_res$method),
                        'aRI' = as.numeric(curr_res$aRI))
    }else{
      curr_dd <- tibble('data' = rep(data_id,n), 
                        'setting' = rep(param_setting,n), 
                        'd_setting' = rep(d_setting, n),
                        'rep' = as.numeric(substr(curr_res$dataID,7,7)),
                        'param_val' = rep(as.numeric(curr_res$gridvalue), 
                                          length(curr_res$gridpoint)/length(curr_res$gridvalue)), 
                        'method' = as.character(curr_res$method),
                        'aRI' = as.numeric(curr_res$aRI))
    }
    
    curr_dd$rep[curr_dd$rep==0] = 10   #hack to correct for rep=10 not correctly parsed
    
    res_tbl <- dplyr::bind_rows(res_tbl, curr_dd)
  }
  
  
  ## convert to factors
  res_tbl$data <- as.factor(res_tbl$data)
  res_tbl$setting <- as.factor(res_tbl$setting)
  res_tbl$d_setting <- as.factor(res_tbl$d_setting)
  res_tbl$method <- as.factor(res_tbl$method)
  
  
  ## save results tibble
  if(!is.null(dir_out)){
    saveRDS(res_tbl, file = file.path(dir_out, paste0(data_id, 
                                                      '_vary_', param_setting,
                                                      '_d_', d_setting,
                                                      '_raw.rds')))
  }
  return(res_tbl)
}


