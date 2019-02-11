



### preliminaries -------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(latex2exp)



### get data
fid_summary <- file.path(getwd(), 'experiments', 'scRNAseq_TCGA', 'results_processed', 
                         'combined_results_ALL_summary_v3.rds')
dd_summary <- readRDS(fid_summary)


## reorder factor levels
dd_summary$method <- factor(dd_summary$method, levels = c('MCAP-PCA',
                                                          'MCAP-PCA-K',
                                                          'MCAP-RP-Gauss',
                                                          'MCAP-RP-Achl',
                                                          'MCAP-RP-Li',
                                                          'MixGLasso',
                                                          'MGL-sub500',
                                                          'K-means',
                                                          'hclust',
                                                          'spectral',
                                                          'mclust',
                                                          'MCAP-PCA-oracle'))


select_method <- c('MCAP-PCA','MCAP-PCA-K', 'MCAP-RP-Gauss','MCAP-RP-Li',
                   'MCAP-RP-Achl', 'K-means', 'hclust', 'spectral', 
                   'MixGLasso','MGL-sub500')        
select_method_mclust <- c('MCAP-PCA','MCAP-PCA-K', 'MCAP-RP-Gauss','MCAP-RP-Li',
                          'MCAP-RP-Achl', 'K-means', 'hclust', 'spectral', 
                          'MixGLasso','MGL-sub500', 'mclust')  
select_method_no_mixglasso <- c('MCAP-PCA','MCAP-PCA-K', 'MCAP-RP-Gauss','MCAP-RP-Li',
                                'MCAP-RP-Achl', 'K-means', 'hclust', 'spectral', 
                                'MGL-sub500')        


### final figures - plotting v.3 ----------------------------------------------
## save all figures to pdf
#pdf('./plots/figures_ALL_JMLR.pdf', width=10, height=5)


## fig: SC - vary n -----------------------------------------------------------
dd <- 'scRNAseq'
pd <- position_dodge(0.030)
fplot <- paste0('./experiments/scRNAseq_TCGA/plots/fig_',dd,'_vary_n_JMLR_v3.eps')
#pdf(paste0('./experiments/scRNAseq_TCGA/plots/fig_',dd,'_vary_n_JMLR_v3.pdf'), width=10, height=5)
print(dd_summary %>% filter(setting %in% c('n'), 
                            data==dd, method %in% select_method) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        facet_wrap(~d_setting, scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,8,19,7,10,6)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,1,1,3,2,2,2,2))+
        theme_light() +
        theme(legend.position='bottom', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x='n', y='adjusted Rand Index'))
#ggsave(fplot, width = 10, height = 5)
#dev.off()   #for pdf plotting
#note: print() is needed when ggplot is inside a for loop



## fig: SC - vary p -----------------------------------------------------------
dd <- 'scRNAseq'
pd <- position_dodge(0.030)
fplot <- paste0('./experiments/scRNAseq_TCGA/plots/fig_',dd,'_vary_p_JMLR_v3.eps')
print(dd_summary %>% filter(setting %in% c('p (log10-scale)'), 
                            d_setting %in% c('original', 'mean centered'),
                            data==dd, method %in% select_method_mclust) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        facet_grid(cols=vars(d_setting), scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,6,8,19,7,10,13)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,3,1,1,3,2,2,2,2,1))+
        theme_light() +
        theme(legend.position='bottom', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x=TeX('log_{10} p'), y='adjusted Rand index'))
#ggsave(fplot, width = 10, height = 5)



## fig: SC - vary d -----------------------------------------------------------
dd <- 'scRNAseq'
pd <- position_dodge(0.01)
fplot <-  paste0('./experiments/scRNAseq_TCGA/plots/fig_',dd,'_vary_d_JMLR_v3.eps')
print(dd_summary %>% filter(setting %in% c('d'), 
                            data==dd, method %in% select_method) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        #facet_wrap(~d_setting, scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,8,19,7,10,6)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,1,1,3,2,2,2,2))+
        theme_light() +
        theme(legend.position='right', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x='d', y='adjusted Rand Index'))
#ggsave(fplot, width = 8, height = 3.6)


## fig: SC - vary k -----------------------------------------------------------
dd <- 'scRNAseq'
pd <- position_dodge(0.1)
fplot <- paste0('./experiments/scRNAseq_TCGA/plots/fig_',dd,'_vary_k_JMLR_v3.eps')
print(dd_summary %>% filter(setting %in% c('K'), 
                            data==dd, method %in% select_method_mclust) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        facet_wrap(~d_setting, scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,6,8,19,7,10,13)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,3,1,1,3,2,2,2,2,1))+
        theme_light() +
        theme(legend.position='bottom', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x='K', y='adjusted Rand Index'))
#ggsave(fplot, width = 10, height = 5)


## fig: sim - vary large p ----------------------------------------------------
pd <- position_dodge(0.030)
fplot <- paste0('./experiments/scRNAseq_TCGA/plots/fig_sim_vary_large_p_JMLR_v3.eps')
print(dd_summary %>% filter(setting %in% c('sim_p'), 
                        method %in% select_method_mclust) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        facet_grid(cols=vars(d_setting), scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,19,7,10,13)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,1,2,2,2,1))+
        theme_light() +
        theme(legend.position='bottom', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x=TeX('log_{10} p'), y='adjusted Rand index'))
#ggsave(fplot, width = 10, height = 5)



## fig: sim (KM) - vary d -----------------------------------------------------
pd <- position_dodge(0.008)
fplot <- paste0('./experiments/scRNAseq_TCGA/plots/fig_sim_vary_d_JMLR_v3.eps')
print(dd_summary %>% filter(setting %in% c('sim_d'), 
                            method %in% select_method_mclust) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        facet_wrap(~d_setting, scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,19,7,10,13)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,1,1,2,2,2,2))+
        theme_light() +
        theme(legend.position='bottom', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x='d', y='adjusted Rand Index'))
#ggsave(fplot, width = 10, height = 5)




## fig: TCGA - vary n ---------------------------------------------------------
dd <- 'TCGA'
pd <- position_dodge(0.06)
fplot <- paste0('./experiments/scRNAseq_TCGA/plots/fig_',dd,'_vary_n_JMLR_v3.eps')
print(dd_summary %>% filter(setting %in% c('n'), 
                            data==dd, method %in% select_method_no_mixglasso) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        facet_wrap(~d_setting, scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,8,19,7,10,6)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,1,1,3,2,2,2,2))+
        theme_light() +
        theme(legend.position='bottom', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x='n', y='adjusted Rand Index'))
#ggsave(fplot, width = 10, height = 5)



## fig: TCGA - vary p -----------------------------------------------------------
dd <- 'TCGA'
pd <- position_dodge(0.030)
fplot <- paste0('./experiments/scRNAseq_TCGA/plots/fig_',dd,'_vary_p_JMLR_v3.eps')
print(dd_summary %>% filter(setting %in% c('p (log10-scale)'), 
                            d_setting %in% c('original', 'mean centered'),
                            data==dd, method %in% select_method_mclust) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        facet_grid(cols=vars(d_setting), scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,6,8,19,7,10,13)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,3,1,1,3,2,2,2,2,1))+
        theme_light() +
        theme(legend.position='bottom', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x=TeX('log_{10} p'), y='adjusted Rand index'))
#ggsave(fplot, width = 10, height = 5)


## fig: TCGA - vary d ---------------------------------------------------------
dd <- 'TCGA'
pd <- position_dodge(0.01)
fplot <- paste0('./experiments/scRNAseq_TCGA/plots/fig_',dd,'_vary_d_JMLR_v3.eps')
print(dd_summary %>% filter(setting %in% c('d'), 
                            data==dd, method %in% c('MCAP-PCA','MCAP-PCA-K', 
                                                    'MCAP-RP-Gauss','MCAP-RP-Li',
                                                    'MCAP-RP-Achl', 'K-means', 
                                                    'hclust', 'MGL-sub500')) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        #facet_wrap(~d_setting, scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,8,19,7,10,6)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,1,1,3,2,2,2,2))+
        theme_light() +
        theme(legend.position='right', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x='d', y='adjusted Rand Index'))
#ggsave(fplot, width = 8, height = 3.6)


## fig: TCGA - vary k -----------------------------------------------------------
dd <- 'TCGA'
pd <- position_dodge(0.12)
fplot <- paste0('./experiments/scRNAseq_TCGA/plots/fig_',dd,'_vary_k_JMLR_v3.eps')
print(dd_summary %>% filter(setting %in% c('K'), 
                            data==dd, method %in% select_method_mclust) %>%
        ggplot(aes(x=param_val, y=aRI_mean, colour=method, label=method, 
                   group=method, shape=method, linetype=method)) +
        facet_wrap(~d_setting, scales='free_x') +
        geom_errorbar(aes(ymin=aRI_mean-aRI_SEM, ymax=aRI_mean+aRI_SEM), 
                      width=0.004, position=pd) + 
        geom_line(position = pd) +
        geom_point(position = pd, size=2) + 
        scale_colour_hue(name='Method',     #Legend label
                         l=60) +            #darker colors: e.g lightness=40
        scale_shape_manual(name='Method',values=c(15,15,17,17,17,6,8,19,7,10,13)) +
        scale_linetype_manual(name='Method', values=c(1,3,3,3,1,3,1,1,3,2,2,2,2,1))+
        theme_light() +
        theme(legend.position='bottom', legend.box = 'horizontal') +
        theme(legend.margin=margin(t = -0.2, unit='cm')) +        #reduce margin 
        theme(plot.margin = margin(.2, .2, .2, .2, "cm")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        ylim(0,1) +
        labs(x='K', y='adjusted Rand Index'))
#ggsave(fplot, width = 10, height = 5)


#dev.off()
