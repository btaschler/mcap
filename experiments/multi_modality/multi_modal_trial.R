
library(mcap)


### prepare data
fid_path <- '//fileserver.dzne.de/taschlerb/DZNE/data/MCAP/generated_data/'
fid_list <- file.path(fid_path, list.files(fid_path))
data = readRDS(fid_list[1])

sel_ds = 10

x_big = data$xx[[sel_ds]]
x_big = x_big[,sample(ncol(data$xx[[sel_ds]]))]
y_big = data$labels[[sel_ds]]


x_big1 = x_big[y_big==1,]
x_big2 = x_big[y_big==2,]

n = nrow(x_big)
p = ncol(x_big)


pa = 50
pb = p-pa

n1 = 100
n2 = 100


x_a = rbind(x_big1[sample(nrow(x_big1), n1), 1:pa], 
            x_big2[sample(nrow(x_big2), n2), 1:pa])
y_a = c(rep(1,n1), rep(2, n2))

x_b = x_big2[sample(nrow(x_big2), n1+n2), (pa+1):p]
y_b = c(rep(1,n1), rep(2, n2))


x = cbind(x_a, x_b)
y = c(rep(1,n1), rep(2, n2))



### fit two blocks separately
mod_fit_a <- MCAPfit(xx = x_a, k = 2, 
                     projection = 'PCA',
                     true_labels = y_a, 
                     #centering_per_group = TRUE,
                     parallel = TRUE)
print(mod_fit_a$fit_gmm$aRI)



mod_fit_b <- MCAPfit(xx = x_b, k = 2, 
                     projection = 'PCA',
                     true_labels = y_b, 
                     #centering_per_group = TRUE,
                     parallel = TRUE)
print(mod_fit_b$fit_gmm$aRI)





### run MCAP on full data matrix
mod_fit_full <- MCAPfit(xx = x, k = 2, 
                        projection = 'PCA',
                        true_labels = y, 
                        parallel = TRUE)

print(mod_fit_full$fit_gmm$aRI)

# q_opt <- OptDimClusterStability(xx = x, k = 2,
#                                 method='PCA', n_grid = 5,
#                                 true_labels = y,
#                                 parallel = TRUE)$q_opt
# 
# mod_fit <- GMMwrapper(GramPCA(x, npc = q_opt)$zz,
#                            k = 2, true_labels = y)
# print(mod_fit$aRI)



### multi-modal optimisation
fit_opt <- OptDimClusterStability_multimodal(xa = x_a, xb = x_b, k = 2,
                                           method='PCA', 
                                           n_grid = 5,
                                           true_labels = y_b,
                                           parallel = TRUE)
qa_opt <- fit_opt$qa_opt
qb_opt <- fit_opt$qb_opt


mod_fit <- GMMwrapper(cbind(GramPCA(x_a, qa_opt)$zz, 
                            GramPCA(x_b, qb_opt)$zz),
                      k = 2, true_labels = y)
print(mod_fit$aRI)

cat('\n opt. proj. dim: ', qa_opt, qb_opt, ' (stab.: ', fit_opt$stab_score, ')',
    '\n oracle: ', fit_opt$qa_oracle, fit_opt$qb_opt, '\n')


## diagnostic
print(fit_opt$full_tbl, n=nrow(fit_opt$full_tbl))





## feed multi-modal solution back into MCAP
mod_fit <- MCAPfit(xx = cbind(GramPCA(x_a, qa_opt)$zz, 
                              GramPCA(x_b, qb_opt)$zz), k = 2, 
                   projection = 'PCA',
                   true_labels = y, 
                   parallel = TRUE)

print(mod_fit$fit_gmm$aRI)




### diagnostics
r=NULL
for(q in 2:15){ 
  mod_fit <- GMMwrapper(xx = GramPCA(cbind(GramPCA(x_a, qa_opt)$zz, 
                                           GramPCA(x_b, qb_opt)$zz), 
                                     npc=6)$zz,
                        k = 2, true_labels = y)
  r <- c(r,mod_fit$aRI)
  
}
print(mod_fit$aRI)


plot(2:15, r)
