library(rio)
source("Rscripts/functions.R")

countdata <- read.table("data/GXEpaper/genecounts/RawCounts_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt", h=T, check.names=F)
covariates = import("data/GXEpaper/genecounts/Covariates_forGEMMA_Jan82021.txt")
covariates$treatment <- as.factor(covariates$treatment)
covariates$RNAlibBatch <- as.factor(covariates$RNAlibBatch)
covariates$RNAseqBatch <- as.factor(covariates$RNAseqBatch)
covariates$egglayBatch <- as.factor(covariates$egglayBatch)
covariates$platingBatch <- as.factor(covariates$platingBatch)

library(edgeR)
library(limma)
library(sva)
library(tictoc)

y <- DGEList(countdata)
y <- calcNormFactors(y)
colnames(covariates)
mod1 <- model.matrix(~0+treatment + 
                     egglayBatch + 
                     RNAseqBatch + 
                     platingBatch + 
                     RNAlibBatch, covariates) 
mod0 = model.matrix(~0+egglayBatch + 
                    RNAseqBatch + 
                    platingBatch + 
                    RNAlibBatch, covariates) 
design <- model.matrix(~0+as.factor(treatment), covariates) 

cpm_data_sva = cpm(y, log=TRUE, prior.count=3)
tic()
cpm_data_sva_svobj = sva(cpm_data_sva, mod1, mod0, n.sv=4); toc()

library(patchwork)
library(purrr)

cpm_data_sva_residuals <- removeBatchEffect(cpm_data_sva, 
                                            design = design, 
                                            covariates = cbind(mod0, cpm_data_sva_svobj$sv))
col = covariates[match(colnames(cpm_data_sva_residuals), covariates$V1), 'treatment'] 
cpm_ctrl = cpm_data_sva_residuals[, col == 1]
cpm_hs = cpm_data_sva_residuals[, col == 2]

write.table(cpm_ctrl, here::here("SBM/rawData/batch/ctrl.tsv"), sep = "\t")
write.table(cpm_hs, here::here("SBM/rawData/batch/hs.tsv"), sep = "\t")


cpm_pca_no.svs = pca(cpm_data_sva_residuals)
col = covariates[match(colnames(cpm_data_sva_residuals), covariates$V1), 'treatment'] 
png(paste0('tmp/CPM-PCA.png'), width = 1080, height = 1080)
print(pca_plot(cpm_pca_no.svs, c("C", "HS")[col]))
dev.off()

pca_outliers = map(list(cpm_ctrl, cpm_hs), pca_plot_outliers)


png(paste0('tmp/CPM-PCA-outliers.png'), width = 2*1080, height = 1080)
pca_outliers[[1]] + ggtitle("Ctrl") + pca_outliers[[2]] + ggtitle("HS")
dev.off()
