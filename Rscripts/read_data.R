library(vcfR) 
library(rio)

source("Rscripts/functions.R")

raw_vcf = read.vcfR("data/GXEpaper/genotypes/VCF/Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind.vcf")


library("SNPRelate")
vcf.fn<-"data/GXEpaper/genotypes/VCF/Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")
ccm_pca<-snpgdsPCA(genofile)

covariates = import("data/GXEpaper/genecounts/Covariates_forGEMMA_Jan82021.txt")
covariates[match(ccm_pca$sample, covariates$V1), 'treatment']

png('test.png', width = 1080, height = 1080)
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],col=covariates[match(ccm_pca$sample, covariates$V1), 'treatment'], pch=19)
dev.off()

x = ccm_pca$eigenval[is.finite(ccm_pca$eigenval)]

png('scree.png', width = 1080, height = 1080)
plot(x/sum(x) * 100, pch=19)
dev.off()

GRM_PCA = import("data/GXEpaper/genotypes/GRM_PCA/Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind.eigenvec", format = "tsv")
colnames(GRM_PCA)

png('test.png', width = 1080, height = 1080)
plot(GRM_PCA$PC1, GRM_PCA$PC2,col=covariates[match(paste(GRM_PCA$FID, GRM_PCA$IID, sep="_"), covariates$V1), 'treatment'], pch=19)
dev.off()


VOOMCounts_CPM1_head_hsctrl_covfree = import("data/GXEpaper/genecounts/VOOMCounts_CPM1_head_hsctrl_covfree_4svs_CORRECT_Jan8.21.txt", format= "tsv")
head(VOOMCounts_CPM1_head_hsctrl_covfree[1:3])
dim(VOOMCounts_CPM1_head_hsctrl_covfree)
rownames(VOOMCounts_CPM1_head_hsctrl_covfree) = VOOMCounts_CPM1_head_hsctrl_covfree[,1]
VOOMCounts_CPM1_head_hsctrl_covfree$V1 = NULL
VOOM_pca = pca(VOOMCounts_CPM1_head_hsctrl_covfree)

col = covariates[match(colnames(VOOMCounts_CPM1_head_hsctrl_covfree), covariates$V1), 'treatment'] 
png('test1.png', width = 1080, height = 1080)
pca_plot(VOOM_pca, c("C", "HS")[col])
dev.off()


countdata <- read.table("data/GXEpaper/genecounts/RawCounts_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt", h=T, check.names=F)
coldata <- read.table("data/GXEpaper/genecounts/Covariates_forGEMMA_Jan82021.txt",h=T)
coldata$treatment <- as.factor(coldata$treatment)
coldata$RNAlibBatch <- as.factor(coldata$RNAlibBatch)
coldata$RNAseqBatch <- as.factor(coldata$RNAseqBatch)
coldata$egglayBatch <- as.factor(coldata$egglayBatch)
coldata$platingBatch <- as.factor(coldata$platingBatch)

y <- DGEList(countdata)
y <- calcNormFactors(y)
design <- model.matrix(~0+coldata$treatment) # 2 groups
v <- voom(y, design)
covariates_design <- model.matrix(~coldata$egglayBatch+coldata$RNAseqBatch+coldata$platingBatch+coldata$RNAlibBatch+coldata$sv1+coldata$sv2+coldata$sv3+coldata$sv4) 
no.svs <- removeBatchEffect(cpm(y, log=TRUE, prior.count=3), design = design, covariates = covariates_design)

VOOM_pca_no.svs = pca(no.svs)
col = covariates[match(colnames(no.svs), covariates$V1), 'treatment'] 
png('test.png', width = 1080, height = 1080)
pca_plot(VOOM_pca_no.svs, c("C", "HS")[col])
dev.off()
