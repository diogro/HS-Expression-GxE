source(here::here("Rscripts/load_install_packages.R"))

################
# Read VCFs
################

vcf.head <- "data/GXEpaper/HEAD/genotypes/VCF/Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind.vcf"
vcf.body <- "data/Genotypes_feb2020/body/body.filtered.IndMiss80.SiteMiss80MAF5.LD8.HWE.vcf.recode.vcf"

seqVCF2GDS(vcf.head, "ccm_head.gds", parallel=6L)
seqVCF2GDS(vcf.body, "ccm_body.gds", parallel=6L)

snps_head = seqOpen("ccm_head.gds") # seqClose("ccm_head.gds")
snps_body = seqOpen("ccm_body.gds") # seqClose("ccm_body.gds")


################
# Read covariates for all samples
################

covariates = bind_rows(head = import("data/GXEpaper/ALLDATA/UsefulPairs_DNA_RNA_heads_noOutliers_Mar3.20.txt"), 
                       body = import("data/GXEpaper/ALLDATA/UsefulPairs_DNA_RNA_body_noOutliers_Feb6.20.txt"),
                       .id = "tissue") |> 
    as_tibble() |>
    mutate(treatment = as.factor(treatment),
           RNAlibBatch = as.factor(RNAlibBatch),
           RNAseqBatch = as.factor(RNAseqBatch),
           egglayBatch = as.factor(egglayBatch),
           platingBatch = as.factor(platingBatch))

###############
# Read count data
###############

countdata = list(head = read.table("data/GXEpaper/ALLDATA/Allgenecounts_DNA_RNA_heads_noOutliers_Mar3.20.txt", 
                                   header=T, check.names=F),
                 body = read.table("data/GXEpaper/ALLDATA/Allgenecounts_DNA_RNA_body_noOutliers_Jun1.20.txt", 
                                   header=T, check.names=F) %>% mutate(Geneid = rownames(.)) %>% relocate(Geneid))

filterGenes <- function(countdata.norm){
    print(paste("Total genes:", nrow(countdata.norm$genes)))
    # Removing genes where max expression is less than 1 cpm
    max_cpm = apply(cpm(countdata.norm), 1, max)
    drop <- which(max_cpm < 1)
    print(length(drop))
    countdata.norm <- countdata.norm[-drop, ]

    # Removing genes with average expression less than 1 cpm
    mean_cpm = apply(cpm(countdata.norm), 1, mean)
    drop <- which(mean_cpm < 1)
    print(length(drop))
    countdata.norm <- countdata.norm[-drop, ]

    # Removing genes with expression less than 1 cpm is more than 20% of the samples
    prop_non_expressed = apply(cpm(countdata.norm), 1, \(x) sum(x < 1)) / nrow(countdata.norm$samples)
    drop <- which(prop_non_expressed > 0.2)
    print(length(drop))
    countdata.norm <- countdata.norm[-drop, ]

    print(paste("Filtered genes:", nrow(countdata.norm$genes)))
    countdata.norm
}

DGEList
y <- map(countdata, DGEList)
y <- map(y, calcNormFactors)
y.filtered <- map(y, filterGenes)


#################
# Set model matrices
#################

setModelMatrices = function(x, y, covariates){
    cov = filter(covariates, tissue == x) |> 
        filter(id %in% rownames(y[[x]]$samples)) |>
        arrange(match(id, rownames(y[[x]]$samples)))
    mod1 <- model.matrix(~0+treatment + 
                        egglayBatch + 
                        RNAseqBatch + 
                        platingBatch + 
                        RNAlibBatch, cov) 
    mod0 <- model.matrix(~0+egglayBatch + 
                        RNAseqBatch + 
                        platingBatch + 
                        RNAlibBatch, cov) 
    design <- model.matrix(~as.factor(treatment), cov) 
    list(counts = y[[x]], covariates = cov, mod1 = mod1, mod0 = mod0, design = design)
}
rnaseq_data = map(c(head = "head", body = "body"), setModelMatrices, y.filtered, covariates)

setCPMVoom <- function(data){
    data$cpm = cpm(data$counts, log=TRUE, prior.count=3)
    data$voom <- voom(data$counts, data$mod1)
    data
}
rnaseq_data = map(rnaseq_data, setCPMVoom)

setSVA = function(data){
    data$sva_cpm = sva(data$cpm, data$mod1, data$mod0, n.sv=4)
    data$sva_voom = sva(data$voom$E, data$mod1, data$mod0, n.sv=4)
    data
}
rnaseq_data = map(rnaseq_data, setSVA)

tic()
voom_svobj = sva(v$E, mod1, mod0, n.sv=5); toc()

for(i in 1:4){
    no.svs <- removeBatchEffect(cpm_data_sva, 
                                design = design, 
                                covariates = cbind(mod0, svobj$sv[,1:i]))
    cpm_pca_no.svs = pca(no.svs)
    col = covariates[match(colnames(no.svs), covariates$V1), 'treatment'] 
    png(paste0('tmp/CPM-PCA-sv', i, '.png'), width = 1080, height = 1080)
    print(pca_plot(cpm_pca_no.svs, c("C", "HS")[col]))
    dev.off()
}

library(patchwork)
library(purrr)
no.svs <- removeBatchEffect(cpm_data_sva, 
                            design = design, 
                            covariates = cbind(mod0, svobj$sv[,1:i]))
col = covariates[match(colnames(cpm_data_sva), covariates$V1), 'treatment'] 
cpm_ctrl = cpm_data_sva[, col == 1]
cpm_hs = cpm_data_sva[, col == 2]

pca_outliers = map(list(cpm_ctrl, cpm_hs), pca_plot_outliers)


png(paste0('tmp/CPM-PCA-outliers.png'), width = 2*1080, height = 1080)
pca_outliers[[1]] + ggtitle("Ctrl") + pca_outliers[[2]] + ggtitle("HS")
dev.off()

for(i in 1:5){
    no.svs <- removeBatchEffect(v$E, 
                                design = design, 
                                covariates = cbind(mod0, voom_svobj$sv[,1:i]))
    cpm_pca_no.svs = pca(no.svs)
    col = covariates[match(colnames(no.svs), covariates$V1), 'treatment'] 
    png(paste0('tmp/VOOM-PCA-sv', i, '.png'), width = 1080, height = 1080)
    print(pca_plot(cpm_pca_no.svs, c("C", "HS")[col]))
    dev.off()
}

col = covariates[match(colnames(VOOMCounts_CPM1_head_hsctrl_covfree), covariates$V1), 'treatment'] 
VOOM_ctrl = VOOMCounts_CPM1_head_hsctrl_covfree[, col == 1]
VOOM_hs = VOOMCounts_CPM1_head_hsctrl_covfree[, col == 2]

pca_outliers = map(list(VOOM_ctrl, VOOM_hs), pca_plot_outliers)


png(paste0('tmp/VOOM-PCA-outliers.png'), width = 2*1080, height = 1080)
pca_outliers[[1]] + ggtitle("Ctrl") + pca_outliers[[2]] + ggtitle("HS")
dev.off()

