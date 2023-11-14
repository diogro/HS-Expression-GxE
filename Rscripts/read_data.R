source(here::here("Rscripts/load_install_packages.R"))

################
# Read VCFs
################

vcf.head <- "data/GXEpaper/HEAD/genotypes/VCF/Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind.vcf"
vcf.body <- "data/Genotypes_feb2020/body/body.filtered.IndMiss80.SiteMiss80MAF5.LD8.HWE.vcf.recode.vcf"

seqVCF2GDS(vcf.head, "cache/ccm_head.gds", parallel=6L)
seqVCF2GDS(vcf.body, "cache/ccm_body.gds", parallel=6L)

snps_head = seqOpen("cache/ccm_head.gds") # seqClose("ccm_head.gds")
snps_body = seqOpen("cache/ccm_body.gds") # seqClose("ccm_body.gds")


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

correctModelMatrix = function(X){
    QR <- qr(crossprod(X))                 # Get the QR decomposition
    vars <- QR$pivot[seq_len(QR$rank)]     # Variable numbers
    X = X[,vars]
    X
}
setModelMatrices = function(current_tissue, y, covariates){
    cov = filter(covariates, tissue == current_tissue) |> 
        filter(id %in% rownames(y[[current_tissue]]$samples)) |>
        arrange(match(id, rownames(y[[current_tissue]]$samples)))
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
    list(counts = y[[current_tissue]], 
         tissue = current_tissue,
         covariates = cov, 
         mod1 = correctModelMatrix(mod1),
         mod0 = correctModelMatrix(mod0), 
         design = design)
}
rnaseq_data = map(c(head = "head", body = "body"), setModelMatrices, y.filtered, covariates)

setCPMVoom <- function(data){
    data$cpm = cpm(data$counts, log=TRUE, prior.count=3)
    data$voom <- voom(data$counts, data$mod1)
    data
}
rnaseq_data = map(rnaseq_data, setCPMVoom)

setSVAnum = function(data){
    data$n.sva <- list(cpm  = num.sv(data$cpm, data$mod1, method="leek"),
                       voom = num.sv(data$voom$E, data$mod1, method="leek"))
    data
}
rnaseq_data = map(rnaseq_data, setSVAnum)

setSVA = function(data){
    data$sva <- list(cpm = sva(data$cpm, data$mod1, data$mod0, 
                               n.sv=data$n.sva$cpm, method = "two-step"),
                     voom = sva(data$voom$E, data$mod1, data$mod0, 
                                n.sv=data$n.sva$voom, method = "two-step"))
    data
}
rnaseq_data = map(rnaseq_data, setSVA)

makeResiduals <- function(data){
    covariates <- data$covariates
    col = pull(covariates, treatment)
    no.svs <- removeBatchEffect(data$cpm, 
                                design = data$design, 
                                covariates = cbind(data$mod0, data$sva$cpm$sv))
    cpm_pca_no.svs = pca(no.svs)
    png(paste0('tmp/CPM-PCA-', data$tissue, '-sva-corrected.png'), width = 1080, height = 1080)
        print(pca_plot(cpm_pca_no.svs, c("C", "HS")[col]))
    dev.off()
    data$cpm_residuals = no.svs
    rownames(data$cpm_residuals) = rownames(data$cpm) = rownames(data$voom$E) = data$counts$genes[,1]
    data
}
rnaseq_data = map(rnaseq_data, makeResiduals)

export(rnaseq_data, "cache/rnaseq_all.rds")
