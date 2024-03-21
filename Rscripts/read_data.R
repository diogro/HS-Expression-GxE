source(here::here("Rscripts/functions.R"))

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

target_gene = "FBgn0038484"

covariates = bind_rows(head = import("data/GXEpaper/ALLDATA/UsefulPairs_DNA_RNA_heads_noOutliers_Mar3.20.txt"), 
                       body = import("data/GXEpaper/ALLDATA/UsefulPairs_DNA_RNA_body_noOutliers_Feb6.20.txt"),
                       .id = "tissue") |> 
    as_tibble() |>
    mutate(treatment = as.factor(treatment-1),
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

filterGenes <- function(countdata.norm, current_tissue, covariates){
    print(paste("Total genes:", nrow(countdata.norm$genes)))
    
    # Removing genes with average expression less than 1 cpm
    mean_cpm = apply(cpm(countdata.norm), 1, mean)
    drop <- which(mean_cpm < 5)
    print(length(drop))
    countdata.norm <- countdata.norm[-drop, ]

    # Removing genes with mean expression less than 10 counts
    mean_count = apply(countdata.norm$counts, 1, mean)
    drop <- which(mean_count < 10)
    print(length(drop))
    if(length(drop) > 0){
        countdata.norm <- countdata.norm[-drop, ]
    }

    cov = filter(covariates, tissue == current_tissue) |> 
              filter(id %in% rownames(countdata.norm$samples)) |>
              arrange(match(id, rownames(countdata.norm$samples)))

    # Removing genes with expression less than 1 cpm in more than 20% of the samples
    prop_non_expressed_hs   = apply(cpm(countdata.norm[,cov$treatment == 1]), 1, \(x) sum(x < 1)) / sum(cov$treatment == 1)
    prop_non_expressed_ctrl = apply(cpm(countdata.norm[, cov$treatment == 0]), 1, \(x) sum(x < 1)) / sum(cov$treatment == 0)
    drop <- which(prop_non_expressed_ctrl > 0.2 & prop_non_expressed_hs > 0.2)
    print(length(drop))
    countdata.norm <- countdata.norm[-drop, ]

    print(paste("Filtered genes:", nrow(countdata.norm$genes)))
    countdata.norm
}

y <- map(countdata, DGEList)
y <- map(y, calcNormFactors)
y.filtered <- map2(y, names(y), filterGenes, covariates)

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
    mod1 <- model.matrix(~1+treatment + 
                        egglayBatch + 
                        RNAseqBatch + 
                        platingBatch + 
                        RNAlibBatch, cov) 
    mod0 <- model.matrix(~1+egglayBatch + 
                        RNAseqBatch + 
                        platingBatch + 
                        RNAlibBatch, cov) 
    design <- model.matrix(~as.factor(treatment), cov) 
    colnames(design)[1:2] <- c("Intercept","hs")
    colnames(mod1)[1:2] <- c("Intercept","hs")
    rownames(y[[current_tissue]]$counts) <- y[[current_tissue]]$genes[,1]
    list(counts = y[[current_tissue]], 
         tissue = current_tissue,
         covariates = cov, 
         mod1 = correctModelMatrix(mod1),
         mod0 = correctModelMatrix(mod0), 
         design = design)
}
rnaseq_data = map(c(head = "head", body = "body"), setModelMatrices, y.filtered, covariates)

setCPMVoom <- function(data){
    data$cpm  <-  cpm(data$counts, log=TRUE)
    data$voom <- voom(data$counts, data$mod1)
    data$l2c  <- log2(data$counts$counts + 1)
    data
}
rnaseq_data = map(rnaseq_data, setCPMVoom)

# setSVAnum = function(data){
#     counts_list = list(cpm = data$cpm, 
#                        voom = data$voom$E, 
#                        l2c = data$l2c)
#     data$n.sva <- future_map(counts_list, num.sv, data$mod1, method="leek")
#     data
# }
# future::plan(list(future::tweak(future::multisession, workers = 2), 
#                   future::tweak(future::multisession, workers = 3)))
# rnaseq_data = future_map(rnaseq_data, setSVAnum)
# rnaseq_data |> map("n.sva")
rnaseq_data$head$n.sva = list(cpm = 2, voom = 2, l2c = 2)
rnaseq_data$body$n.sva = list(cpm = 0, voom = 0, l2c = 1)

setSVA = function(data){
     counts_list = list(cpm = data$cpm, 
                        voom = data$voom$E, 
                        l2c = data$l2c)
     
     data$sva <- future_map2(counts_list[data$n.sva!=0], data$n.sva[data$n.sva!=0],
                             \(x, y) sva(x, data$mod1, data$mod0, n.sv=y, method = "two-step"))
     data
}
options(future.globals.maxSize = 1e8 * 1024)
rnaseq_data = map(rnaseq_data, setSVA)

getBatchResiduals <- function(x, label, tissue, design, mod0, sva = NULL, treatment){
    no.batch <- removeBatchEffect(x, 
                                  design = design, 
                                  covariates = cbind(mod0, sva))
    pca_no.batch = pca(no.batch)
    png(paste0('tmp/PCA-batch-corrected-', tissue, '-', label, '.png'), width = 1080, height = 1080)
        print(pca_plot(pca_no.batch, c("C", "HS")[treatment]))
    dev.off()
    rownames(no.batch) = data$counts$genes[,1]
    no.batch
}
makeResiduals <- function(data){
    covariates <- data$covariates
    counts_list = list(cpm = data$cpm, 
                       voom = data$voom$E, 
                       l2c = data$l2c)
    data$batch_residuals <- future_map2(counts_list, names(counts_list), 
                            \(x, y) getBatchResiduals(x, y, data$tissue, 
                                                      data$design, data$mod0,
                                                      treatment = pull(covariates, treatment)))  
    data
}
options(future.globals.maxSize = 1e8 * 1024)
rnaseq_data = future_map(rnaseq_data, makeResiduals)

setMouthwash = function(data){
     counts_list = list(cpm = data$cpm, 
                        voom = data$voom$E, 
                        l2c = data$l2c)
     
     data$mwash <- future_map2(counts_list[data$n.sva!=0], data$n.sva[data$n.sva!=0], 
                             \(x, y) mouthwash(Y = t(x), X = data$mod1, k = y, cov_of_interest = 2,
                                               include_intercept = FALSE))
     data
}
rnaseq_data = future_map(rnaseq_data, setMouthwash)

makeResidualsMwash <- function(data){
    covariates <- data$covariates
    counts_list = list(cpm = data$cpm, 
                       voom = data$voom$E, 
                       l2c = data$l2c)
    data$mwash_residuals <- future_map2(counts_list[data$n.sva!=0], names(counts_list)[data$n.sva!=0], 
                            \(x, y) getBatchResiduals(x, paste0(y, "-mwash"), data$tissue, 
                                                      data$design, data$mod0, sva = data$mwash[[y]]$Zhat,
                                                      treatment = pull(covariates, treatment)))  
    data
}
rnaseq_data = future_map(rnaseq_data, makeResidualsMwash)

export(rnaseq_data, affix_date("cache/rnaseq_all.rds"))
export(covariates, affix_date("cache/covariates.rds"))

rnaseq_data <- import("cache/rnaseq_all.rds")
covariates  <- import("cache/covariates.rds")
