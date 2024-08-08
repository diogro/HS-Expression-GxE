library(dglm)
library(snpStats)
library(dplyr)
library(rio)
library(poolr)
library(qvalue)

source(here::here("eQTLmapping/scripts/snpPosClassify.R"))
source(here::here("Rscripts/functions.R"))

tissue = snakemake@wildcards[["tissue"]]
current_gene = snakemake@wildcards[["gene"]]

# setwd(here::here("eQTLmapping"))
# tissue = "body"
# current_gene = "FBgn0000008"

Xp = read.plink(paste0("bed_files/", tissue))
X = as(Xp$genotypes,'numeric') 
rownames(X) = Xp$fam$member
colnames(X) = Xp$map$snp.name

covariates = import(paste0("covariates/", tissue, ".tsv")) |>
     mutate(egglayBatch = as.character(egglayBatch), 
            RNAseqBatch = as.character(RNAseqBatch),
            platingBatch = as.character(platingBatch),
            RNAlibBatch = as.character(RNAlibBatch))
genes = import(paste0("phenotypes/", tissue, ".genes.txt"), header = FALSE)[,1]

# current_gene = genes[1]

global_formulas <- list(
          head = 
      "~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2 + snp",
          body = 
      "~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + snp"
)

local_formula = global_formulas[[tissue]]

rundGLM = function(current_gene, tissue, covariates){

    y = t(import(paste0("phenotypes/", tissue, ".tsv"), 
                 skip = which(genes == current_gene),
                 nrows = 1))
    data = covariates |>
          mutate(y = y) |>
          as.data.frame()

    run_dGLM_SNP = function(current_snp){
        data = mutate(data, snp = X[,current_snp])
        m = dglm(as.formula(paste("y", local_formula)), 
                 as.formula(local_formula),
            data = data, 
            dlink = "log",
            verbose = FALSE)
        s = summary(m)
        df = data.frame(rbind(s$coefficients["snp",] , 
                                s$dispersion.summary$coefficient["snp",]), 
                name = c("mean", "dispersion"), 
                Geneid = current_gene,
                SNP = colnames(X)[current_snp])
        names(df) = c("estimate", "std.error", "z.value", "p.value", "variable", "Geneid", "SNP")
        if(df$p.value[2] < 0.01){
            return(df)
        } else {
            #print(paste("no: ", colnames(X)[current_snp]))
            NULL
        }
    }
    map_df(seq(1, ncol(X)), run_dGLM_SNP, .progress = 'SNPs')
}

dGLM_results = rundGLM(current_gene, tissue, covariates)
export(dGLM_results,  snakemake@output[[1]])