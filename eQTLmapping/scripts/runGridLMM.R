library(GridLMM)
library(snpStats)
library(dplyr)
library(rio)
library(poolr)
library(qvalue)

source(here::here("eQTLmapping/scripts/snpPosClassify.R"))

tissue = snakemake@wildcards[["tissue"]]
current_gene = snakemake@wildcards[["gene"]]

# setwd(here::here("eQTLmapping"))
# tissue = "body"

Xp = read.plink(paste0("bed_files/", tissue))
X = as(Xp$genotypes,'numeric') 
rownames(X) = Xp$fam$member
colnames(X) = Xp$map$snp.name

covariates = import(paste0("covariates/", tissue, ".tsv")) |>
     mutate(egglayBatch = as.character(egglayBatch), 
            RNAseqBatch = as.character(RNAseqBatch),
            platingBatch = as.character(platingBatch),
            RNAlibBatch = as.character(RNAlibBatch))
GRM = import(paste0("GRMs/", tissue, ".cXX.txt"), header = FALSE)
colnames(GRM) = rownames(GRM) = covariates$id
genes = import(paste0("phenotypes/", tissue, ".genes.txt"), header = FALSE)[,1]

# current_gene = genes[100]

global_formulas <- list(
          head = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2 + (1|id),
          body = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 +         (1|id)
)

runGmodel = function(current_gene, tissue, covariates, GRM){
     cache_folder = paste0('cache/',
                           tissue, '/',
                           current_gene)
     y = t(import(paste0("phenotypes/", tissue, ".tsv"), 
                  skip = which(genes == current_gene), 
                  nrows = 1))
     data = covariates |>
          mutate(y = y) |>
          as.data.frame()
     rownames(data) = data$id
     V_setup = import(paste0(cache_folder, '/V_setup.rds'))
     gxe_gwas = GridLMM_GWAS(formula = global_formulas[[tissue]],
                             test_formula =  ~ 1,
                             reduced_formula = ~ 0,
                             data = data,
                             X = X,
                             X_ID = 'id',
                             relmat = list(id = list(K = GRM)),
                             V_setup = V_setup,
                             method = 'REML',
                             mc.cores = 1,
                             verbose = T)
     results = gxe_gwas$results
     if(tissue == "head"){
            out_file = results |>
                dplyr::mutate(Trait = current_gene) |>
                dplyr::rename(snp = X_ID,
                        p_main = p_value_REML,
                        effect = beta.19) |>
                dplyr::select(Trait, snp, effect, p_main) |> 
                as_tibble()
        } else if (tissue == "body"){
            out_file = results |>
                dplyr::mutate(Trait = current_gene) |>
                dplyr::rename(snp = X_ID,
                        p_main = p_value_REML,
                        effect = beta.11) |>
                dplyr::select(Trait, snp, effect, p_main) |> 
                as_tibble()
    }
     qvalues = qvalue(out_file$p_main, fdr.level = 0.01)
     out_file = out_file |> 
          dplyr::mutate(q_main = qvalues$qvalues, 
                 cis = as.numeric(snp %in% getCisSnps(current_gene)))
     pvalue_filtered = out_file |>  
          dplyr::filter(p_main < 0.01) 
     export(pvalue_filtered,  snakemake@output[[1]])
     return(0)
}
runGmodel(current_gene, tissue, covariates, GRM)
cat("Done!")