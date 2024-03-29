library(GridLMM)
library(snpStats)
library(dplyr)
library(rio)
library(poolr)
library(qvalue)

tissue = snakemake@wildcards[["tissue"]]
current_gene = snakemake@wildcards[["gene"]]

# setwd("eQTLmapping")
# tissue = "body"

Xp = read.plink(paste0("bed_files/", tissue))
X = as(Xp$genotypes,'numeric') 
rownames(X) = Xp$fam$member
colnames(X) = Xp$map$snp.name

covariates = import(paste0("covariates/", tissue, ".tsv"))
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
                  skip = which(genes == current_gene)-1, nrows = 1))
     data = covariates |>
          mutate(y = y) |>
          as.data.frame()
     rownames(data) = data$id
     V_setup = import(paste0(cache_folder, '/V_setup.rds'))
     gxe_gwas = GridLMM_GWAS(formula = global_formulas[[tissue]],
                             test_formula =  ~1,
                             reduced_formula = ~0,
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
                mutate(Trait = current_gene) |>
                rename(snp = X_ID,
                        p_main = p_value_REML,
                        effect = beta.9) |>
                select(Trait, snp, effect, p_main) |> 
                as_tibble()
        } else if (tissue == "body"){
            out_file = results |>
                mutate(Trait = current_gene) |>
                rename(snp = X_ID,
                        p_main = p_value_REML,
                        effect = beta.7) |>
                select(Trait, snp, effect, p_main) |> 
                as_tibble()
    }
     qvalues = qvalue(out_file$p_main, fdr.level = 0.01)
     out_file = out_file |> 
          mutate(q_main = qvalues$qvalues) 
     fdr_filtered = out_file |>  
          filter(qvalues$significant) 
     out_file = filter(out_file, q_main < 0.1) |>
          as_tibble()
     export(out_file,  snakemake@output[[1]])
     if(nrow(fdr_filtered) > 0){
          detection_file = paste0("detections/eqtl/", tissue, "/", current_gene, ".tsv")
          export(fdr_filtered, detection_file)
     } 
     return()
}
runGmodel(current_gene, tissue, covariates, GRM)
cat("Done!")