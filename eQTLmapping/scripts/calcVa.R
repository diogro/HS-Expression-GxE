library(GridLMM)
library(rio)
library(dplyr)

tissue = snakemake@wildcards[["tissue"]]

covariates = import(paste0("covariates/", tissue, ".tsv"))
GRM = import(paste0("GRMs/", tissue, ".cXX.txt"), header = FALSE)
colnames(GRM) = rownames(GRM) = covariates$id
phenos = import(paste0("phenotypes/", tissue, ".tsv"))
genes = import(paste0("phenotypes/", tissue, ".genes.txt"), header = FALSE)[,1]

rownames(phenos) = genes
phenos = t(phenos)
     
runNullmodel = function(current_gene, tissue, covariates, phenos, GRM){
     cache_folder = paste0('cache/',
                           tissue, '/',
                           current_gene)
     Va_file = paste0(cache_folder, '/Va.txt')
     if (file.exists(Va_file)) {
          Va = read.table(Va_file, header = FALSE)[1,1]
          return(Va)
     }
     if (!dir.exists(cache_folder)) {
          dir.create(cache_folder, recursive = TRUE)
     }
     data = covariates |>
          mutate(y = phenos[,current_gene]) |>
          as.data.frame()
     rownames(data) = data$id
     null_model = GridLMM_ML(formula = global_formulas[[tissue]],
                         data = data,
                         relmat = list(id = list(K = GRM)),
                         REML = T,
                         save_V_folder = cache_folder,
                         tolerance = 1e-3)
     Va = null_model$results[,c('id.REML')]
     V_setup = null_model$setup
     export(V_setup, paste0(cache_folder, '/V_setup.rds'))
     write.table(Va, Va_file,
                 row.names = FALSE, col.names = FALSE)
     return(Va)
}

library(furrr)
plan(multisession, workers = 32)

global_formulas = list(
     head = 
y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2 + (1|id),
     body = 
y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + (1|id))

va_list = future_map(genes, runNullmodel, 
                     tissue = tissue, 
                     covariates = covariates,
                     phenos = phenos,
                     GRM = GRM)

Va_df = data.frame(gene = genes,
                   Va = unlist(va_list))
export(Va_df, snakemake@output[["Va"]])
