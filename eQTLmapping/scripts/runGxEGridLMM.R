library(GridLMM)
library(snpStats)
library(tidyverse)
library(rio)
library(poolr)

tissue = snakemake@wildcards[["tissue"]]

setwd("eQTLmapping")
tissue = "head"

Xp = read.plink(paste0("bed_files/", tissue))
X = as(Xp$genotypes,'numeric') 
rownames(X) = Xp$fam$member
colnames(X) = Xp$map$snp.name

covariates = import(paste0("covariates/", tissue, ".tsv"))
GRM = import(paste0("GRMs/", tissue, ".cXX.txt"), header = FALSE)
colnames(GRM) = rownames(GRM) = covariates$id
phenos = import(paste0("phenotypes/", tissue, ".tsv"))
genes = import(paste0("phenotypes/", tissue, ".genes.txt"), header = FALSE)[,1]

rownames(phenos) = genes
phenos = t(phenos)

global_formulas <- list(
          head = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2 + (1|id),
          body = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 +         (1|id)
)

runGxEmodel = function(current_gene, tissue, covariates, phenos, GRM){
     cache_folder = paste0('cache/',
                           tissue, '/',
                           current_gene)
     results_file = paste0(cache_folder, '/results.tsv')
     if (file.exists(results_file)) {
          results = import(results_file)
          return(results)
     }
     if (!dir.exists(cache_folder)) {
          dir.create(cache_folder, recursive = TRUE)
     }
     data = covariates |>
          mutate(y = phenos[,current_gene]) |>
          as.data.frame()
     rownames(data) = data$id
     V_setup = import(paste0(cache_folder, '/V_setup.rds'))
     gxe_gwas = GridLMM_GWAS(formula = global_formulas[[tissue]],
                              test_formula =  ~1 + treatment,
                              reduced_formula = ~1,
                              data = data,
                              X = X,
                              X_ID = 'id',
                              relmat = list(id = list(K = GRM)),
                              V_setup = V_setup,  # This tells the function to reuse existing matrices
                              method = 'REML',
                              mc.cores = 2,
                              verbose = T)
     out_file = gxe_gwas$results |>
          mutate(Trait = current_gene) |>
          rename(snp = X_ID,
                 p_main = p_value_REML.1,
                 p_gxe = p_value_REML.2) |>
          select(Trait, snp,  p_main, p_gxe) |>
          mutate(p.adj_gxe = p.adjust(p_gxe, method = 'BY')) |>
          filter(p.adj_gxe < 0.1) |>
          as_tibble()
     export(out_file, results_file)
     if(nrow(out_file) > 0){
          detection_file = paste0("detections/", current_gene, ".tsv")
          export(out_file, detection_file)
     } 
     return(gxe_gwas$results)
}
results = runGxEmodel(genes[2], tissue, covariates, phenos, GRM)

library(furrr)
plan(multisession, workers = 32)

results = future_map(genes, runGxEmodel, 
                     tissue = tissue, 
                     covariates = covariates,
                     phenos = phenos,
                     GRM = GRM)
