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
genes = import(paste0("phenotypes/", tissue, ".genes.txt"), header = FALSE)[,1]

global_formulas <- list(
          head = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2 + (1|id),
          body = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 +         (1|id)
)

runGxEmodel = function(current_gene, tissue, covariates, y, GRM){
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
     y = t(import(paste0("phenotypes/", tissue, ".tsv"), skip = which(genes == current_gene)-1, nrows = 1))
     data = covariates |>
          mutate(y = y) |>
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
                             V_setup = V_setup,
                             method = 'REML',
                             mc.cores = 2,
                             verbose = T)
     results = gxe_gwas$results
     out_file = results |>
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
          detection_file = paste0("detections/", tissue, "/", current_gene, ".tsv")
          export(out_file, detection_file)
     } 
     return(results)
}

runGxEmodel(genes[3], tissue, covariates, y, GRM)
library(furrr)
plan(multisession, workers = 32)
options(future.globals.maxSize = 1e8 * 1024)
results = future_map(genes[1:64], runGxEmodel, 
                     tissue = tissue, 
                     covariates = covariates,
                     GRM = GRM)
