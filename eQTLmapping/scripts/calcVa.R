library(GridLMM)
library(rio)
library(dplyr)

# Load data replacing head with tissue of interest
load_data = function(tissue){
     covariates = import(paste0("eQTLmapping/covariates/", tissue, ".tsv"))
     GRM = import(paste0("eQTLmapping/GRMs/", tissue, ".cXX.txt"), header = FALSE)
     colnames(GRM) = rownames(GRM) = covariates$id
     phenos = import(paste0("eQTLmapping/phenotypes/", tissue, ".tsv"))
     genes = import(paste0("eQTLmapping/phenotypes/", tissue, ".genes.txt"), header = FALSE)[,1]
     
     rownames(phenos) = genes
     phenos = t(phenos)
     
     return(list(covariates = covariates,
                 GRM = GRM,
                 phenos = phenos,
                 genes = genes))
}

runNullmodel = function(current_gene, tissue, covariates, phenos, GRM){
     cache_folder = paste0('eQTLmapping/cache/V_folder/',
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
     write.table(Va, Va_file,
                 row.names = FALSE, col.names = FALSE)
     return(Va)
}

library(furrr)

plan(multisession, workers = 10)

global_formulas = list(
     head = 
y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2 + (1|id),
     body = 
y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + (1|id))

d = load_data('head')
head_va = future_map(d$genes, runNullmodel, 
                     tissue = 'head', 
                     covariates = d$covariates,
                     phenos = d$phenos,
                     GRM = d$GRM)

b = load_data('body')
body_va = future_map(b$genes, runNullmodel, 
                     tissue = 'body', 
                     covariates = b$covariates,
                     phenos = b$phenos,
                     GRM = b$GRM)

Va_df = data.frame(tissue = 'head',
                   gene = d$genes,
                   Va = unlist(head_va)) |>
     bind_rows(data.frame(tissue = 'body',
                          gene = b$genes,
                          Va = unlist(body_va)))
export(Va_df, 'eQTLmapping/Va_cache/Va.tsv')
