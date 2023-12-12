library(GridLMM)
library(rio)

# Load data
covariates = import("eQTLmapping/covariates/head.tsv")
GRM = import("eQTLmapping/GRMs/head.cXX.txt", header = FALSE)
colnames(GRM) = rownames(GRM) = covariates$id
phenos = import("eQTLmapping/phenotypes/head.tsv")
genes = import("eQTLmapping/phenotypes/head.genes.txt", header = FALSE)[,1]

rownames(phenos) = genes
phenos = t(phenos)

data = covariates |>
     mutate(y = phenos[,2]) |>
     filter(id %in% rownames(GRM)) |>
     as.data.frame()
rownames(data) = data$id
global_formula = y ~ 1 + egglayBatch + RNAseqBatch +  platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2 + (1|id)
null_model = GridLMM_ML(formula = global_formula,
                        data = data,
                        relmat = list(id = list(K = GRM)),
                        REML = T,
                        save_V_folder = 'V_folder',
                        tolerance = 1e-3)

null_model$results[,c('id.REML')]
