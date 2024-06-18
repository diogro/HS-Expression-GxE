source(here::here("Rscripts/functions.R"))
#pak::pkg_install("PhenotypeSimulator")
library(PhenotypeSimulator)

rnaseq = import(here::here("cache/rnaseq_all_2024-03-21.rds"))

source(here::here("eQTLmapping/scripts/snpPosClassify.R"))
             
tissue = 'head'

Xp = read.plink(here::here(paste0("eQTLmapping/bed_files/", tissue)))
X = as(Xp$genotypes,'numeric') 
rownames(X) = Xp$fam$member
colnames(X) = Xp$map$snp.name

covariates = import(here::here(paste0("eQTLmapping/covariates/", tissue, ".tsv"))) |>
     mutate(egglayBatch = as.character(egglayBatch), 
            RNAseqBatch = as.character(RNAseqBatch),
            platingBatch = as.character(platingBatch),
            RNAlibBatch = as.character(RNAlibBatch))            
GRM = import(here::here(paste0("eQTLmapping/GRMs/", tissue, ".cXX.txt")), header = FALSE)
colnames(GRM) = rownames(GRM) = covariates$id
genes = import(here::here(paste0("eQTLmapping/phenotypes/", tissue, ".genes.txt")), header = FALSE)[,1]

# current_gene = gxe_genes[5]
# current_gene = "FBgn0000579"

# results = vector("list", length = length(gxe_genes))
# k = 1

# for(i in k:length(gxe_genes)){
#     print(paste0("Current gene: ", gxe_genes[i], " Percent: ", 100*round(i/length(gxe_genes), 3), "%"))
#     results[[k]] = runGxEmodel(gxe_genes[i], tissue, covariates, GRM)
#     k = k + 1
# }
# export(results, here::here(paste0("cache/eQTL_detections_gxe-", tissue, ".rds")))
results = import(here::here(paste0("cache/eQTL_detections_gxe-", tissue, ".rds"))) 

x = results[[1]]
any(x$q_gxe < 0.1)


signCut <- 0.1

# Make a manhattan plot
res = results[[1]]
res |> arrange(q_gxe)
gene = gxe_genes[1]

plotManhattan(res, gene, signCut, tissue)
# map2(results, gxe_genes, plotManhattan, signCut, tissue, .progress = "Manhattan Plots")
