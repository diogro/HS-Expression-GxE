library(rio)
source("Rscripts/functions.R")

rnaseq_data = import("cache/rnaseq_all_2024-03-21.rds", trust = TRUE)

samples = import("eQTLmapping/sample_list/head.txt", header = F)[,1]

covariates = import("eQTLmapping/covariates/head.tsv") 
stopifnot(covariates$id == samples)

treatment = pull(covariates, treatment)
expr_df = t(rnaseq_data$head$mwash_residuals$l2c[,samples])

ID_C = covariates[covariates$treatment == 0,"id"]
ID_HS = covariates[covariates$treatment == 1,"id"]

Z_HS = as.numeric(rownames(expr_df) %in% ID_HS)
Z_C = as.numeric(rownames(expr_df) %in% ID_C)

GRM = as.matrix(read.table("eQTLmapping/GRMs/head.cXX.txt"))
rownames(GRM) = colnames(GRM) = samples

x = expr_df[ID_C,]
out_df= data.frame(fam = NA, ID = rownames(x), x)
write.table(out_df, 
            "gcta/log2Counts_head_ctrl_covfree.txt", 
            quote = FALSE,
            col.names = FALSE, row.names = FALSE)

x = expr_df[ID_HS,]
out_df= data.frame(fam = NA, ID = rownames(x), x)
write.table(out_df, 
            "gcta/log2Counts_head_hs_covfree.txt", 
            quote = FALSE,
            col.names = FALSE, row.names = FALSE)

x = expr_df
out_df= data.frame(fam = NA, ID = rownames(x), x)
gxe = data.frame(fam = NA, ID = rownames(x), c("ctrl", "HS")[covariates$treatment+1])
write.table(out_df, 
            "gcta/log2Counts_head_hsctrl_covfree.txt",
            quote = FALSE,
            col.names = FALSE, row.names = FALSE)
write.table(gxe, 
            "gcta/HSC.gxe",
            quote = FALSE,
            col.names = FALSE, row.names = FALSE)

write.table(data.frame(colnames(expr_df)), file = "gcta/gene_list.csv", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

if(!require(genio)){install.packages("genio"); library(genio)}
write_grm("gcta/grm", kinship = GRM, 
          fam = data.frame(id = rownames(expr_df), fam = NA))
