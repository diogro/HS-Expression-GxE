library(rio)
library(tidyverse)

# rna_seq_file = "cache/rnaseq_all_2023-11-17.rds"
rna_seq_file = snakemake@input[["rna_seq"]]
all_data = import(rna_seq_file)

sample_list_file = snakemake@input[["sample_list"]]

tissue = snakemake@wildcards[["tissue"]]
data = all_data[[tissue]]

samples = read.table(sample_list_file, header = FALSE, sep = "\t")[,1]

colnames(data$mwash$l2c$Zhat) = paste0("Zhat", 1:dim(data$mwash$l2c$Zhat)[2])

input_df = cbind(data$covariates, 
                 data$mwash$l2c$Zhat) |>
               filter(id %in% samples)
input_df = input_df[match(samples, input_df$id),]

stopifnot(all(input_df$id == samples))
phenos = data$l2c[,samples]
write.table(rownames(phenos), snakemake@output[["gene_list"]], 
          row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
export(input_df, snakemake@output[["covariates"]])
export(phenos, snakemake@output[["phenotypes"]])

