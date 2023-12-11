library(rio)
source("Rscripts/functions.R")

rnaseq_data = import("cache/rnaseq_all_2023-11-17.rds")

treatment = pull(rnaseq_data$head$covariates, treatment)
head_ctrl = rnaseq_data$head$mwash_residuals$l2c[,treatment == 0]
head_hs   = rnaseq_data$head$mwash_residuals$l2c[,treatment == 1]

treatment = pull(rnaseq_data$body$covariates, treatment)
body_ctrl = rnaseq_data$body$mwash_residuals$l2c[,treatment == 0]
body_hs   = rnaseq_data$body$mwash_residuals$l2c[,treatment == 1]

write.table(head_ctrl, here::here("SBM/rawData/batch/head-ctrl.tsv"), sep = "\t")
write.table(head_hs, here::here("SBM/rawData/batch/head-hs.tsv"), sep = "\t")

write.table(body_ctrl, here::here("SBM/rawData/batch/body-ctrl.tsv"), sep = "\t")
write.table(body_hs, here::here("SBM/rawData/batch/body-hs.tsv"), sep = "\t")
