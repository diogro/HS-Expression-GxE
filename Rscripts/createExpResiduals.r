library(rio)
source("Rscripts/functions.R")

rnaseq_data = import("cache/rnaseq_all.rds")

treatment = pull(rnaseq_data$head$covariates, treatment)
head_ctrl = rnaseq_data$head$cpm_residuals[,treatment == 1]
head_hs = rnaseq_data$head$cpm_residuals[,treatment == 2]

treatment = pull(rnaseq_data$body$covariates, treatment)
body_ctrl = rnaseq_data$body$cpm_residuals[,treatment == 1]
body_hs = rnaseq_data$body$cpm_residuals[,treatment == 2]

write.table(head_ctrl, here::here("SBM/rawData/batch/head-ctrl.tsv"), sep = "\t")
write.table(head_hs, here::here("SBM/rawData/batch/head-hs.tsv"), sep = "\t")

write.table(body_ctrl, here::here("SBM/rawData/batch/body-ctrl.tsv"), sep = "\t")
write.table(body_hs, here::here("SBM/rawData/batch/body-hs.tsv"), sep = "\t")



