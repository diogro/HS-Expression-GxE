library(rio)
source("Rscripts/functions.R")

rnaseq_data = import("cache/rnaseq_all_2025-08-05.rds")

treatment = pull(rnaseq_data$head$covariates, treatment)
head_ctrl = rnaseq_data$head$mwash_residuals$l2c[,treatment == 0]
head_hs   = rnaseq_data$head$mwash_residuals$l2c[,treatment == 1]

treatment = pull(rnaseq_data$body$covariates, treatment)
body_ctrl = rnaseq_data$body$mwash_residuals$l2c[,treatment == 0]
body_hs   = rnaseq_data$body$mwash_residuals$l2c[,treatment == 1]

# if necessary, create directory for output files
if(!dir.exists(here::here("SBM/rawData/batch"))){
    dir.create(here::here("SBM/rawData/batch"), recursive = TRUE)
}
if(!dir.exists(here::here("SBM/rawData/layered/body/"))){
    dir.create(here::here("SBM/rawData/layered/body/"), recursive = TRUE)
}   
if(!dir.exists(here::here("SBM/rawData/layered/head/"))){
    dir.create(here::here("SBM/rawData/layered/head/"), recursive = TRUE)
}
write.table(head_ctrl, here::here("SBM/rawData/batch/head-ctrl.tsv"), sep = "\t")
write.table(head_hs, here::here("SBM/rawData/batch/head-hs.tsv"), sep = "\t")

write.table(body_ctrl, here::here("SBM/rawData/batch/body-ctrl.tsv"), sep = "\t")
write.table(body_hs, here::here("SBM/rawData/batch/body-hs.tsv"), sep = "\t")

write.table(head_ctrl, here::here("SBM/rawData/layered/head/ctrl.tsv"), sep = "\t")
write.table(head_hs, here::here("SBM/rawData/layered/head/hs.tsv"), sep = "\t")

write.table(body_ctrl, here::here("SBM/rawData/layered/body/ctrl.tsv"), sep = "\t")
write.table(body_hs, here::here("SBM/rawData/layered/body/hs.tsv"), sep = "\t")

