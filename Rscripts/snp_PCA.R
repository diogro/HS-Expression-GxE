source(here::here("Rscripts/functions.R"))
source(here::here("Rscripts/read_data.R"))

snp_pca <- list(head = snpgdsPCA(snps_head),
                body = snpgdsPCA(snps_body))