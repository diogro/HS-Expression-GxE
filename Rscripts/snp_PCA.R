source(here::here("Rscripts/functions.R"))
# source(here::here("Rscripts/read_data.R"))

covariates = import("cache/covariates.rds")

snps_head = seqOpen("cache/ccm_head.gds") # seqClose("ccm_head.gds")
snps_body = seqOpen("cache/ccm_body.gds") # seqClose("ccm_body.gds")

snp_pca <- list(head = snpgdsPCA(snps_head),
                body = snpgdsPCA(snps_body))


pdf('tmp/snp_pca.pdf')
    treatment = filter(covariates, tissue == "head") |> 
        arrange(match(id, snp_pca$head$sample.id)) |>
        pull(treatment)
    plot(snp_pca$head, eig=1:5,  
        col=treatment,
        pch=20, cex=0.5)
    treatment = filter(covariates, tissue == "body") |> 
        arrange(match(id, snp_pca$body$sample.id)) |>
        pull(treatment)
    plot(snp_pca$body, eig=1:5,  
        col=treatment,
        pch=20, cex=0.5)
    par(mfrow=c(1, 2))
    pc.percent = na.omit(snp_pca$head$varprop*100)
    plot(pc.percent, pch=19)
    pc.percent = na.omit(snp_pca$body$varprop*100)
    plot(pc.percent, pch=19)
dev.off()