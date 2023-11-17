source(here::here("Rscripts/functions.R"))

rnaseq_data = import("cache/rnaseq_all.rds")

data = rnaseq_data$body

DE_limma = function(data){
    design <- model.matrix(~0+treatment + 
                        egglayBatch + 
                        RNAseqBatch + 
                        platingBatch + 
                        RNAlibBatch, data$covariates) |> 
              cbind(data$mwash$voom$Zhat)
    colnames(design)[1:2] <- c("ctrl","hs")
    missing = which(colnames(design) == "")
    if(length(missing) > 0){
        n = length(missing)
        colnames(design)[missing] = paste0("sv", 1:n)
    }
    contrast <- makeContrasts(HSvsC=hs-ctrl, levels = design)
    fit <- lmFit(data$voom, design)
    fitc <- contrasts.fit(fit, contrasts = contrast)
    fit2 <- eBayes(fitc)
    summary(decideTests(fit2))
    res <- topTable(fit2, number = Inf)
    significance = as.data.frame(decideTests(fit2))
    significance$Geneid = rownames(significance)
    res = left_join(res, significance)
    list(table = res, lmfit = fit, eBayes = fit2, summary = summary(decideTests(fit2, p.value = 0.01)))
}
diff_expression = future_map(rnaseq_data, DE_limma)
diff_expression |> map("summary")

table = diff_expression |> 
    map("table") |> 
    list_rbind(names_to = "tissue")
table = left_join(table,
        bitr(table$Geneid, 
             fromType="FLYBASE", 
             toType = "SYMBOL",
             OrgDb = org.Dm.eg.db, 
             drop = FALSE),
        by = c(Geneid = "FLYBASE"))
names(table)
write_csv(table, affix_date("output/HS-ctrl-DE_limma-table.csv"))
 
# from here I get the mean expression logcpm per group fit$coefficients
# and the stderror fit$stdev.unscaled * fit$sigma

