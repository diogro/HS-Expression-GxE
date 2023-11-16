source(here::here("Rscripts/functions.R"))

rnaseq_data = import("cache/rnaseq_all.rds")

data = rnaseq_data$body

DE_limma = function(data){
    design <- cbind(data$mod1, data$sva$voom$sv)
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
diff_expression = map(rnaseq_data, DE_limma)
diff_expression |> map("summary")

table = diff_expression |> 
    map("table") |> 
    list_rbind(names_to = "tissue")
table = left_join(table
        bitr(table$Geneid, 
             fromType="FLYBASE", 
             toType = "SYMBOL",
             OrgDb = org.Dm.eg.db, 
             drop = FALSE),
        by = c(Geneid = "FLYBASE"))
names(table)
write_csv(table, affix_date("output/HS-ctrl-DE-table.csv"))
 
# from here I get the mean expression logcpm per group fit$coefficients
# and the stderror fit$stdev.unscaled * fit$sigma

str(diff_expression$head$table)
names(diff_expression$head$table)
{ggplot(table, 
        aes(logFC, -log10(adj.P.Val))) +
    geom_point(size = 0.1, alpha = 0.3, aes(color = as.factor(HSvsC))) + 
    geom_text_repel(data = table |> filter(-log10(adj.P.Val) > 50), 
                    aes(label = SYMBOL), 
                    size = 2, max.overlaps = 30, segment.size = 0.1) + 
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~tissue) + theme_bw() + 
    theme(legend.position='none') } %>%
{save_plot("tmp/DE_volcano.png", ., base_width = 7)}
