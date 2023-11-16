source(here::here("Rscripts/functions.R"))

rnaseq_data = import("cache/rnaseq_all.rds")

data = rnaseq_data$head

library("BiocParallel")
register(MulticoreParam(4))

DE_DESeq2 <- function(data){
    colnames(data$mod1)[1:2] <- c("ctrl","hs")
    rownames(data$counts$counts) <- data$counts$genes[,1]
    dds <- DESeqDataSetFromMatrix(countData = data$counts$counts,
                                colData = data$covariates,
                                design = data$mod1)
    dds <- estimateSizeFactors(dds)

    dat <- vst(dds, blind = FALSE)
    dat_noBatch <- removeBatchEffect(assay(dat), 
                                design = data$design, 
                                covariates = data$mod0)
    norm_counts_pca = pca(dat_noBatch)
    col = pull(data$covariates, treatment)
    png(paste0('tmp/vstCounts-PCA-', data$tissue, '-forDESeq2.png'), width = 1080, height = 1080)
        print(pca_plot(norm_counts_pca, c("C", "HS")[col]))
    dev.off()

    dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(4))
    res <- results(dds, contrast=list("hs","ctrl"))
    resultsNames(dds)

    list(table = res, DESeq = dds)
}
diff_expression = map(rnaseq_data, DE_DESeq2)

table = diff_expression |> 
    map("table") |> 
    map(as.data.frame) %>%
    map(\(x) mutate(x, Geneid = rownames(x))) |>
    map(as_tibble) |>
    list_rbind(names_to = "tissue") |>
    relocate(tissue, Geneid)
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
