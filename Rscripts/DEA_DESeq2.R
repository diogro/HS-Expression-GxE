source(here::here("Rscripts/functions.R"))

rnaseq_data = import("cache/rnaseq_all.rds")

data = rnaseq_data$head

library("BiocParallel")
register(MulticoreParam(4))

DE_DESeq2 <- function(data){
    design <- data$mod1
    design[,1] = abs(design[,2]-1)
    colnames(design)[1:2] <- c("ctrl","hs")
    rownames(data$counts$counts) <- data$counts$genes[,1]
    dds <- DESeqDataSetFromMatrix(countData = data$counts$counts,
                                colData = data$covariates,
                                design = design)
    dds <- estimateSizeFactors(dds)
    mod0 = data$mod0
    if(data$tissue == "head"){
        normalized_counts <- counts(dds, normalized=TRUE)
        mwash <- mouthwash(Y = t(normalized_counts), X = data$mod1, k = 2, cov_of_interest = 2,
                   include_intercept = FALSE)
        design(dds) <- cbind(design, mwash$Zhat)
        mod0 =  cbind(mod0, mwash$Zhat)
    }

    dat <- vst(dds, blind = FALSE)
    dat_noBatch <- removeBatchEffect(assay(dat), 
                                design = design, 
                                covariates = mod0)
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
export(diff_expression, affix_date("cache/DE_DESeq2.rds"))
diff_expression = import("cache/DE_DESeq2_2023-11-17.rds")

table = diff_expression |> 
    map("table") |> 
    map(as.data.frame) %>%
    map(\(x) mutate(x, Geneid = rownames(x))) |>
    map(as_tibble) |>
    list_rbind(names_to = "tissue") |>
    relocate(tissue, Geneid)
table = left_join(table,
        bitr(table$Geneid, 
             fromType="FLYBASE", 
             toType = "SYMBOL",
             OrgDb = org.Dm.eg.db, 
             drop = FALSE),
        by = c(Geneid = "FLYBASE"))

names(table)
write_csv(table, affix_date("output/HS-ctrl-DE_DESeq2-table.csv"))
 

table = table |> 
    dlply("tissue") |> 
    map(\(x) mutate(x, qvalue = qvalue(pvalue)$qvalue)) |> 
    list_rbind()

label_head = table |> 
    filter(tissue == "head",
           -log10(qvalue) > 11)
label_body = table |> 
    filter(tissue == "body",
           -log10(qvalue) > 15)
labels = bind_rows(label_head, label_body )
names(table)
{p2 = ggplot(table, 
        aes(log2FoldChange, -log10(qvalue))) +
    geom_point(size = 0.1, alpha = 0.1, color = "gray") + 
    geom_point(data = table |> filter(padj < 0.05), size = 0.1, alpha = 0.5, color = "blue") + 
    geom_text_repel(data = labels, 
                    aes(label = SYMBOL), 
                    size = 2, max.overlaps = 30, segment.size = 0.1) + 
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~tissue) + theme_bw() + 
    theme(legend.position='none') }
save_plot("tmp/DE_volcano_DESeq2.png", p2, base_width = 7)
