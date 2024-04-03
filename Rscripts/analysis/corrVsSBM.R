source(here::here("Rscripts/functions.R"))
pak::pkg_install("hexbin")

upt = function(x) x[upper.tri(x, diag = FALSE)]

DiffExp = DiffExpAll$vicar %>% filter(tissue == "head") %>% as_tibble()

DiffExp$Geneid[DiffExp$Geneid == "FBgn0085200"]

hs = import("SBM/snakemake/cache/blockSummary/fdr-1e-3/head-hs/gene_block.csv")
ctrl = import("SBM/snakemake/cache/blockSummary/fdr-1e-3/head-ctrl/gene_block.csv")

res_exp = list(ctrl = read.table(here::here("SBM/rawData/batch/head-ctrl.tsv"), sep = "\t"),
               hs = read.table(here::here("SBM/rawData/batch/head-hs.tsv"), sep = "\t"))

corrMatrices = list(ctrl = cor(t(res_exp$ctrl), method = "spearman"),
                    hs = cor(t(res_exp$hs), method = "spearman"))

long_corrMatrices = lapply(corrMatrices, function(X) {
    ind <- which(upper.tri(X, diag = TRUE), arr.ind = TRUE)
    nn <- dimnames(X)
    data.frame(row = nn[[1]][ind[, 1]],
               col = nn[[2]][ind[, 2]],
               val = X[ind]) %>% filter(val != 1.)
})
names(long_corrMatrices) = c("ctrl", "hs")

df_corrs = data.frame(long_corrMatrices$ctrl[,1:2], 
                      ctrl = long_corrMatrices$ctrl$val, 
                      hs = long_corrMatrices$hs$val) %>%
                      rename("Gene1" = "row", "Gene2" = "col") %>%
                      as_tibble()

classifyCorrBlock = function(x, blockSummary, level){   
    block_col = paste("B", level, sep = "")
    gene1_block = blockSummary[[block_col]][match(x$Gene1, blockSummary$Gene)]
    gene2_block = blockSummary[[block_col]][match(x$Gene2, blockSummary$Gene)]
    return(gene1_block == gene2_block)
}
df_corrs$B2 = classifyCorrBlock(df_corrs, ctrl, 2)

# p = ggplot(df_corrs, aes(x = ctrl, y = hs)) + 
#     geom_point(alpha = 0.1) + geom_abline(intercept = 0, slope = 1) + 
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     geom_vline(xintercept = 0, linetype = "dashed") +
#     facet_wrap(~B2) +
#     theme_classic(12) + labs(x = "Correlation (ctrl)", y = "Correlation (hs)")
# save_plot(here::here("output/corrVsSBM.png"), p, base_width = 10, base_height = 6)

subsample = sample(nrow(df_corrs), 1000000)

p = ggplot(df_corrs, aes(x = ctrl, y = hs)) + 
        geom_hex(bins = 200) +
        scale_fill_viridis() +
        geom_abline(intercept = 0, slope = 1) + 
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        facet_wrap(~B2, labeller = labeller(B2 = c("FALSE" = "Between level-2 blocks", 
                                                   "TRUE" = "Within level-2 blocks"))) +
        theme_classic(12) + labs(x = "Correlation (ctrl)", y = "Correlation (hs)")
    
save_plot(here::here("output/corrVsSBM.png"), p, base_width = 10, base_height = 6)
