source(here::here("Rscripts/functions.R")) 

DiffExpAll = list(limma = import(here::here("output/HS-ctrl-DE_limma-table_2024-03-21.csv")),
                  vicar = import(here::here("output/HS-ctrl-DE_vicar-table_2024-03-21.csv")),
                  DESeq = import(here::here("output/HS-ctrl-DE_DESeq2-table_2024-03-21.csv"))) 

DiffExp = DiffExpAll$vicar %>% filter(tissue == "head") %>% as_tibble()

head = import(here::here("SBM/snakemake/cache/blockSummary/fdr-1e-3/layered/head/gene_block.csv"))
block_summary_clip = import(here::here("cache/clip/blockSummary_clip_fdr-1e2_head.csv"))

block_summary = import(here::here("SBM/snakemake/cache/GO/fdr-1e-3/layered/head/blockSummary.csv")) 

block_summary_clip = import(here::here("cache/clip/blockSummary_clip_fdr-1e2_head.csv"))
block_summary_clip = block_summary_clip |> 
    filter(Nested_Level == 1, Direction == "decohere") |>
    arrange(Nested_Level, desc(Total_degree))
head(block_summary_clip, 20)

DiffSBM = inner_join(DiffExp, select(head, Gene, B2), by = join_by(Geneid == Gene)) %>% 
    inner_join(block_summary, by = join_by(B2 == Block))

block_median = inner_join(DiffExp, select(head, Gene, B2), by = join_by(Geneid == Gene)) %>%   
   group_by(B2) %>% dplyr::summarize(median = median(-log(lfdr))) %>% arrange(desc(median))
block_median$B2 = factor(block_median$B2, levels = block_median$B2)   
DiffSBM$B2 = factor(DiffSBM$B2, levels = block_median$B2)
DE_blocks = block_median[block_median$median > 5,]

clip_corrs = import(here::here("cache/long_clip.csv"))
all_corrs = import(here::here("cache/long_corrs.csv")) 

sortGenes <- function(x){
    Gene1 = with(x, pmin(Gene1, Gene2))
    Gene2 = with(x, pmax(Gene1, Gene2))
    x$Gene1 = Gene1
    x$Gene2 = Gene2
    return(x)
}
clip_corrs = sortGenes(clip_corrs)
all_corrs = sortGenes(all_corrs)

clip_corrs <- left_join(clip_corrs, all_corrs, by = c("Gene1", "Gene2")) 

DiffSBM_DE_blocks = DiffSBM %>% filter(B2 %in% DE_blocks$B2)
DE_genes = DiffSBM_DE_blocks$Geneid
gene1_DE = clip_corrs$Gene1 %in% DE_genes
gene2_DE = clip_corrs$Gene2 %in% DE_genes
clip_corrs$DE = gene1_DE & gene2_DE
DE_corrs = filter(clip_corrs, DE) |> select(-DE)




classifyCorrBlock = function(x, blockSummary, level){   
    block_col = paste("B", level, sep = "")
    gene1_block = blockSummary[[block_col]][match(x$Gene1, blockSummary$Gene)]
    gene2_block = blockSummary[[block_col]][match(x$Gene2, blockSummary$Gene)]
    return(gene1_block == gene2_block)
}
DE_corrs$B2 = classifyCorrBlock(DE_corrs, head, 2)
p = pivot_longer(DE_corrs, ctrl:hs, names_to = "treatment", values_to =  "correlation") |>
    ggplot(aes(B2, correlation, fill = treatment, group = interaction(treatment, B2))) + 
    geom_boxplot() 
save_plot(here::here("output/corrTreatmentSBM.png"), p, base_width = 10, base_height = 6)

p = DiffSBM %>%
    ggplot(aes(B2, -log(lfdr), group = B2)) + 
    geom_boxplot(aes(fill = is_enriched, color = is_enriched)) +
    geom_point(data = block_median, aes(B2, median)) + 
    geom_line(data = block_median, aes(B2, median)) + 
    theme_classic(12) + labs(x = "SBM Block", y = "Differential Expression (-log(lfdr))") + theme(axis.text.x = element_text(angle = 45), legend.position = c(0.8, 0.9)) + 
    geom_hline(yintercept = 5, linetype = "dashed", color = 2)

save_plot(here::here("output/DEvsSBM.png"), p, base_width = 10, base_height = 6)


p = ggplot() + 
    geom_point(data = block_median, aes(B2, median)) + 
    geom_line(data = block_median, aes(B2, median)) + 
    theme_classic() + labs(x = "SBM Block", y = "median -log(lfdr)") 
save_plot(here::here("output/DEvsSBM_line.png"), p)
