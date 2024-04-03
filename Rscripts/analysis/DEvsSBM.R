source(here::here("Rscripts/functions.R"))

DiffExpAll = list(limma = import(here::here("output/HS-ctrl-DE_limma-table_2023-11-17.csv")),
                  vicar = import(here::here("output/HS-ctrl-DE_vicar-table_2023-11-17.csv")),
                  DESeq = import(here::here("output/HS-ctrl-DE_DESeq2-table_2023-11-18.csv"))) 

DiffExp = DiffExpAll$vicar %>% filter(tissue == "head") %>% as_tibble()

DiffExp$Geneid[DiffExp$Geneid == "FBgn0085200"]

hs = import("SBM/snakemake/cache/blockSummary/fdr-1e-3/head-hs/gene_block.csv")
ctrl = import("SBM/snakemake/cache/blockSummary/fdr-1e-3/head-ctrl/gene_block.csv")

block_summary = import("SBM/snakemake/cache/GO/fdr-1e-3/head-ctrl/block_summary.csv") %>% filter(Nested_Level == 2) %>% mutate(is_enriched = n_enrich != 0)

DiffSBM = inner_join(DiffExp, select(ctrl, Gene, B2), by = join_by(Geneid == Gene)) %>% 
    inner_join(block_summary, by = join_by(B2 == Block))

block_median = inner_join(DiffExp, select(ctrl, Gene, B2), by = join_by(Geneid == Gene)) %>%   
   group_by(B2) %>% dplyr::summarize(median = median(-log(lfdr))) %>% arrange(desc(median))
block_median$B2 = factor(block_median$B2, levels = block_median$B2)   
DiffSBM$B2 = factor(DiffSBM$B2, levels = block_median$B2)


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
