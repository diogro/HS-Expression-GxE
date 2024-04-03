pak::pkg_install("ggalluvial")
library(ggalluvial)
library(tidyverse)
library(rio)
library(patchwork)
library(cowplot)

ctrl = import("SBM/snakemake/cache/blockSummary/fdr-1e-3/head-ctrl/gene_block.csv")
hs = import("SBM/snakemake/cache/blockSummary/fdr-1e-3/head-hs/gene_block.csv")



gene_B3 =inner_join(select(ctrl, Gene, B3), 
           select(hs, Gene, B3), by = "Gene") %>% 
           as_tibble() %>%
           count(B3.x, B3.y) 
  mutate(gene = factor(Gene, levels = Gene)) %>%
  ggplot(aes(axis1 = ctrl, axis2 = hs)) +
  geom_alluvium(aes(fill = Gene)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() +
  theme(legend.position = "none")

ggplot(gene_B3, aes(y = n, axis1 = B3.x, axis2 = B3.y)) +
  geom_alluvium(width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey")
dev.off()
