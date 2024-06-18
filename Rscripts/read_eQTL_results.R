source(here::here("Rscripts/functions.R"))
#pak::pkg_install("PhenotypeSimulator")
library(PhenotypeSimulator)

rnaseq = import(here::here("cache/rnaseq_all_2024-03-21.rds"))

source(here::here("eQTLmapping/scripts/snpPosClassify.R"))
             
tissue = 'head'

Xp = read.plink(here::here(paste0("eQTLmapping/bed_files/", tissue)))
X = as(Xp$genotypes,'numeric') 
rownames(X) = Xp$fam$member
colnames(X) = Xp$map$snp.name

covariates = import(here::here(paste0("eQTLmapping/covariates/", tissue, ".tsv"))) |>
     mutate(egglayBatch = as.character(egglayBatch), 
            RNAseqBatch = as.character(RNAseqBatch),
            platingBatch = as.character(platingBatch),
            RNAlibBatch = as.character(RNAlibBatch))            
GRM = import(here::here(paste0("eQTLmapping/GRMs/", tissue, ".cXX.txt")), header = FALSE)
colnames(GRM) = rownames(GRM) = covariates$id
genes = import(here::here(paste0("eQTLmapping/phenotypes/", tissue, ".genes.txt")), header = FALSE)[,1]

# results_files = list(head = dir(here::here("eQTLmapping/results/gxe/head"), 
#                                 full.names = TRUE),
#                      body = dir(here::here("eQTLmapping/results/gxe/body"), 
#                                 full.names = TRUE)) 
# results_df = map(results_files[[tissue]], 
#                  import, select = c("Trait", "snp", "p_gxe", "cis"), 
#                  .progress = 'reading files') |> 
#                  bind_rows() |> 
#                  as_tibble()

# # countCisTests = lapply(genes, getCisSnps) |> 
# #       map_dbl(~length(.))
# # export(countCisTests, here::here(paste0("cache/countCisTests_", tissue, ".rds")))
# countCisTests = import(here::here(paste0("cache/countCisTests_", tissue, ".rds")))
# total_cis_test = sum(countCisTests)     
# total_snps = ncol(X)
# total_trans_tests = sum(total_snps - countCisTests)

# results_df$p_gxe_adj = NA
# results_df$p_gxe_adj[results_df$cis == 1] = p.adjust(results_df$p_gxe[results_df$cis == 1], method = "BH", n = total_cis_test) 
# results_df$p_gxe_adj[results_df$cis == 0] = p.adjust(results_df$p_gxe[results_df$cis == 0], method = "BH", n = total_trans_tests) 

# sig_results_df = results_df |>
#      filter(p_gxe_adj < 0.05)
# enclosing_genes_df = map(unique(sig_results_df$snp), getEnclosingGene, .progress = "Mapping SNPs to genes") |> bind_rows() 
# enclosing_genes_df = as_tibble(enclosing_genes_df)

# sig_results_df = results_df |>
#      filter(p_gxe_adj < 0.05) |>
#      left_join(enclosing_genes_df, by = "snp") |>
#      select(Trait, snp, p_gxe, p_gxe_adj, cis, index, enclosingGene) |>
#      as_tibble() 
# export(sig_results_df, here::here(paste0("output/significant_gxe_eqtl_", tissue, ".tsv")))
sig_results_df = import(here::here(paste0("output/significant_gxe_eqtl_", tissue, ".tsv")))

gemma_sig_results = import("/Genomics/ayroleslab2/lamaya/bigProject/GXEpaper/HEAD/eQTLmapping/mapping/gxe_covfree_allgenes/resultstable/resultsGEMMA.GxE.CISTRANS.fdr5.headctrlhs.feb21.2021.txt") |> mutate(chr = gsub("23", "X", chr)) |> as_tibble()

enclosing_genes_df = map(unique(gemma_sig_results$rs), getEnclosingGene, window = 0, .progress = "Mapping SNPs to genes") |> bind_rows() 
enclosing_genes_df = as_tibble(enclosing_genes_df)

enclosing_genes_df |> filter(enclosingGene == "FBgn0264908")

# Manual count os SNPs in FBgn0264908
gemma_sig_results |>
     filter(chr == "3L") |>
     filter(ps > 15848578, ps < 15899029) |>
     select(chr, rs, ps, gene, enclosingGene) |>
     group_by(rs) |>
     dplyr::count() 

gemma_sig_results = inner_join(gemma_sig_results, enclosing_genes_df, by = c("rs" = "snp")) |> 
     as_tibble() 

gemma_sig_results |> filter(chr == "4") |> select(chr, rs, id, gene, enclosingGene) |> 
     print(n = 300) 


enclosing_genes_df |>
     group_by(enclosingGene) |>
     dplyr::count() |> 
     arrange(desc(n)) |>
     filter(n>2) |>
     print(n = 300)

filter(gemma_sig_results,enclosingGene == "FBgn0028704") |>
     select(rs, p_wald, gene) |>
     print(n = 300)

inner_join(gemma_sig_results, enclosing_genes_df, by = c("rs" = "snp")) |> 
     as_tibble()  |>
     filter(enclosingGene == "FBgn0032151") |>
     select(rs, p_wald, gene) |>
     print(n = 300)

gemma_sig_results

joint_gxe_snps = inner_join(sig_results_df, gemma_sig_results, by = c("snp" = "rs")) |> 
     select(Trait, gene, snp, p_gxe, p_wald) 

png(here::here("tmp/joint_gxe_snps.png"), width = 10, height = 10, units = "in", res = 300)
ggplot(joint_gxe_snps, aes(x = -log10(p_gxe), y = -log10(p_wald))) +
     geom_point() +
     geom_abline(intercept = 0, slope = 1)
dev.off()

sig_results_df |>
     group_by(Trait, cis) |>
     dplyr::count() |>
     arrange(desc(cis)) |> 
     print(n = 300)
length(unique(sig_results_df$Trait))

sig_results_df |>
     group_by(Trait, cis) |>
     dplyr::count() |>
     arrange(desc(cis)) |> 
     filter(Trait == "FBgn0015035")

x = as_tibble(gene_locations_GR) |>
     filter(seqnames %in% paste0("chr", chrs))
x$gene_id = unlist(x$gene_id)
x = x[!duplicated(x$gene_id),]

tx.git  = GNCList(gene_locations_GR)

sig_results_df = results_df |>
     filter(p_gxe_adj < 0.05) |>
     left_join(enclosing_genes_df, by = "snp") |>
     select(Trait, snp, p_gxe, p_gxe_adj, cis, index, enclosingGene) |>
     as_tibble() 

sig_results_df = gemma_sig_results %>%
     rename(Trait = gene, 
               snp = rs, 
               p_gxe = p_wald, 
               p_gxe_adj = fdr) |>
     select(Trait, chr, snp, p_gxe, p_gxe_adj)

gemma_sig_results

# Bind all results, filter for significant qvalues at 1%
qtl_pos = sig_results_df |>
     separate(snp, c("chr", "bp"), sep = "_") |>
     select(Trait:p_gxe_adj) |>
     left_join(x, by = c("Trait" = "gene_id")) |>
     select(Trait, seqnames, start, end, chr, bp, p_gxe) |>
     mutate(seqnames = gsub("chr", "", seqnames)) |>
     mutate(chr = factor(chr, levels = chrs),
            seqnames = factor(seqnames, levels = chrs),
            bp = as.numeric(bp), 
            log10p = -log10(p_gxe)) |>
     filter(!is.na(seqnames))

 data_cum_snp <- qtl_pos %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(bp)) %>% 
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(chr, bp_add)

 data_cum <- qtl_pos %>% 
    group_by(seqnames) %>% 
    summarise(max_start = max(start)) %>% 
    mutate(start_add = lag(cumsum(max_start), default = 0)) %>% 
    select(seqnames, start_add)

data_cum$start_add = map2_dbl(data_cum$start_add, data_cum_snp$bp_add, max)
data_cum_snp$bp_add = map2_dbl(data_cum$start_add, data_cum_snp$bp_add, max)

qtl_pos <- qtl_pos %>% 
    inner_join(data_cum_snp, by = "chr") %>% 
    mutate(bp_cum = bp + bp_add)

# count row per chr
qtl_pos |>
     group_by(seqnames) |>
     dplyr::count()
# count row per chr
qtl_pos |>
     group_by(chr) |>
     dplyr::count()

qtl_pos <- qtl_pos %>% 
    inner_join(data_cum, by = "seqnames") %>% 
    mutate(start_cum = start + start_add) |>
    dplyr::select(-bp_add, -start_add)

eqtl_GR <- GRanges(
    seqnames = paste0("chr", qtl_pos$chr),
    ranges = IRanges(qtl_pos$bp, end = qtl_pos$bp, 
                     names = paste(qtl_pos$chr, qtl_pos$bp, sep = "_"), 
                     trait = qtl_pos$Trait)
                     )
overlap_git = findOverlaps(eqtl_GR, tx.git)
qtl_gene_overlap = eqtl_GR[overlap_git@from,]
qtl_gene_overlap$gene_id = tx.git[overlap_git@to]$gene_id
qtl_gene_overlap = unique(qtl_gene_overlap)

qtl_g_over_df = as_tibble(qtl_gene_overlap) |>
     mutate(seqnames = gsub("chr", "", seqnames), 
            gene_id = unlist(gene_id))

filter(qtl_g_over_df, gene_id == "FBgn0264908")
filter(qtl_g_over_df, gene_id == "FBgn0000114") |> select(trait, gene_id) |> unique()

hotspot_count = as_tibble(qtl_gene_overlap) |>
     mutate(seqnames = gsub("chr", "", seqnames), 
            gene_id = unlist(gene_id)) |>
     group_by(seqnames, trait, gene_id) |>  
     dplyr::count() |>
     arrange(desc(n)) |>
     ungroup() |>
     group_by(gene_id) |>
     dplyr::count() |>
     arrange(desc(n)) |>
     filter(n > 2)

gene_pos = as_tibble(gene_locations_GR) |>
     filter(seqnames %in% paste0("chr", chrs))
gene_pos$gene_id = unlist(x$gene_id)
gene_pos = gene_pos[!duplicated(gene_pos$gene_id),]
gene_pos$gene_id = unlist(gene_pos$gene_id)

hotspot = gene_pos |> 
     filter(gene_id %in% unlist(hotspot_count$gene_id)) |> 
     mutate(seqnames = gsub("chr", "", seqnames)) 
hotspot$cum_start = hotspot$start +  data_cum_snp[match(hotspot$seqnames, data_cum_snp$chr),]$bp_add
hotspot$cum_end = hotspot$end +  data_cum_snp[match(hotspot$seqnames, data_cum_snp$chr),]$bp_add
hotspot = unique(hotspot)

png(here::here("tmp/eQTL_GxE_pos.png"), width = 10, height = 10, units = "in", res = 300)
ggplot(qtl_pos, aes(bp_cum, start_cum)) +
    geom_point(size = 0.5) +
    geom_vline(data = data_cum, aes(xintercept = start_add), linetype = "dashed", linewidth = 0.2) +
    geom_hline(data = data_cum, aes(yintercept = start_add), linetype = "dashed", linewidth = 0.2) +
    #geom_vline(data = hotspot, aes(xintercept= cum_start), linetype = "dashed", linewidth = 0.2, color = 2) +
    #geom_vline(data = hotspot, aes(xintercept= cum_end), linetype = "dashed", linewidth = 0.2, color = 2) +
     #annotate("text", x = hs_start + 10e5, y = 1000 , label = "pHCl-1", vjust = 1, hjust = 0) +
    theme_classic() +
    labs(x = "eQTL position", y = "Gene position") +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.2, alpha = 0.5) 
dev.off()


filter(qtl_pos, bp_cum > hs_start, bp_cum < hs_end) |>
     arrange(log10p) |>
     group_by(Trait) |>
     dplyr::count()

# SNP validation plot

rnaseq[[tissue]] |> names()

data_set = "l2c"

gemma_sig_results_df = gemma_sig_results %>%
     rename(Trait = gene, 
               snp = rs, 
               p_gxe = p_wald, 
               p_gxe_adj = fdr) |>
     select(Trait, chr, snp, p_gxe, p_gxe_adj)

top_snps = sig_results_df %>% 
  mutate(chr = snp) |>
  separate(chr, c("chr", "pos"), sep = "_") |>
  group_by(Trait, chr, cis) %>%
  slice(which.min(p_gxe)) |>
  arrange(p_gxe) |>
  select(Trait, chr, pos, snp, p_gxe, cis, enclosingGene) |>
  print(n = 300)

current_snp = top_snps$snp[1]
plot_GxE_eQTL(current_snp)
walk(top_snps$snp, plot_GxE_eQTL)


## Gemma grid comparison:

grid = runGxEmodel('FBgn0004396', tissue, covariates, GRM)
gemma_results = import(here::here("eQTLmapping/cache/gemma/head/FBgn0004396/output/gemma.assoc.txt"))

head(gemma_results)
head(grid)
gemma_gridlmm = inner_join(grid, gemma_results, by = c("snp" = "rs")) |>
     select(snp, p_score, p_gxe) |>
     arrange(p_gxe) |> as_tibble()
with(gemma_gridlmm, cor(p_score, p_gxe))
p = ggplot(gemma_gridlmm, aes(x = p_score, y = p_gxe)) +
     geom_point() +
     geom_abline(intercept = 0, slope = 1) +
     scale_x_log10() +
     scale_y_log10()
save_plot(here::here("tmp/gemma_gridlmm.png"), p, base_width = 7)

# # Run GxE model GEMMA
# runGxEmodelGEMMA = function(current_gene, tissue, covariates, GRM){
#      old_dir = getwd()
#      cache_folder = here::here(paste0('eQTLmapping/cache/gemma/'),
#                            tissue, '/',
#                            current_gene)
#      #create cache folder
#      if(!dir.exists(cache_folder)){
#         dir.create(cache_folder, recursive = TRUE)
#      }
#      y = t(import(here::here(paste0("eQTLmapping/phenotypes/", tissue, ".tsv")), 
#                   skip = which(genes == current_gene), nrows = 1))
#      data = covariates |>
#           mutate(y = y) |>
#           as.data.frame()
#      mod1 <- model.matrix(global_formulas_gemma[[tissue]], covariates) 
#      write.table(mod1, paste0(cache_folder, "/mod1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
#      write.table(data$y, paste0(cache_folder, "/pheno.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
#      write.table(data$treatment, paste0(cache_folder, "/env.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
#      setwd(cache_folder)
#      system("gemma -bfile ../../../../bed_files/head -p pheno.txt -c mod1.txt -gxe env.txt -d /Genomics/argo/users/damelo/projects/HS-Expression-GxE/eQTLmapping/GRMs/head.eigenD.txt -u /Genomics/argo/users/damelo/projects/HS-Expression-GxE/eQTLmapping/GRMs/head.eigenU.txt -lmm 4 -o gemma")
#      gemma_results = import("output/gemma.assoc.txt") |>
#           rename(snp = "rs",
#                  p_gxe = "p_lrt") |>
#           select(snp, p_gxe) |>
#           mutate(Trait = current_gene) |>
#           as_tibble()
#      setwd(old_dir)
#      results = gxe_gwas$results
#      out_file = results |>
#           mutate(Trait = current_gene) |>
#           rename(snp = X_ID,
#                  p_main = p_value_REML.1,
#                  p_gxe = p_value_REML.2) |>
#           select(Trait, snp,  p_main, p_gxe) |> 
#           as_tibble()
#      qvalues = qvalue(out_file$p_gxe, fdr.level = 0.01)
#      out_file = out_file |> 
#           mutate(q_gxe = qvalues$qvalues) 
#      return(out_file)
# }

library(GenomicRanges)
df1 = data.frame(chr = c("chr1",  "chr12"), start = c(10000,  10000), end = c(20000, 20000))
df2 = data.frame(chr = rep("chr1", 4), posn = c(100, 12000, 15000, 250000), x = rep(1, 4), y = rep(2, 4), z = rep(3, 4))

df1.ir = IRanges(start = df1$start, end = df1$end, names = df1$chr)
df2.ir = IRanges(start = df2$posn, end = df2$posn, names = df2$chr) 

df1.it = NCList(df1.ir)

overlap_it = findOverlaps(df2.ir, df1.it)
print(overlap_it)

# define a genomic interval tree (git) for faster search 

df1.gr = GRanges (IRanges(start = df1$start, end = df1$end), seqnames=df1$chr) 
df2.gr = GRanges(IRanges(start=df2$posn, end = df2$posn), seqnames = df2$chr) 

df1.git  = GNCList(df1.gr)
t1 = Sys.time ()
overlap_git = findOverlaps(df2.gr, df1.git)
t2 = Sys.time()
sub1 = difftime(t2, t1, tz, units = c("auto"))
print (sub1)
print (overlap_git)
