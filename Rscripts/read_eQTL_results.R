source(here::here("Rscripts/functions.R"))
#pak::pkg_install("PhenotypeSimulator")
library(PhenotypeSimulator)

rnaseq = import(here::here("cache/rnaseq_all_2024-03-21.rds"))

# pak::pkg_install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
tx <- transcripts(txdb, columns = "gene_id")

chrs <- c("2L", "2R", "3L", "3R", "4", "X")

detection_files = list(head = dir(here::here("eQTLmapping/detections/gxe/head")),
                          body = dir(here::here("eQTLmapping/detections/gxe/body")))              
gxe_genes = lapply(detection_files, function(x) gsub(".tsv", "", x)) 

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
gxe_genes = gxe_genes[[tissue]]

# current_gene = gxe_genes[5]
# current_gene = "FBgn0000579"

global_formulas <- list(
          head = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2 + (1|id),
          body = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 +         (1|id)
)
global_formulas_gemma <- list(
          head = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2,
          body = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat
)

runGxEmodel = function(current_gene, tissue, covariates, GRM){
     old_dir = getwd()
     setwd(here::here("eQTLmapping"))
     cache_folder = here::here(paste0('eQTLmapping/cache/'),
                           tissue, '/',
                           current_gene)
     y = t(import(here::here(paste0("eQTLmapping/phenotypes/", tissue, ".tsv")), 
                  skip = which(genes == current_gene), nrows = 1))
     data = covariates |>
          mutate(y = y) |>
          as.data.frame()
     rownames(data) = data$id
     V_setup = import(paste0(cache_folder, '/V_setup.rds'))
     gxe_gwas = GridLMM_GWAS(formula = global_formulas[[tissue]],
                             test_formula =  ~1 + treatment,
                             reduced_formula = ~1,
                             data = data,
                             X = X,
                             X_ID = 'id',
                             relmat = list(id = list(K = GRM)),
                             V_setup = V_setup,
                             method = 'REML',
                             mc.cores = 1,
                             verbose = F)
     setwd(old_dir)
     results = gxe_gwas$results
     out_file = results |>
          mutate(Trait = current_gene) |>
          rename(snp = X_ID,
                 p_main = p_value_REML.1,
                 p_gxe = p_value_REML.2) |>
          select(Trait, snp,  p_main, p_gxe) |> 
          as_tibble()
     qvalues = qvalue(out_file$p_gxe, fdr.level = 0.1)
     out_file = out_file |> 
          mutate(q_gxe = qvalues$qvalues) 
     return(out_file)
}

# results = vector("list", length = length(gxe_genes))
# k = 1

# for(i in k:length(gxe_genes)){
#     print(paste0("Current gene: ", gxe_genes[i], " Percent: ", 100*round(i/length(gxe_genes), 3), "%"))
#     results[[k]] = runGxEmodel(gxe_genes[i], tissue, covariates, GRM)
#     k = k + 1
# }
# export(results, here::here(paste0("cache/eQTL_detections_gxe-", tissue, ".rds")))
results = import(here::here(paste0("cache/eQTL_detections_gxe-", tissue, ".rds"))) 

x = results[[1]]
any(x$q_gxe < 0.1)


signCut <- 0.1

# Make a manhattan plot
res = results[[1]]
res |> arrange(q_gxe)
gene = gxe_genes[1]
plotManhattan <- function(res, gene, signCut, tissue){
    qvalues = qvalue(res$p_gxe, fdr.level = 0.1)
    res = res |> 
          mutate(q_gxe = qvalues$qvalues) |>
          arrange(q_gxe)
    res = res |>
        mutate(log_p_gxe = -log10(p_gxe)) |>
        separate(snp, c("chr", "bp"), sep = "_") |>
        mutate(chr = factor(chr, levels = chrs),
                bp = as.numeric(bp)) |>
        arrange(q_gxe) 

    gwas_data_load = res
    sig_data <- gwas_data_load %>% 
    subset(q_gxe < signCut)
    notsig_data <- gwas_data_load %>% 
    subset(q_gxe >= signCut) %>%
    group_by(chr) %>% 
    sample_frac(0.5)
    gwas_data <- bind_rows(sig_data, notsig_data)

    data_cum <- gwas_data %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(bp)) %>% 
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(chr, bp_add)

    gwas_data <- gwas_data %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(bp_cum = bp + bp_add)

    axis_set <- gwas_data %>% 
    filter(chr != 4) %>%
    group_by(chr) %>% 
    summarize(center = mean(bp_cum))

    ylim <- gwas_data %>% 
    filter(p_gxe == min(p_gxe)) %>% 
    mutate(ylim = abs(floor(log10(p_gxe))) + 2) %>% 
    pull(ylim)

    manhplot <- ggplot(filter(gwas_data, q_gxe < signCut, chr != 4), aes(x = bp_cum, y = -log10(p_gxe), 
                                    color = chr)) +
    geom_hline(yintercept = -log10(signCut), color = "grey40", linetype = "dashed") + 
    geom_point(alpha = 0.7, size = 0.005) +
    geom_point(data = filter(gwas_data, q_gxe > signCut), alpha = 0.75, size = 0.005, color = "grey") +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim), breaks = seq(0, 10, 1)) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosome", 
        y = expression("-log10(" ~ p[GxE] ~ ")")) + 
    theme_classic() +
    theme( 
        legend.position = "none",
        #panel.grid.major.x = element_blank(),
        #panel.grid.minor.x = element_blank(),
        #axis.title.y = element_markdown(),axis.line = element_line(color = 'black'), 
        plot.title = element_text(size = 10), 
        axis.text = element_text(size = 7), 
        axis.text.x = element_text(size = 7, vjust = 0.5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)
    )
    plot_folder = here::here(file.path("tmp/eQTL_GxE/plots", tissue))
    # Check if exists and create plot_folder
    if(!dir.exists(plot_folder)){
    dir.create(plot_folder, recursive = TRUE)
    }
    save_plot(paste0(plot_folder, "/", gene, ".png"), manhplot, base_width = 7)
    return(manhplot)
}
plotManhattan(res, gene, signCut, tissue)
# map2(results, gxe_genes, plotManhattan, signCut, tissue, .progress = "Manhattan Plots")

x = as_tibble(tx) |>
     filter(seqnames %in% paste0("chr", chrs))
x$gene_id = unlist(x$gene_id)
x = x[!duplicated(x$gene_id),]

# Bind all results, filter for significant qvalues at 1%
qtl_pos = bind_rows(results) |>
    filter(q_gxe < 0.005) |>
     arrange(q_gxe) |>
     separate(snp, c("chr", "bp"), sep = "_") |>
     select(Trait:p_main) |>
     left_join(x, by = c("Trait" = "gene_id")) |>
     select(Trait, seqnames, start, chr, bp, p_main) |>
     mutate(seqnames = gsub("chr", "", seqnames)) |>
     mutate(chr = factor(chr, levels = chrs),
            seqnames = factor(seqnames, levels = chrs),
            bp = as.numeric(bp), 
            log10p = -log10(p_main)) |>
     filter(!is.na(seqnames))

 data_cum_snp <- qtl_pos %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(bp)) %>% 
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(chr, bp_add)

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

 data_cum <- qtl_pos %>% 
    group_by(seqnames) %>% 
    summarise(max_start = max(start)) %>% 
    mutate(start_add = lag(cumsum(max_start), default = 0)) %>% 
    select(seqnames, start_add)

qtl_pos <- qtl_pos %>% 
    inner_join(data_cum, by = "seqnames") %>% 
    mutate(start_cum = start + start_add) |>
    dplyr::select(-bp_add, -start_add)

gr <- GRanges(
    seqnames = paste0("chr", qtl_pos$chr),
    ranges = IRanges(qtl_pos$bp, end = qtl_pos$bp, 
                     names = paste(qtl_pos$chr, qtl_pos$bp, sep = "_"), 
                     trait = qtl_pos$Trait)
                     )
genes_gr = subsetByOverlaps(tx, gr)
matched_genes = table(unlist(genes_gr$gene_id)) |> sort()
matched_genes[matched_genes > 15]

chr3L = genes_gr[seqnames(genes_gr) == "chr3L"]
matched_genes = table(unlist(chr3L$gene_id)) |> sort()


hotspot = x |> 
     filter(gene_id %in% names(matched_genes)[matched_genes > 15]) |> 
     mutate(seqnames = gsub("chr", "", seqnames)) 
hotspot$cum_start = hotspot$start +  data_cum_snp[match(hotspot$seqnames, data_cum_snp$chr),]$bp_add
hotspot$cum_end = hotspot$end +  data_cum_snp[match(hotspot$seqnames, data_cum_snp$chr),]$bp_add
hotspot = unique(hotspot)

png("test.png", width = 10, height = 10, units = "in", res = 300)
ggplot(qtl_pos, aes(bp_cum, start_cum, color = chr)) +
    geom_point(size = 0.1, aes(alpha = log10p)) +
    geom_vline(data = hotspot, aes(xintercept= cum_start), linetype = "dashed", linewidth = 0.2) +
     geom_vline(data = hotspot, aes(xintercept= cum_end), linetype = "dashed", linewidth = 0.2) +
     #annotate("text", x = hs_start + 10e5, y = 1000 , label = "pHCl-1", vjust = 1, hjust = 0) +
    theme_classic() +
    labs(x = "eQTL position", y = "Gene position") +
    geom_abline(intercept = 0, slope = 1) 
dev.off()


filter(qtl_pos, bp_cum > hs_start, bp_cum < hs_end) |>
     arrange(log10p) |>
     group_by(Trait) |>
     dplyr::count()

# Create a GRanges object

## Convert chr to numeric


# gemma = import(here::here("eQTLmapping/cache/gemma/head/FBgn0004396/output/gemma.assoc.txt"))

# head(gemma_results)
# head(results)
# gemma_gridlmm = inner_join(results, gemma_results, by = c("X_ID" = "rs")) |>
#      select(X_ID, p_score, p_value_REML.2) |>
#      arrange(p_value_REML.2) |> as_tibble()
# p = ggplot(gemma_gridlmm, aes(x = p_score, y = p_value_REML.2)) +
#      geom_point() +
#      geom_abline(intercept = 0, slope = 1) +
#      scale_x_log10() +
#      scale_y_log10()
# save_plot(here::here("tmp/gemma_gridlmm.png"), p, base_width = 7)

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