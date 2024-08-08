source(here::here("Rscripts/load_install_packages.R"))

future::plan(
  list(future::tweak(future::multisession, workers = 2), 
        future::tweak(future::multisession, workers = 3))
  )

pca <- function(x, space = c("rows", "columns"),
                center = TRUE, scale = FALSE) {
  space <- match.arg(space)
  if (space == "columns") {
    x <- t(x)
  }
  x <- t(scale(t(x), center = center, scale = scale))
  x <- x / sqrt(nrow(x) - 1)
  s <- svd(x)
  loading <- s$u
  colnames(loading) <- paste0("Loading", 1:ncol(loading))
  rownames(loading) <- rownames(x)
  pc <- diag(s$d) %*% t(s$v)
  rownames(pc) <- paste0("PC", 1:nrow(pc))
  colnames(pc) <- colnames(x)
  pve <- s$d^2 / sum(s$d^2)
  if (space == "columns") {
    pc <- t(pc)
    loading <- t(loading)
  }
  return(list(pc = pc, loading = loading, pve = pve, d = s$d^2))
}
pca_plot_outliers = function(resids, color = NULL){
  pca_resid = pca(resids)
  results = t(pca_resid$pc)
  mean = data.frame(matrix(colMeans(results),
                    ncol = dim(results)[2]))
  colnames(mean) = colnames(results)
  if(is.null(color)){
    rpca_resid = PcaGrid(t(resids), 10, crit.pca.distances = 0.99)
    results = data.frame(results, outlier = !rpca_resid@flag)
  } else {
    results = data.frame(results, outlier = color)
  }
  colors = c("black",  "tomato3", viridis(length(unique(color))))
  plot <-  ggplot(results, aes(PC1, PC2, color = outlier)) +
    geom_point() +
    stat_ellipse(level = 0.99, color = "black") +
    stat_ellipse(level = 0.99, type = "norm", color = "black", linetype = 2) +
    geom_point(data= mean, color = "tomato3", size = 7) +
    coord_fixed(ratio=1) +
    labs(x = paste0("PC1 (", round(pca_resid$pve[1]*100, 2), "%)"),
         y = paste0("PC2 (", round(pca_resid$pve[2]*100, 2), "%)")) +
    theme_cowplot() +
    scale_color_manual(values = colors) +
    theme(legend.position = "none")
  plot
}
pca_plot = function(pca_resid, treatment = NULL){
#   pca_resid = pca(resids)
  results = t(pca_resid$pc)
  mean = data.frame(matrix(colMeans(results),
                    ncol = dim(results)[2]))
  colnames(mean) = colnames(results)
  
  results = data.frame(results[, 1:2], treatment = as.factor(treatment))
  colors = c("black",  "tomato3")
  plot <-  ggplot(results, aes(PC1, PC2, color = treatment, group = treatment)) +
    geom_point() +
    stat_ellipse(level = 0.99, color = "black") +
    stat_ellipse(level = 0.99, type = "norm", color = "black", linetype = 2) +
    coord_fixed(ratio=1) +
    labs(x = paste0("PC1 (", round(pca_resid$pve[1]*100, 2), "%)"),
         y = paste0("PC2 (", round(pca_resid$pve[2]*100, 2), "%)")) +
    theme_cowplot() +
    scale_color_manual(values = colors[1:length(unique(treatment))]) 
  plot
}
gt2num <- function(genomat) {
  
  if (!is.matrix(genomat)) {
    stop("genomat parameter must be a matrix")
  }

  ## Define genotypes
  refs <- c("0/0", "0|0")
  alts <- c("1/1", "1|1")
  hets <- c("0/1", "1/0", "0|1", "1|0")
  misses <- c("./.", ".|.")
  
  ## Perform numeric substitution
  genomat[genomat %in% refs] <- 0
  genomat[genomat %in% alts] <- 2
  genomat[genomat %in% hets] <- 1
  genomat[genomat %in% misses] <- NA
  
  ## Record loci with > 2 alleles; remove from matrix
  biallelic <- apply(genomat, 1, function(x) { all(x %in% c(NA, "0", "1", "2")) })
  if (all(biallelic)) {
    removed_loci <- NULL
  } else {
    removed_loci <- rownames(genomat)[!biallelic]
    genomat <- genomat[biallelic, ]
  }
  
  ## Construct and return output list
  class(genomat) <- "numeric"
  return(list("genomat" = genomat, "removed_loci" = removed_loci))
}

# pak::pkg_install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
gene_locations_GR <- transcripts(txdb, columns = "gene_id")
gene_locations_GR$gene_id = unlist(gene_locations_GR$gene_id)
chrs <- c("2L", "2R", "3L", "3R", "4", "X")

# snp_pos = data.frame(rs = snp_locations) |>
#   mutate(pos = rs) |>
#   separate(pos, into = c("chr", "bp"), sep = "_") |>
#   mutate(chr = paste0("chr", as.character(chr))) |>
#   mutate(bp = as.numeric(bp))
# export(snp_pos, file = here::here("cache/snp_pos.rds"))
snp_pos = import(here::here("cache/snp_pos.rds"), trust = TRUE)

getCisSnps = function(transcript, window = 10000, snps = snp_pos){
  snp_GR <- GRanges(
      seqnames = snps$chr,
      ranges = IRanges(snps$bp, end = snps$bp, 
                      names = snps$rs))
  gene_location = gene_locations_GR[which(gene_locations_GR$gene_id == transcript)] + window
  return(snps$rs[unique(findOverlaps(snp_GR, gene_location)@from)])
}

getEnclosingGene <- function(snp, all_snps = snp_pos, window = 0){
  c_snp_pos = all_snps |>
    filter(rs == snp)
  c_snp_GR = GRanges(
    seqnames = c_snp_pos$chr,
    ranges = IRanges(c_snp_pos$bp, end = c_snp_pos$bp, 
                     names = c_snp_pos$rs))
  gene_location = gene_locations_GR + window
  overlap = findOverlaps(c_snp_GR, gene_location)@to
  if(length(overlap) == 0){
    return(NULL)
  } else {
    enclosingGene = unique(gene_location[overlap]$gene_id)
    return(data.frame(snp = snp, index = 1:length(enclosingGene), enclosingGene = enclosingGene))
  }
}

#   c_snp_GR = GRanges(
#     seqnames = snp_pos$chr,
#     ranges = IRanges(snp_pos$bp, end = snp_pos$bp, 
#                      names = snp_pos$rs))
#   gene_location = gene_locations_GR + 0
#   overlap = findOverlaps(c_snp_GR, gene_location)@to
#   enclosingGene = gene_location[overlap]$gene_id

#   snps_with_genes = data.frame(snp = snp_pos$rs[findOverlaps(c_snp_GR, gene_location)@from], enclosingGene = enclosingGene) |> unique()
# export(snps_with_genes, file = here::here("cache/snps_with_genes.tsv"))

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
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1
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
     V_setup = import(paste0(cache_folder, '/V_setup.rds'), trust = TRUE)
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

# convert flybase id to gene symbol
fly2sym = function(flybase){
bitr(flybase, 
             fromType="FLYBASE", 
             toType = "SYMBOL",
             OrgDb = org.Dm.eg.db, 
             drop = FALSE)
}

plot_GxE_eQTL <- function(current_snp, gene = NULL, snp_table = sig_results_df){

     if(is.null(gene)){
          current_gene = snp_table |> 
               filter(snp == current_snp) |> 
               select(Trait) |> 
               pull()
     } else {
          current_gene = gene
     }
     if(length(current_gene) > 1){
         out_plots = vector("list", length = length(current_gene))
         for(i in seq_along(current_gene)){
             out_plots[[i]] = plot_GxE_eQTL(current_snp, gene = current_gene[i], snp_table)
         }
         return(out_plots)
     }

     pheno = rnaseq[[tissue]][["mwash_residuals"]][[data_set]][current_gene, covariates$id]
     data_rna = covariates
     data_rna$y = pheno
     m1 = lm(global_formulas_gemma[[tissue]], data = data_rna) 
     data_rna$residuals = residuals(m1)
     data_rna$snp = as.factor(X[covariates$id,current_snp])
     data_rna$treatment = c("ctrl", "hs")[data_rna$treatment + 1]

     out_folder = here::here(paste0("output/snp_validation/", tissue, "/", data_set))
     if(!dir.exists(out_folder)){
          dir.create(out_folder, recursive = TRUE)
     }
     
     enclosing_gene = filter(snp_table, snp == current_snp) |> 
          select(enclosingGene) |> 
          pull()
     if(!grepl("FBgn", enclosing_gene[[1]])){
          enclosing_gene = "No gene"
     } else {
          enclosing_gene = paste("Enclosing gene:", fly2sym(enclosing_gene)$SYMBOL[[1]])
     }
     current_symbol = tryCatch(fly2sym(current_gene)$SYMBOL[[1]], error = function(e) "No symbol")

     cis = filter(snp_table, snp == current_snp, Trait == current_gene) |> 
          select(cis) |> 
          pull()
     cis = c("trans", "cis")[cis + 1]
     if(length(cis) > 1){
          cis = cis[1]
     }
     out_file = paste0(out_folder, "/", current_snp, "_", current_gene, "_", cis, ".png")

     # ggplot(data_eqtl, aes(x = snp, y = residuals, color = treatment, group = interaction(snp, treatment))) +
     #      geom_jitter(position=position_jitterdodge(dodge.width=0.75, jitter.width = 0.1)) +
     #      geom_boxplot(color = "black", fill = "transparent", outlier.shape = NA) +
     #      labs(x = paste("SNP -", current_snp), y = paste("Expression -", current_gene)) +
     #      theme_classic() + ggtitle("log2 counts") + 
     out_plot = ggplot(data_rna, 
                       aes(x = treatment, y = residuals, 
                           fill = snp, 
                           group = interaction(snp, treatment))) +
          geom_jitter(position=position_jitterdodge(dodge.width=0.75, jitter.width = 0.2), shape = 21, size = 1) +
          geom_boxplot(color = "black", fill = "transparent", outlier.shape = NA) +
          scale_fill_manual(values = c("#edae49", "#d1495b", "#00798c"), 
                             name = paste0("SNP: ", current_snp, "\n", 
                                           enclosing_gene, 
                                           "\nIn ", cis)) +
          labs(x = "Condition", y = paste("Expression of", current_gene, "-", current_symbol)) +
          theme_cowplot(12)
     save_plot(out_file, out_plot, base_width = 7, base_height = 5)
     return(out_plot)
}

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