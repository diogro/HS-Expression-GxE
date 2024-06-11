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
snp_pos = import(here::here("cache/snp_pos.rds"))

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