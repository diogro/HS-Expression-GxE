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
