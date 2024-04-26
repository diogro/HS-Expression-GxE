source(here::here("Rscripts/functions.R"))

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

covariates = import(here::here(paste0("eQTLmapping/covariates/", tissue, ".tsv")))
GRM = import(here::here(paste0("eQTLmapping/GRMs/", tissue, ".cXX.txt")), header = FALSE)
colnames(GRM) = rownames(GRM) = covariates$id
genes = import(here::here(paste0("eQTLmapping/phenotypes/", tissue, ".genes.txt")), header = FALSE)[,1]
gxe_genes = gxe_genes[[tissue]]

current_gene = genes[100]

global_formulas <- list(
          head = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 + Zhat2 + (1|id),
          body = 
     y ~ 1 + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + treatment + Zhat1 +         (1|id)
)

runGxEmodel = function(current_gene, tissue, covariates, GRM){
     old_dir = getwd()
     setwd(here::here("eQTLmapping"))
     cache_folder = here::here(paste0('eQTLmapping/cache/'),
                           tissue, '/',
                           current_gene)
     y = t(import(here::here(paste0("eQTLmapping/phenotypes/", tissue, ".tsv")), 
                  skip = which(genes == current_gene)-1, nrows = 1))
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
     qvalues = qvalue(out_file$p_gxe, fdr.level = 0.01)
     out_file = out_file |> 
          mutate(q_gxe = qvalues$qvalues) 
     return(out_file)
}

results = vector("list", length = length(gxe_genes))
k = 1

for(i in k:length(gxe_genes)){
    print(paste0("Current gene: ", gxe_genes[i], " Percent: ", 100*round(i/length(gxe_genes), 3), "%"))
    results[[k]] = runGxEmodel(gxe_genes[i], tissue, covariates, GRM)
    k = k + 1
}
export(results, here::here(paste0("cache/eQTL_detections_gxe-", tissue, ".rds")))
results = import(here::here(paste0("cache/eQTL_detections_gxe-", tissue, ".rds"))) 

x = results[[1]]
any(x$q_gxe < 0.1)


signCut <- 1e-4

# Make a manhattan plot
res = results[[1]]
gene = gxe_genes[1]
plotManhattan <- function(res, gene, signCut, tissue){
    res = res |>
        mutate(log_p_gxe = -log10(p_gxe)) |>
        separate(snp, c("chr", "bp"), sep = "_") |>
        mutate(chr = factor(chr, levels = chrs),
                bp = as.numeric(bp)) |>
        arrange(p_gxe) 

    gwas_data_load = res
    sig_data <- gwas_data_load %>% 
    subset(p_gxe < signCut)
    notsig_data <- gwas_data_load %>% 
    subset(p_gxe >= signCut) %>%
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

    manhplot <- ggplot(filter(gwas_data, p_gxe < signCut, chr != 4), aes(x = bp_cum, y = -log10(p_gxe), 
                                    color = chr)) +
    geom_hline(yintercept = -log10(signCut), color = "grey40", linetype = "dashed") + 
    geom_point(alpha = 0.7, size = 0.005) +
    geom_point(data = filter(gwas_data, p_gxe > signCut), alpha = 0.75, size = 0.005, color = "grey") +
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
map2(results, genes, plotManhattan, signCut, tissue, .progress = "Manhattan Plots")





