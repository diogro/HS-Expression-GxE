source(here::here("Rscripts/functions.R"))

if(!require('dglm')) {
  install.packages('dglm')
  library('dglm')
}

rnaseq = import(here::here("cache/rnaseq_all_2024-03-21.rds"))

t = 'head'

phenos = t(rnaseq[[t]]$l2c)
rownames(phenos)
covariates = rnaseq[[t]]$covariates

colnames(rnaseq[[t]]$mwash$l2c$Zhat) = paste0("Zhat", 1:dim(rnaseq[[t]]$mwash$l2c$Zhat)[2])
input_df = cbind(covariates, 
                 rnaseq[[t]]$mwash$l2c$Zhat)

stopifnot(rownames(phenos) == input_df$id)

current_gene = 764
rundGLM = function(current_gene){
    data = input_df |>
    mutate(y = phenos[,current_gene]) |>
    as.data.frame()

    m = dglm(y ~ 1 + treatment + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + Zhat1 + Zhat2, 
        ~ 1 + treatment + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + Zhat1 + Zhat2,
        data = data, 
        dlink = "log",
        verbose = FALSE)
    s = summary(m)
    df = data.frame(rbind(s$coefficients[2,] , 
                    s$dispersion.summary$coefficient[2,]), 
            name = c("mean", "dispersion"), 
            Geneid = colnames(phenos)[current_gene])
    names(df) = c("estimate", "std.error", "z.value", "p.value", "name", "Geneid")
    df
}
safe_runDGLM = possibly(.f = rundGLM, otherwise = NULL)
dGLM_head = map_df(1:ncol(phenos), safe_runDGLM, .progress = 'dGLM')

t = 'body'

phenos = t(rnaseq[[t]]$l2c)
rownames(phenos)
covariates = rnaseq[[t]]$covariates

colnames(rnaseq[[t]]$mwash$l2c$Zhat) = paste0("Zhat", 1:dim(rnaseq[[t]]$mwash$l2c$Zhat)[2])
input_df = cbind(covariates, 
                 rnaseq[[t]]$mwash$l2c$Zhat)

stopifnot(rownames(phenos) == input_df$id)

current_gene = 765
rundGLM = function(current_gene){
    data = input_df |>
    mutate(y = phenos[,current_gene]) |>
    as.data.frame()

    m = dglm(y ~ 1 + treatment + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + Zhat1, 
        ~ 1 + treatment + egglayBatch + RNAseqBatch + platingBatch + RNAlibBatch + Zhat1,
        data = data, 
        dlink = "log",
        verbose = FALSE)
    s = summary(m)
    df = data.frame(rbind(s$coefficients[2,] , 
                    s$dispersion.summary$coefficient[2,]), 
            name = c("mean", "dispersion"), 
            Geneid = colnames(phenos)[current_gene])
    names(df) = c("estimate", "std.error", "z.value", "p.value", "name", "Geneid")
    df
}
safe_runDGLM = possibly(.f = rundGLM, otherwise = NULL)
dGLM_body = map_df(1:ncol(phenos), safe_runDGLM, .progress = 'dGLM')
length(unique(dGLM_body$Geneid)) 

dGLM = bind_rows(list(head = dGLM_head, body = dGLM_body), .id = 'tissue', )  |> 
    as_tibble() |>
    rename(Geneid = trait) |>
    select (tissue, Geneid, name, everything()) 
export(dGLM, affix_date(here::here("output/dGLM.csv")))

results = list(limma = import(here::here("output/HS-ctrl-DE_limma-table_2024-03-21.csv")),
               vicar = import(here::here("output/HS-ctrl-DE_vicar-table_2024-03-21.csv")),
               DESeq = import(here::here("output/HS-ctrl-DE_DESeq2-table_2024-03-21.csv"))) 

limma_dglm = results$limma |> 
    select(tissue, Geneid, SYMBOL, logFC, P.Value) |>
    as_tibble() |>
    inner_join(dGLM |> filter(name == "mean"), by = c("tissue", "Geneid")) 

p = ggplot(limma_dglm, aes(x = estimate, y = logFC)) + 
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, color = "red") + 
    facet_wrap(~tissue) + 
    theme_minimal() + 
    labs(x = "dGLM estimate", y = "limma logFC")
save_plot(affix_date(here::here("output/plots/limma_dglm.png")), p)


vicar_dglm = results$vicar |> 
    select(tissue, Geneid, SYMBOL, betahat, lfdr) |>
    as_tibble() |>
    inner_join(dGLM |> filter(name == "mean"), by = c("tissue", "Geneid")) 

p = ggplot(vicar_dglm, aes(x = estimate, y = betahat)) + 
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, color = "red") + 
    facet_wrap(~tissue) + 
    theme_minimal() + 
    labs(x = "dGLM estimate", y = "vicar logFC")
save_plot(affix_date(here::here("output/plots/vicar_dglm.png")), p)

DESeq_dglm = results$DESeq |> 
    select(tissue, Geneid, SYMBOL, log2FoldChange, padj) |>
    as_tibble() |>
    inner_join(dGLM |> filter(name == "mean"), by = c("tissue", "Geneid"))

p = ggplot(DESeq_dglm, aes(x = estimate, y = log2FoldChange)) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, color = "red") + 
    facet_wrap(~tissue) + 
    theme_minimal() + 
    labs(x = "dGLM estimate", y = "DESeq logFC")
save_plot(affix_date(here::here("output/plots/DESeq_dglm.png")), p)
