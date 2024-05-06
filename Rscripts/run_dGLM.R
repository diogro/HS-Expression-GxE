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
            trait = colnames(phenos)[current_gene])
    names(df) = c("estimate", "std.error", "z.value", "p.value", "name", "trait")
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
            trait = colnames(phenos)[current_gene])
    names(df) = c("estimate", "std.error", "z.value", "p.value", "name", "trait")
    df
}
safe_runDGLM = possibly(.f = rundGLM, otherwise = NULL)
dGLM_body = map_df(1:ncol(phenos), safe_runDGLM, .progress = 'dGLM')
length(unique(dGLM_body$trait)) 

bind_rows(list(head = dGLM_head, body = dGLM_body), .id = 'tissue', )  |> 
    as_tibble() |>
    select (tissue, trait, name, everything()) |> 
    export(affix_date(here::here("output/dGLM.csv")))
