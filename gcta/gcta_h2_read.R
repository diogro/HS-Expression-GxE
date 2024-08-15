library(tidyverse)
library(cowplot)
library(rio)


getGCTA_output_single = function(current_file, file){
    variance_table = matrix(NA, 4, 2)
    variance_table[1,] = as.numeric(current_file[5:6])
    variance_table[2,] = as.numeric(current_file[8:9])
    variance_table[3,] = as.numeric(current_file[11:12])
    variance_table[4,] = as.numeric(current_file[14:15])
    rownames(variance_table) = current_file[c(4, 7, 10, 13)]
    colnames(variance_table) = c("Estimate", "SE")
    stats = c(logL = current_file[17],
      logL0 = current_file[19],
      LRT = current_file[21], 
      df = current_file[23], 
     Pval = current_file[25], 
     n = current_file[27])
    list(var = variance_table, stats = stats)
}
getGCTA_output_all = function(current_file, file){
    variance_table = matrix(NA, 7, 2)
    variance_table[1,] = as.numeric(current_file[5:6])
    variance_table[2,] = as.numeric(current_file[8:9])
    variance_table[3,] = as.numeric(current_file[11:12])
    variance_table[4,] = as.numeric(current_file[14:15])
    variance_table[5,] = as.numeric(current_file[17:18])
    variance_table[6,] = as.numeric(current_file[20:21])
    variance_table[7,] = as.numeric(current_file[25:26])
    rownames(variance_table) = c(current_file[c(4, 7, 10, 13, 16, 19)], "Sum.V(G)/Vp")
    colnames(variance_table) = c("Estimate", "SE")
    stats = c(logL = current_file[28],
              logL0 = current_file[30],
              LRT = current_file[32], 
              df = current_file[34], 
              Pval = current_file[36], 
              n = current_file[38])
    list(var = variance_table, stats = stats)
}
getGCTA_output = function(gene, treatment = "ctrl", 
                          results_folder = "/Genomics/argo/users/damelo/projects/HS-Expression-GxE/output/gcta"){
    file = r_meta_df_gcta[r_meta_df_gcta$gene == gene & 
                          r_meta_df_gcta$treatment == treatment,
                          "file"]
    current_file = scan(file.path(results_folder, file), character())
    if(treatment %in% c("hs", "ctrl")) return(getGCTA_output_single(current_file, file))
    if(treatment == "hsctrl") return(getGCTA_output_all(current_file, file))
}


results_folder = "/Genomics/argo/users/damelo/projects/HS-Expression-GxE/output/gcta"
result_files = dir(results_folder, pattern = ".hsq")
results_meta = strsplit(result_files, "_")

index = as.numeric(sapply(results_meta, `[[`, 1))
gene = sapply(results_meta, `[[`, 2)
treatment = sapply(strsplit(sapply(results_meta, `[[`, 3), "\\."), `[[`, 1)
r_meta_df_gcta = data.frame(index, gene, treatment, file = result_files)

models_ctrl = lapply(unique(r_meta_df_gcta$gene[r_meta_df_gcta$treatment=="ctrl"]), getGCTA_output)
models_HiSu = lapply(unique(r_meta_df_gcta$gene[r_meta_df_gcta$treatment=="hs"]), getGCTA_output, "hs")
models_HSxC = lapply(unique(r_meta_df_gcta$gene[r_meta_df_gcta$treatment=="hsctrl"]), getGCTA_output, "hsctrl")
names(models_ctrl) = unique(r_meta_df_gcta$gene[r_meta_df_gcta$treatment=="ctrl"])
names(models_HiSu) = unique(r_meta_df_gcta$gene[r_meta_df_gcta$treatment=="hs"])
names(models_HSxC) = unique(r_meta_df_gcta$gene[r_meta_df_gcta$treatment=="hsctrl"])

length(models_ctrl)
length(models_HiSu)
length(models_HSxC)

p = map_df(models_HSxC, function(x) x$var[,1]) %>%
      mutate(.id = names(models_HSxC)) |>
      select(.id, 'V(G)', 'V(GxE)') %>%
      pivot_longer(cols = 2:ncol(.)) %>%
      ggplot(aes(name, value)) + 
      geom_boxplot() + theme_minimal_grid() + 
      labs(y = "Estimates", x = "Variance Component")
save_plot("test.png", p)

p = map_df(models_HSxC, function(x) x$var[,1]) %>%
      mutate(.id = names(models_HSxC)) |>
      select(.id, 'V(G)/Vp', 'V(GxE)/Vp') %>%
      pivot_longer(cols = 2:ncol(.)) %>%
      ggplot(aes(name, value)) + 
      geom_boxplot() + theme_minimal_grid() + 
      labs(y = "Estimates", x = "Variance Component")
save_plot("test.png", p)

models_ctrl[[1]]$var
h2_by_condition = 
    map_df(models_ctrl, function(x) x$var[,1]) |>
        mutate(Gene = names(models_ctrl), Condition = "ctrl") |> 
    bind_rows(
        map_df(models_HiSu, function(x) x$var[,1]) |>
            mutate(Gene = names(models_ctrl), Condition = "hs")) |>
    rename(h2  = `V(G)/Vp`) |>
    mutate(Tissue = "head") |>
    relocate(Gene, Tissue, Condition)
export(h2_by_condition, "output/h2_head_hs-ctrl.tsv")

