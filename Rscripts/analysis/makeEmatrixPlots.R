source(here::here("Rscripts/functions.R"))

getBlockSizedf = function(level, block_df, all = FALSE, draw_level=level){
  upper = paste0("B", level+1)
  if(all){
    s_block_df = unique(block_df[,paste0("B", draw_level:(level+1))])
  } else
    s_block_df = unique(block_df[,paste0("B", level:(level+1))])

  t_upper = table(s_block_df[,upper])
  block_sizes = t_upper[match(unique(as.data.frame(s_block_df)[,upper]), names(t_upper))]
  sum_block_sizes = cumsum(block_sizes)
  b_size_mat = matrix(NA, nrow(block_sizes), 2)
  for(i in 1:nrow(block_sizes)){
    if(i == 1)
      b_size_mat[i,] = c(0.5, sum_block_sizes[i]+0.5)
    else
      b_size_mat[i,] = c(sum_block_sizes[i-1]+0.5, sum_block_sizes[i]+0.5)
  }
  b_size_df = data.frame(start = b_size_mat[,1],
                         end = b_size_mat[,2])
  return(b_size_df)
}
t = "head"
d = "decohere"
fdr = "1e2"
nfdr = 0.01
makeEmatrixPlots = function(t, d, fdr, 
                            data_path = here::here("cache/clip"),
                            plot_path = here::here("tmp"),
                            levels = 5){
  block_df = import(file.path(data_path, paste0("clip_fdr-", fdr,"_", t, ".csv"))) |> filter(direction == d)
  block_summary =  import(file.path(data_path, paste0("blockSummary_clip_fdr-", fdr,"_", t, ".csv"))) |> filter(Direction == d) %>%
    filter(Nested_Level == 1) %>% select(Block, Internal_degree, Assortativity) %>%
    rename(B1 = Block)
  block_df = inner_join(block_df, block_summary)
  block_df = block_df %>%
    arrange(B6, B5, B4, B3, B2, desc(Assortativity)) %>%
    select(-1)

  e_mats = vector("list", levels)
  for (i in 1:levels){
    e_matrix = read_csv(file.path(data_path, paste0("Ematrices/fdr-0.01/", t, "/", d, "_E_matrix_level",i-1,".csv")))
    e_matrix[,1] = NULL
    e_matrix = as.matrix(e_matrix)
    rownames(e_matrix) = colnames(e_matrix) = as.character(colnames(e_matrix))
    blocks = as.character(unique(as.data.frame(block_df)[,paste0("B", i)]))
    e_matrix= e_matrix[blocks, blocks]
    #e_matrix = e_matrix[-61, -61]
    diag(e_matrix) = diag(e_matrix)/2
    e_mats[[i]] = e_matrix
  }
  makePlotE = function(level, e_mats, block_df){
    e_matrix = e_mats[[level]]
    colnames(e_matrix) = rownames(e_matrix) = paste0("b", rownames(e_matrix))
    melted_cormat <- reshape2::melt(e_matrix)
    melted_cormat[melted_cormat == 0] = NA
    #plot heatmap
    plot = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
      geom_tile() +
      scale_fill_viridis_c(alpha = 1)  + 
      theme_cowplot() +
      labs(y = "Level-1 Blocks", x = paste0("Level-", level," Blocks")) + 
      ggtitle(paste("Level", level)) +
      scale_y_discrete(label = gsub("b", "", colnames(e_matrix))) + 
      scale_x_discrete(label = gsub("b", "", colnames(e_matrix)), guide = guide_axis(n.dodge = 2)) + 
      theme(axis.text.x = element_text(4, angle=35, vjust=0.3), 
            legend.position = "bottom", 
            legend.key.width= unit(22, 'mm'), 
            legend.key.height= unit(2.2, 'mm'), 
            legend.title = element_blank())
    if(level < levels){
      for(i in level:(levels-1)){
        b_size_df = getBlockSizedf(i, block_df, all = TRUE, level)
        plot = plot + geom_rect(data = b_size_df, color = "tomato3", alpha = 0, size = 0.5,
                                aes(x = NULL, y = NULL, fill = NULL, xmin=start, xmax=end,
                                    ymin=start, ymax=end))
      }
    }
    return(plot)
  }
  makePlotDeg = function(level, block_df){
    block_df[[paste0("B", level)]] = factor(block_df[[paste0("B", level)]], 
                                            level = unique(block_df[[paste0("B", level)]]))
    plot = ggplot(data = block_df, aes(x=B1, y=Degree)) + geom_boxplot() + 
        scale_y_continuous(breaks = c(0, 1000, 2000), minor_breaks = c(500, 1500)) +
        coord_flip() + theme_cowplot() + background_grid(minor = "xy") + easy_remove_axes() + 
        theme(axis.title.x = element_text(8), axis.text.x = element_text(4))
    return(plot)
  }
  plot_list = lapply(seq_along(1:levels), makePlotE, e_mats, block_df)
  plot_list_deg = lapply(seq_along(1:levels), makePlotDeg, block_df)
  layout = 
"
AA
AA
BC"
  all_plots = plot_list[[1]] + plot_list[[2]] + plot_list[[3]] +  plot_layout(design = layout)
  size = 6
  save_plot(file.path(plot_path, paste0("E_matrices-", fdr, "-", t, "-", d, ".png")), all_plots,
            base_height = size*4, base_width = size*3)
  return(list(df = block_df, E = e_mats, plots = all_plots, plot_list = plot_list, plot_list_deg = plot_list_deg))
}

for(d in c("decohere", "integrate")){
  makeEmatrixPlots("head", d, "1e2", levels = 5)
  makeEmatrixPlots("body", d, "1e2", levels = 5)
}

out_fdr_1e3_body = makeEmatrixPlots("body_weights-spearman_fdr-1e-03_mcmc_mode", levels = 5)
