results = list(vicar = import(here::here("output/HS-ctrl-DE_vicar-table_2023-11-17.csv")),
               vicar2 = import(here::here("output/HS-ctrl-DE_vicar-head2-body2-table_2023-11-17.csv"))) 

names(results$vicar)

table = list_rbind(list(
    vicar = results$vicar |>
        rename(effect = PosteriorMean) |>
        select(tissue, Geneid, SYMBOL, effect, qvalue),
    vicar2 = results$vicar |>
        rename(effect = PosteriorMean) |>
        select(tissue, Geneid, SYMBOL, effect, qvalue)), names_to = "method") %>%
        filter(tissue == "body")

wide_table = table |>
    pivot_wider(values_from = c("effect", "qvalue"), names_from = "method") 

color_points = wide_table |>
    filter(tissue == "body") %>%
    filter(-log10(qvalue_vicar) > 1.5)

png("tmp/effects_vicar-vicar2.png", width = 2000, height = 2000)
p1 = ggplot(wide_table, aes(effect_vicar, effect_vicar2)) + 
    geom_point(alpha= 0.3) + 
    geom_point(data = color_points, 
               alpha = 1, color = 2) + 
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~tissue) +
    theme_cowplot() + theme(legend.position = "none") 
p2 = ggplot(wide_table, aes(-log10(qvalue_vicar), -log10(qvalue_vicar2))) + 
    geom_point(alpha= 0.3) + 
    geom_point(data = color_points, 
               alpha = 1, color = 2) + 
    scale_color_manual(values = c("black", 2)) +
    facet_wrap(~tissue, scales = "free_y") +
    theme_cowplot() + theme(legend.position = "none") 
 (p1 + ggtitle("Effects")) / (p2 + ggtitle("qvalues"))
dev.off()
