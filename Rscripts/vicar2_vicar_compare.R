results = list(vicar = import(here::here("output/HS-ctrl-DE_vicar-table_2023-11-17.csv")),
               vicar2 = import(here::here("output/HS-ctrl-DE_vicar_body-2ruv-table_2023-11-17.csv"))) 

names(results$vicar)

table = list_rbind(list(
    vicar = results$vicar |>
        rename(effect = PosteriorMean) |>
        select(tissue, Geneid, SYMBOL, effect, qvalue),
    vicar2 = results$vicar2 |>
        rename(effect = PosteriorMean) |>
        select(tissue, Geneid, SYMBOL, effect, qvalue)), names_to = "method") %>%
        filter(tissue == "body")

wide_table = table |>
    pivot_wider(values_from = c("effect", "qvalue"), names_from = "method") 

color_points = wide_table |>
    filter(-log10(qvalue_vicar) < 4 & -log10(qvalue_vicar2) > 4)


p1 = ggplot(wide_table, aes(effect_vicar, effect_vicar2)) + 
    geom_point(alpha= 0.3, size = 0.3) + 
    geom_point(data = color_points, size = 0.6,
               alpha = 0.5, color = 2) + 
    geom_abline(intercept = 0, slope = 1) +
    theme_bw() + theme(legend.position = "none") 
p2 = ggplot(wide_table, aes(-log10(qvalue_vicar), -log10(qvalue_vicar2))) + 
    geom_point(alpha= 0.3, size = 0.3) + 
    geom_point(data = color_points, size = 0.6,
               alpha = 0.5, color = 2) + 
    scale_color_manual(values = c("black", 2)) +
    theme_bw() + theme(legend.position = "none") 
panel = (p1 + ggtitle("Effects")) / (p2 + ggtitle("qvalues"))
save_plot("tmp/effects_vicar-vicar2.png",panel , base_width = 5, base_height = 7)


table = results$vicar |> 
    filter(tissue == "body")
labels = table  %>%
    filter(NegativeProb > 0.999 | PositiveProb > 0.999) |>
    select(-NegativeProb, -PositiveProb) |>
    arrange(qvalue) %>% 
    filter(-log10(qvalue) > 4)

p1 = ggplot(table,
        aes(PosteriorMean, -log10(qvalue))) +
    geom_point(size = 0.1, alpha = 0.1, color = "gray") + 
    geom_point(data = table |> filter(NegativeProb > 0.999 | PositiveProb > 0.999), 
               size = 0.1, alpha = 0.5, color = "blue") + 
    geom_text_repel(data = labels, 
                    aes(label = SYMBOL), 
                    size = 2, max.overlaps = 30, segment.size = 0.1) + 
    scale_color_brewer(palette = "Set1") +
    theme_bw() + 
    theme(legend.position='none')

table = results$vicar2 |> 
    filter(tissue == "body")
labels = table  %>%
    filter(NegativeProb > 0.999 | PositiveProb > 0.999) |>
    select(-NegativeProb, -PositiveProb) |>
    arrange(qvalue) %>% 
    filter(-log10(qvalue) > 4)

p2 = ggplot(table,
        aes(PosteriorMean, -log10(qvalue))) +
    geom_point(size = 0.1, alpha = 0.1, color = "gray") + 
    geom_point(data = table |> filter(NegativeProb > 0.999 | PositiveProb > 0.999), 
               size = 0.1, alpha = 0.5, color = "blue") + 
    geom_text_repel(data = labels, 
                    aes(label = SYMBOL), 
                    size = 2, max.overlaps = 30, segment.size = 0.1) + 
    scale_color_brewer(palette = "Set1") +
    theme_bw() + 
    theme(legend.position='none') 

 panel = (p1 + ggtitle("Body - vicar 1 RUV")) / (p2 + ggtitle("Body - vicar 2 RUV"))
save_plot("tmp/DE_volcano_vicar2-vs-vicar-body.png",panel , base_width = 5, base_height = 7)
