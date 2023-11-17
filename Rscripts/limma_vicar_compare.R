results = list(limma = import(here::here("output/HS-ctrl-DE_limma-table_2023-11-17.csv")),
               vicar = import(here::here("output/HS-ctrl-DE_vicar-head2-body2-table_2023-11-17.csv"))) 

names(results$vicar)
names(results$limma)
table = list_rbind(list(
    limma = results$limma |>
        dlply("tissue") |>
        map(\(x) mutate(x, qvalue = qvalue(P.Value)$qvalue)) |> list_rbind() |>
        rename( effect = logFC) |>
        select(tissue, Geneid, SYMBOL, effect, qvalue),
    vicar = results$vicar |>
        rename(effect = PosteriorMean) |>
        select(tissue, Geneid, SYMBOL, effect, qvalue)), names_to = "method")

wide_table = table |>
    pivot_wider(values_from = c("effect", "qvalue"), names_from = "method") 

color_points = wide_table |> 
    filter(tissue == "head") %>%
    filter(-log10(qvalue_vicar) > 10) |>
bind_rows(wide_table |> 
    filter(tissue == "body") %>%
    filter(-log10(qvalue_vicar) > 1.5))

png("tmp/effects_vicar-limma.png", width = 2000, height = 2000)
p1 = ggplot(wide_table, aes(effect_limma, effect_vicar)) + 
    geom_point(alpha= 0.3) + 
    geom_point(data = color_points, 
               alpha = 1, color = 2) + 
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~tissue) +
    theme_cowplot() + theme(legend.position = "none") 
p2 = ggplot(wide_table, aes(-log10(qvalue_limma), -log10(qvalue_vicar))) + 
    geom_point(alpha= 0.3) + 
    geom_point(data = color_points, 
               alpha = 1, color = 2) + 
    scale_color_manual(values = c("black", 2)) +
    facet_wrap(~tissue, scales = "free_y") +
    theme_cowplot() + theme(legend.position = "none") 
 (p1 + ggtitle("Effects")) / (p2 + ggtitle("qvalues"))
dev.off()


table = results$vicar
label_head = table |> 
    filter(tissue == "head") %>%
    filter(NegativeProb > 0.95 | PositiveProb > 0.95) |>
    select(-NegativeProb, -PositiveProb) |>
    arrange(qvalue) %>% 
    filter(-log10(qvalue) > 10)
label_body = table |> 
    filter(tissue == "body") %>%
    filter(NegativeProb > 0.999 | PositiveProb > 0.999) |>
    select(-NegativeProb, -PositiveProb) |>
    arrange(qvalue) %>% 
    filter(-log10(qvalue) > 4)
labels = bind_rows(label_head, label_body )

{p1 = ggplot(table,
        aes(PosteriorMean, -log10(qvalue))) +
    geom_point(size = 0.1, alpha = 0.1, color = "gray") + 
    geom_point(data = table |> filter(NegativeProb > 0.999 | PositiveProb > 0.999), 
               size = 0.1, alpha = 0.5, color = "blue") + 
    geom_text_repel(data = labels, 
                    aes(label = SYMBOL), 
                    size = 2, max.overlaps = 30, segment.size = 0.1) + 
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~tissue, scales = "free_y") + theme_bw() + 
    theme(legend.position='none') } 
{save_plot("tmp/DE_volcano_vicar.png", p1, base_width = 7)}

table = results$limma |> dlply("tissue") |> map(\(x) mutate(x, qvalue = qvalue(P.Value)$qvalue)) |> list_rbind()
label_head = table |> 
    filter(tissue == "head",
           -log10(qvalue) > 11)
label_body = table |> 
    filter(tissue == "body",
           -log10(qvalue) > 15)
labels = bind_rows(label_head, label_body )
names(table)
{p2 = ggplot(table, 
        aes(logFC, -log10(qvalue))) +
    geom_point(size = 0.1, alpha = 0.1, color = "gray") + 
    geom_point(data = table |> filter(HSvsC != 0), size = 0.1, alpha = 0.5, color = "blue") + 
    geom_text_repel(data = labels, 
                    aes(label = SYMBOL), 
                    size = 2, max.overlaps = 30, segment.size = 0.1) + 
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~tissue) + theme_bw() + 
    theme(legend.position='none') }
{save_plot("tmp/DE_volcano_limma.png", p2, base_width = 7)}

panel = (p1 + ggtitle("Vicar")) / (p2 + ggtitle("Limma/Voom"))
save_plot("tmp/DE_volcano_limma-vs-viacar.png",panel , base_width = 7, base_height = 7)


