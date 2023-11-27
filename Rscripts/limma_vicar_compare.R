source(here::here("Rscripts/functions.R"))

results = list(limma = import(here::here("output/HS-ctrl-DE_limma-table_2023-11-17.csv")),
               vicar = import(here::here("output/HS-ctrl-DE_vicar-table_2023-11-17.csv")),
               DESeq = import(here::here("output/HS-ctrl-DE_DESeq2-table_2023-11-18.csv"))) 

table = list_rbind(list(
    limma = results$limma |>
        dlply("tissue") |>
        map(\(x) mutate(x, qvalue = qvalue(P.Value)$qvalue)) |> list_rbind() |>
        rename( effect = logFC) |>
        select(tissue, Geneid, SYMBOL, effect, qvalue),
    vicar = results$vicar |>
        rename(effect = PosteriorMean) |>
        select(tissue, Geneid, SYMBOL, effect, qvalue),
    DESeq = results$DESeq |>
        dlply("tissue") |>
        map(\(x) mutate(x, qvalue = qvalue(pvalue)$qvalue)) |> list_rbind() |>
        rename( effect = log2FoldChange) |>
        select(tissue, Geneid, SYMBOL, effect, qvalue)), names_to = "method")

wide_table = table |>
    pivot_wider(values_from = c("effect", "qvalue"), names_from = "method") 

color_points = wide_table |> 
    filter(tissue == "head") %>%
    mutate(lqv = -log10(qvalue_vicar),
           lql = -log10(qvalue_limma) ) %>%
    mutate(type = case_when(
            lql > 15 & lqv < 5 ~ "limma",
            lql < 15 & lqv > 5 ~ "vicar",
            lql > 15 & lqv > 5 ~ "both",
            TRUE ~ NA
  )) %>%
bind_rows(wide_table |> 
    filter(tissue == "body") %>%
    mutate(lqv = -log10(qvalue_vicar),
           lql = -log10(qvalue_limma) ) %>%
    mutate(type = case_when(
            lql > 15 & lqv < 1.5 ~ "limma",
            lql < 15 & lqv > 1.5 ~ "vicar",
            lql > 15 & lqv > 1.5 ~ "both",
            TRUE ~ NA
  ))) %>% filter(!is.na(type))

p1 = ggplot(wide_table, aes(effect_limma, effect_vicar)) + 
    geom_point(alpha= 0.3, size = 0.3) + 
        geom_point(data = color_points, size = 0.6,
               alpha = 0.5, aes(color = type)) + 
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~tissue) +
    theme_bw() + theme(legend.position = "none") 
p2 = ggplot(wide_table, aes(-log10(qvalue_limma), -log10(qvalue_vicar))) + 
    geom_point(alpha= 0.3, size = 0.3) + 
    geom_point(data = color_points, size = 0.6,
               alpha = 0.5, aes(color = type)) + 
    facet_wrap(~tissue, scales = "free_y") +
    theme_bw() + theme(legend.position = "none") 
panel = (p1 + ggtitle("Effects")) / (p2 + ggtitle("qvalues"))
save_plot("tmp/effects_vicar-limma.png",panel , base_width = 5, base_height = 7, nrow = 2, ncol = 2)

x = "DESeq" 
y = "vicar"
t = "head"
plotEffects <- function(x, y, t){
    data  = filter(wide_table, tissue == t)
    effect_x = rlang::sym(paste0("effect_", x))
    effect_y = rlang::sym(paste0("effect_", y))
    p = ggplot(data, aes(!!effect_x, !!effect_y)) + 
        geom_point(alpha = 0.3, size = 0.5) +
        geom_abline(intercept = 0, slope = 1) +
        theme_bw() + theme(legend.position = "none") 
    p
}

p = plotEffects("limma", "vicar", "body") + 
    plotEffects("DESeq", "vicar", "body") + 
    plotEffects("DESeq", "limma", "body")
save_plot("tmp/test.png",p , base_width = 7, base_height = 7, ncol = 3)
p2 = plotEffects("limma", "vicar", "head") + 
    plotEffects("DESeq", "vicar", "head") +
    plotEffects("DESeq", "limma", "head")
save_plot("tmp/test.png",p2 , base_width = 7, base_height = 7, ncol = 3)

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
save_plot("tmp/DE_volcano_limma-vs-vicar.png",panel , base_width = 7, base_height = 7)




x = results$vicar |>
    dlply("tissue") |> 
    map(\(x) slice_min(x, lfdr, n = 500)) |> 
    map(\(x) mutate(x, n = 1:nrow(x))) |>
    list_rbind() |>
    select(n, Geneid, SYMBOL, tissue, lfdr)

p = ggplot(x, aes(n, y = lfdr)) + 
    geom_point() + 
    facet_wrap(~tissue, scales = "free_y") +
    theme_bw() #+ coord_cartesian(ylim = c(0, 1))
save_plot("tmp/test.png",p , base_width = 7, base_height = 7, ncol = 2)
