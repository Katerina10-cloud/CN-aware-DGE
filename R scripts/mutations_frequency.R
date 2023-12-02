#Mutations & CNV frequency analysis CLL (public data Nature Cancer)

library(tidyverse)
library(ggplot2)
library(CNAgc)

load("./gene.RData")

colnames(gene)
n_samples = nrow(gene %>% dplyr::filter(!is.na(del13q),!is.na(trisomy12),!is.na(IGHV)))

table_n =
  gene %>%
  dplyr::filter(!is.na(del13q),!is.na(trisomy12),!is.na(IGHV)) %>%
  group_by(del13q, trisomy12, IGHV) %>%
  reframe(n = n()) %>%
  mutate(
    group = case_when(
      del13q == 0 & trisomy12 == 0 ~ "WT",
      del13q == 1 & trisomy12 == 0 ~ "del13q",
      del13q == 0 & trisomy12 == 1 ~ "trisomy12",
      del13q == 1 & trisomy12 == 1 ~ "both"
    )
  ) %>%
  group_by(group) %>%
  reframe(N = sum(n), group, n, IGHV)

plot_n = table_n %>%
  # tidyr::pivot_longer(cols = c(IGHV), values_to = “nn”)
  ggplot() +
  geom_bar(aes(
    y = group,
    x = n,
    fill = factor(IGHV, levels = c(0, 1))
  ), stat = "identity") +
  CNAqc:::my_ggplot_theme() +
  scale_fill_manual(
    values = c("deepskyblue1", "goldenrod1"),
    labels = c("IGHV-UM", "IGHV-M")
  ) +
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.y =  element_blank(),
        axis.title.y = element_blank())
table_genotype = gene %>%
  dplyr::as_tibble(rownames = "sample") %>%
  mutate(
    group = case_when(
      del13q == 0 & trisomy12 == 0 ~ "WT",
      del13q == 1 & trisomy12 == 0 ~ "del13q",
      del13q == 0 & trisomy12 == 1 ~ "trisomy12",
      del13q == 1 & trisomy12 == 1 ~ "both"
    )
  ) %>%
  dplyr::filter(!is.na(del13q),!is.na(trisomy12)) %>%
  select(-c(del13q, trisomy12)) %>%
  reshape2::melt(id = c('sample', 'group')) %>%
  as_tibble() %>%
  # tidyr::pivot_longer(cols = colnames(.), names_to = “genotype”, values_to = “n”, values_drop_na = T) %>%
  group_by(group, variable) %>%
  reframe(n = sum(value, na.rm = T)) %>%
  unique() %>%
  full_join(table_n %>% select(group, N)) %>%
  mutate(freq = n/N) %>%
  unique()
  # select(group, variable, freq)
  # tidyr::pivot_wider(names_from = variable, values_from = freq, values_fill = 0)

plot_genotype = table_genotype %>%
  group_by(variable) %>%
  reframe(freq_var = sum(freq), n_var = sum(n), group, n, N, freq) %>%
  filter(freq_var >= .15 | n_var >= 10) %>%
  ggplot()+
  geom_tile(aes(y = group, x = variable, fill = freq))+
  scale_fill_distiller(palette= 'Spectral')+
  CNAqc:::my_ggplot_theme()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill = guide_colorbar(barwidth = unit(2.5, 'cm')))+
  xlab("Mutation")

library(patchwork)
plot_final = plot_genotype+plot_n+plot_layout(ncol = 2, guides = 'collect', widths = c(4.5,1))+
  plot_annotation(title = "Data from Nat. Canc.",
                  subtitle = paste0("Samples: ", n_samples)) & theme(legend.position  = 'bottom')


