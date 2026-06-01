# Supplemental Figure 10i: Contractile VIC marker expression by distance to nodule edge (LOESS smoothed)

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

source("../../utils.R")

atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")

genes_of_interest <- c("ACTA2", "MYH11", "ACTB", "CNN1")
gene_cols <- c(
  "ACTA2" = "#D55E00",  # vermillion
  "MYH11" = "#0072B2",  # blue
  "ACTB"  = "#009E73",  # green
  "CNN1"  = "#CC79A7"
)
expr <- GetAssayData(atlas_lim_nodular_tissue, layer = "data")[genes_of_interest, ]

df <- data.frame(
  cell_id = colnames(atlas_lim_nodular_tissue),
  dist    = atlas_lim_nodular_tissue$dist_to_nodule_edge,
  celltype = atlas_lim_nodular_tissue$annotations_level1,
  t(expr)
) %>%
  tidyr::pivot_longer(
    cols = all_of(genes_of_interest),
    names_to = "gene",
    values_to = "expr"
  ) %>%
  filter(!is.na(dist), !is.na(expr))

df <- df %>%
  mutate(
    dist_bin = cut(
      dist,
      breaks = seq(0, 850, by = 50),
      include.lowest = TRUE
    )
  )

expr_bin_ct <- df %>%
  group_by(gene, celltype, dist_bin) %>%
  summarise(
    mean_expr = mean(expr),
    n_cells   = n(),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 10)

expr_bin_ct <- expr_bin_ct %>%
  mutate(
    dist_mid = as.numeric(str_extract(dist_bin, "\\d+")) + 25
  )

plot_df <- expr_bin_ct %>%
  filter(celltype == "Valvular interstitial cell", gene %in% genes_of_interest) %>%
  group_by(gene) %>%
  mutate(expr_scaled = as.numeric(scale(mean_expr))) %>%
  ungroup()

plt <- ggplot(plot_df, aes(dist_mid, expr_scaled, colour = gene)) +
  scale_colour_manual(values = gene_cols) +
  geom_smooth(se = FALSE, method = "loess", span = 0.75, linewidth = 0.25) +
  scale_x_continuous(labels = comma, expand = c(0, 0)) +
  ggpubr::theme_pubr() +
  theme(
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(.1, "cm"),
    axis.text.x = element_text(colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    plot.margin = unit(c(1, 3, 1, 1), "mm")
  )

ggsave(plot = plt,
       filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10i.svg"),
       width = 2.5,
       height = 2.5,
       units = "cm")

write.csv(plot_df, file = paste0(source_data_output_path, "Supplemental_Figure_10i.csv"), row.names = TRUE)
