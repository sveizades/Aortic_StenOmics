# Figure 4e: Xenium gene expression by distance to nodule edge

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

source("../../utils.R")
atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")

genes_of_interest <- c("FAP", "COL11A1", "IBSP", "CHAD", "BGLAP", "COL10A1")  # start small
gene_cols <- c(
  "FAP"      = "#D81B60", 
  "COL11A1"  = "#5E60CE",  
  "IBSP"     = "#1E88E5",  
  "CHAD"     = "#7CB342",  
  "BGLAP"    = "#F9A825",  
  "COL10A1"  = "#8E24AA"   
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
  mutate(expr_scaled = as.numeric(scale(mean_expr))) %>%   # z-score within gene
  ungroup()

plt <- ggplot(plot_df, aes(dist_mid, expr_scaled, colour = gene))+
  scale_colour_manual(values = gene_cols)  +
  geom_smooth(se = FALSE, method = "loess", span = 0.75, linewidth = 0.25) +
  scale_x_continuous(labels = comma, expand = c(0, 0)) +
  ggpubr::theme_pubr() +
  theme(
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none", axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),axis.ticks.length = unit(.1, "cm"),
    axis.text.x = element_text(colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    plot.margin = unit(c(1, 3, 1, 1), "mm"))
ggsave(plot = plt, 
       filename = paste0(main_figure_output_path, "Figure_4e.svg"), 
       width = 2.5, 
       height = 2.5, 
       units = "cm")
write.csv(plot_df, file = paste0(source_data_output_path, "Figure_4e.csv"), row.names = TRUE)


