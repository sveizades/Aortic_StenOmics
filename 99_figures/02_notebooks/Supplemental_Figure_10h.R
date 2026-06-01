# Supplemental Figure 10h: Distance to nodule edge by VIC subtype (Xenium spatial analysis)

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(rstatix)

source("../../utils.R")

atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")

df <- data.frame(
  cell_id = colnames(atlas_lim_nodular_tissue),
  distance_to_nodule_edge = atlas_lim_nodular_tissue@meta.data[["dist_to_nodule_edge"]],
  annotations_level2 = atlas_lim_nodular_tissue$annotations_level2
) %>%
  filter(
    !is.na(distance_to_nodule_edge),
    !is.na(annotations_level2),
    str_starts(annotations_level2, "VIC")
  )

# order by median distance
vic_order <- df %>%
  group_by(annotations_level2) %>%
  summarise(median_dist = median(distance_to_nodule_edge), .groups = "drop") %>%
  arrange(median_dist) %>%
  pull(annotations_level2)

df <- df %>%
  mutate(annotations_level2 = factor(annotations_level2, levels = vic_order))

# global test
kw_stats <- df %>%
  kruskal_test(distance_to_nodule_edge ~ annotations_level2)

# pairwise only vs VIC-osteogenic
pairwise_vs_osteo <- df %>%
  pairwise_wilcox_test(
    distance_to_nodule_edge ~ annotations_level2,
    comparisons = list(c("VIC-FAP-osteogenic", "VIC-contractile")),
    p.adjust.method = "BH"
  ) %>%
  mutate(
    p_label = ifelse(
      p.adj < 0.001,
      "p<0.001",
      paste0("p=", signif(p.adj, 2))
    )
  ) %>%
  add_xy_position(x = "annotations_level2")

pairwise_vs_osteo <- pairwise_vs_osteo %>%
  mutate(y.position = 2500)

write.csv(kw_stats, file = paste0(source_data_output_path, "Supplemental_Figure_10h_distance_to_nodule_edge_VIC_kw_stats.csv"), row.names = FALSE)

plt <- ggplot(df, aes(x = annotations_level2, y = distance_to_nodule_edge, fill = annotations_level2)) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    linewidth = 0.3,
    colour = "#3A3A3A"
  ) +
  stat_pvalue_manual(
    pairwise_vs_osteo,
    label = "p_label",
    tip.length = 0.01,
    size = 2
  ) +
  scale_x_discrete(labels = VIC_celltype_names) +
  scale_y_continuous(limits = c(0, 3000)) +
  scale_fill_manual(values = vic_cols_black_background, drop = FALSE) +
  ggpubr::theme_pubr() +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Distance to nodule edge"
  ) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text.x = element_text(
      colour = "black",
      size = 5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(colour = "black", size = 5),
    axis.title.y = element_text(size = 6),
    plot.margin = unit(c(2, 3, 1, 1), "mm")
  )

ggsave(
  plot = plt,
  filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10h.svg"),
  width = 5,
  height = 4,
  units = "cm"
)
write.csv(df, file = paste0(source_data_output_path, "Supplemental_Figure_10h.csv"), row.names = TRUE)
