# Supplemental Figure 3j: VIC signature score boxplots across UBC snRNA-seq atlas cell types with pairwise significance.
library(Seurat)
library(ggplot2)
library(scico)
library(scales)
library(ggpubr)
library(dplyr)
library(tidyr)

source("../../utils.R")

atlas <- readRDS("../../02_UBC_snRNA_seq/03_output/11_UBC_snRNAseq_annotation_combining/atlas.rds")

atlas$annotations_level1 <- factor(atlas$annotations_level1, levels = c("vics", "myeloid", "B_cell", "T_cell", "NK_cell", "endothelial", "adipocytes", "mast", "proliferating"))
x_order <- levels(atlas@meta.data$annotations_level1)

pvals_df <- read.csv("../../02_UBC_snRNA_seq/03_output/11_UBC_snRNAseq_annotation_combining/Supplementary_Table_CellType_Signature_Pairwise.csv") %>%
  rename(pvalue_label = p_label)

pvals_vics <- pvals_df %>%
  filter(group1 == "vics" | group2 == "vics")

pvals_vics <- pvals_vics %>%
  mutate(
    other_group = ifelse(group1 == "vics", group2, group1),
    x_position = match(other_group, x_order)
  ) %>%
  arrange(x_position)

y_max <- max(atlas@meta.data$Cluster1, na.rm = TRUE)

pvals_vics <- pvals_vics %>%
  mutate(
    y.position = y_max * (0.08 * (row_number() - 1)) + 0.2
  )

UBC_atlas_colours <- c(vics = "#332288", myeloid = "#CC6677", B_cell = "#DDCC77", T_cell = "#117733", NK_cell = "#88CCEE", endothelial = "#882255", adipocytes = "#44AA99", mast = "#999933", proliferating = "#AA4499")
UBC_atlas_colours_fill <- colorspace::lighten(UBC_atlas_colours, amount = 0.5)

atlas@meta.data$annotations_level1 <- factor(
  atlas@meta.data$annotations_level1,
  levels = names(UBC_atlas_colours)
)

plt <- ggplot(atlas@meta.data, aes(x = annotations_level1, y = Cluster1)) +

  stat_boxplot(
    aes(colour = annotations_level1),
    geom = "errorbar",
    width = 0.5,
    alpha = 1
  ) +

  geom_boxplot(
    aes(fill = annotations_level1, colour = annotations_level1),
    outlier.shape = NA,
    width = 0.6
  ) +

  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.10))) +

  scale_colour_manual(values = UBC_atlas_colours, drop = FALSE) +
  scale_fill_manual(values = UBC_atlas_colours_fill, drop = FALSE) +
  scale_x_discrete(labels = c("vics" = "VICs", "myeloid" = "Myeloid", "B_cell" = "B cells", "T_cell" = "T cells", "NK_cell" = "NK cells", "endothelial" = "Endothelial cells", "adipocytes" = "Adipocytes", "mast" = "Mast cells")) +

  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, colour = "black"),
    axis.text.y = element_text(colour = "black", size = 5),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 7),
    plot.margin = unit(c(2, 2, 2, 2), "mm"),
    legend.position = "none"
  ) +
  stat_pvalue_manual(
    pvals_vics,
    label = "pvalue_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 1.6
  )

ggsave(
  file.path(supplemental_figure_output_path, "Supplemental_Figure_3j.svg"),
  plt,
  width = 40, height = 60, units = "mm"
)

write.csv(atlas@meta.data[, c("annotations_level1_readable", "Cluster1")], file = paste0(source_data_output_path, "Supplemental_Figure_3i.csv"), row.names = TRUE)
