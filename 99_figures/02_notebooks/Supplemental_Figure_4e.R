# Supplemental Figure 4e: UMAP of VIC MAGIC-imputed CNN1 expression with Contractile cell density contour (UBC snRNA-seq).
library(dplyr)
library(ggplot2)
library(Seurat)
library(scico)

source("../../utils.R")

vic_seurat <- readRDS("../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/vics_magic.rds")

vic_df <- data.frame(UMAP1 = vic_seurat@reductions$umap@cell.embeddings[, "umap_1"], UMAP2 = vic_seurat@reductions$umap@cell.embeddings[, "umap_2"], annotations_level2_readable = vic_seurat$annotations_level2_readable, CNN1 = GetAssayData(vic_seurat, assay = "MAGIC", layer = "data")["CNN1", ], row.names = colnames(vic_seurat))
vic_df <- vic_df %>%
  arrange(CNN1)
cluster_cellsVICn5 <- vic_df[vic_df$annotations_level2_readable == "Contractile", ]

plt <- ggplot(vic_df, aes(x = UMAP1, y = UMAP2, color = CNN1)) +
  ggrastr::rasterise(geom_point(size = 0.15, stroke = 0, shape = 16), dpi = 600) +
  scale_color_scico(palette = "batlow", direction = -1) +
  stat_density_2d(
    data = cluster_cellsVICn5,
    aes(x = UMAP1, y = UMAP2),
    colour = "black",
    linetype = "dashed",
    linewidth = 0.6,
    breaks = 0.02,
    inherit.aes = FALSE
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    legend.position = "none"
  ) +
  coord_equal()

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_4e.svg"), width = 4, height = 4, units = "cm")

write.csv(vic_df, file = paste0(source_data_output_path, "Supplemental_Figure_4e.csv"), row.names = TRUE)
