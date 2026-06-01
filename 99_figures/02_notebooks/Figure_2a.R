# Figure 2a: UMAP of VIC subpopulations coloured by level-2 annotation (UBC snRNA-seq).
library(Seurat)
library(ggplot2)

source("../../utils.R")

# Load data
vics <- readRDS("../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/vics.rds")

atlas_umap_df <- data.frame(UMAP1 = vics@reductions$umap@cell.embeddings[, "umap_1"], UMAP2 = vics@reductions$umap@cell.embeddings[, "umap_2"], annotations_level2_readable = vics$annotations_level2_readable, row.names = colnames(vics))

plt <- ggplot(atlas_umap_df, aes(x = UMAP1, y = UMAP2, color = annotations_level2_readable)) +
  ggrastr::rasterise(geom_point(size = 0.4, stroke = 0, shape = 16), dpi = 600) +
  scale_color_manual(values = UBC_atlas_colours_vics) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        legend.position = "none") +
  coord_equal()

ggsave(plot = plt, filename = paste0(main_figure_output_path, "Figure_2a.svg"), width = 6, height = 6, units = "cm")
write.csv(atlas_umap_df, file = paste0(source_data_output_path, "Figure_2a.csv"), row.names = TRUE)
