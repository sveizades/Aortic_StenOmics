# Supplemental Figure 9d: UMAP of endothelial subclusters coloured by level-3 annotation from SALTIRE3 scRNA-seq
library(Seurat)
library(ggplot2)

source("../../utils.R")

# Load data
atlas_SALTIRE3_sc_endo <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/04_S3_scRNAseq_endothelial_annotation/endothelial.rds")

atlas_endo_umap_df <- data.frame(UMAP1 = atlas_SALTIRE3_sc_endo@reductions$umap@cell.embeddings[, "umap_1"], UMAP2 = atlas_SALTIRE3_sc_endo@reductions$umap@cell.embeddings[, "umap_2"], annotations_level3 = atlas_SALTIRE3_sc_endo$annotations_level3_readable, row.names = colnames(atlas_SALTIRE3_sc_endo))
plt <- ggplot(atlas_endo_umap_df, aes(x = UMAP1, y = UMAP2, color = annotations_level3)) +
  ggrastr::rasterise(geom_point(size = 0.2, stroke = 0, shape = 16), dpi = 600) +
  scale_color_manual(values = A3_celltype_cols) +
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
ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_9d.svg"), width = 4, height = 4, units = "cm")
write.csv(atlas_endo_umap_df, file = paste0(source_data_output_path, "Supplemental_Figure_9d.csv"), row.names = TRUE)
