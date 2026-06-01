# Figure 1h: UMAP of snRNA-seq atlas coloured by bulk RNA-seq disease score projection (UBC cohort).
library(Seurat)
library(ggplot2)
library(scico)

source("../../utils.R")

# Load data
atlas <- readRDS("../../02_UBC_snRNA_seq/03_output/11_UBC_snRNAseq_annotation_combining/atlas.rds")

atlas_umap_df <- data.frame(UMAP1 = atlas@reductions$umap@cell.embeddings[, "umap_1"], UMAP2 = atlas@reductions$umap@cell.embeddings[, "umap_2"], bulk_score = atlas$Cluster1, row.names = colnames(atlas))

plt <- ggplot(atlas_umap_df, aes(x = UMAP1, y = UMAP2, color = bulk_score)) +
  ggrastr::rasterise(geom_point(size = 0.2, stroke = 0, shape = 16), dpi = 600) +
  scale_color_scico(palette = "vik", direction = 1, midpoint = 0, limits = c(min(atlas_umap_df$bulk_score), quantile(atlas_umap_df$bulk_score, 0.99))) +
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

write.csv(atlas_umap_df, file = paste0(source_data_output_path, "Figure_1h.csv"), row.names = TRUE)

ggsave(plot = plt, filename = paste0(main_figure_output_path, "Figure_1h.svg"), width = 6, height = 6, units = "cm")
