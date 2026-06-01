# Figure 1g: UMAP of snRNA-seq atlas coloured by cell type annotation (UBC cohort).
library(Seurat)
library(ggplot2)

source("../../utils.R")

# Load data
atlas <- readRDS("../../02_UBC_snRNA_seq/03_output/11_UBC_snRNAseq_annotation_combining/atlas.rds")

atlas_umap_df <- data.frame(UMAP1 = atlas@reductions$umap@cell.embeddings[, "umap_1"], UMAP2 = atlas@reductions$umap@cell.embeddings[, "umap_2"], annotations_level1_readable = atlas$annotations_level1_readable, row.names = colnames(atlas))

plt <- ggplot(atlas_umap_df, aes(x = UMAP1, y = UMAP2, color = annotations_level1_readable)) +
  ggrastr::rasterise(geom_point(size = 0.2, stroke = 0, shape = 16), dpi = 600) +
  scale_color_manual(values = UBC_atlas_colours) +
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

write.csv(atlas_umap_df, file = paste0(source_data_output_path, "Figure_1g.csv"), row.names = TRUE)

ggsave(plot = plt, filename = paste0(main_figure_output_path, "Figure_1g.svg"), width = 6, height = 6, units = "cm")
