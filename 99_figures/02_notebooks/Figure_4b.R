# Figure 4b: SALTIRE3 scRNA-seq VIC sub-state UMAP

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
source("../../utils.R")

atlas_SALTIRE3_sc <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/merged_objects.rds")

paga_coords <- read.csv("../../03_SALTIRE3_scRNA_seq/01_data/merged_objects_paga_umap_coordinates_annotations_level3_readable.csv",row.names = 1,check.names = FALSE)

umap_df <- paga_coords %>%
  rownames_to_column("cell_id") %>%
  left_join(
    atlas_SALTIRE3_sc@meta.data %>%
      rownames_to_column("cell_id") %>%
      select(cell_id, annotations_level2_readable),
    by = "cell_id"
  )

plt <-ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = annotations_level2_readable)) +
  ggrastr::rasterise(geom_point(size = 0.2, stroke = 0, shape = 16), dpi = 600) +
  scale_color_manual(values = A2_celltype_cols) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position = 'none')+
  coord_equal()
ggsave(plot = plt, filename = paste0(main_figure_output_path, "Figure_4b.svg"), width = 6, height = 6, units = "cm")
write.csv(umap_df, file = paste0(source_data_output_path, "Figure_4b.csv"), row.names = TRUE)
