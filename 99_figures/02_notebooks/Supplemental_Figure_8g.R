# Supplemental Figure 8g: Stacked bar chart of cell type proportions per donor in SALTIRE3 scRNA-seq atlas
library(Seurat)
library(ggplot2)

source("../../utils.R")

# Load data
atlas_SALTIRE3_sc <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/merged_objects.rds")

atlas_SALTIRE3_sc$donor <- sapply(strsplit(atlas_SALTIRE3_sc$donor_id, "_"), `[`, 2)
atlas_SALTIRE3_sc$donor <- factor(atlas_SALTIRE3_sc$donor, levels = sort(unique(atlas_SALTIRE3_sc$donor)))

sample_proprotions <- as.data.frame(table(atlas_SALTIRE3_sc$donor, atlas_SALTIRE3_sc$annotations_level1_readable)); colnames(sample_proprotions) <- c("Sample", "Cluster", "Proportion")

p <- ggplot(sample_proprotions, aes(fill = Cluster, y = Proportion, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values = A1_celltype_cols) +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "line"),
        plot.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, margin = margin(t = -2)),
        axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + NoLegend()
ggsave(plot = p, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_8g.svg"), width = 3, height = 3, units = "cm")
write.csv(sample_proprotions, file = paste0(source_data_output_path, "Supplemental_Figure_8g.csv"), row.names = TRUE)
