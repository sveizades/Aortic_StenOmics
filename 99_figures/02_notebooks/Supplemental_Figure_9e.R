# Supplemental Figure 9e: Dot plot of endothelial subtype marker genes across endothelial subclusters from SALTIRE3 scRNA-seq
library(Seurat)
library(ggplot2)
library(scico)

source("../../utils.R")

# Load data
endothelial <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/04_S3_scRNAseq_endothelial_annotation/endothelial.rds")

endothelial$annotations_level3_readable <- factor(endothelial$annotations_level3_readable, levels = rev(c("VEC-venous", "VEC-arterial", "VEC-capillary", "VEC-lymphatic", "VEC-activated", "VEC-endocardial")))

endothelial_markers <- c(
  # Venous
  "ACKR1", "COL15A1", "PLVAP",
  # Arterial
  "SEMA3G", "HEY1",
  # Capillary
  "MAFB", "RGCC", "PPARG",
  # Lymphatic
  "PROX1", "PDPN", "CCL21",
  # Activated
  "SELE", "ICAM1", "CCL2", "IL6",
  # Endocardial
  "BMPER", "GATA4", "CDH11"
)
plt <- DotPlot(endothelial, features = endothelial_markers, group.by = "annotations_level3_readable", col.min = 0, col.max = 1.5, dot.scale = 2) +
  scale_colour_gradientn(colours = scico(10, palette = "oslo", direction = -1)) +
  coord_equal() +
  scale_y_discrete(labels = c("VEC-venous" = "Vascular-venous", "VEC-arterial" = "Vascular-arterial", "VEC-capillary" = "Vascular-capillary", "VEC-endocardial" = "Endocardial", "VEC-lymphatic" = "Vascular-lymphatic", "VEC-activated" = "Vascular-activated")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 1)),
        axis.text.y = element_text(size = 5, colour = "black", margin = margin(r = 1)),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"), legend.position = "none")

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_9e.svg"), width = 6, height = 3, units = "cm")
write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_9e.csv"), row.names = TRUE)
