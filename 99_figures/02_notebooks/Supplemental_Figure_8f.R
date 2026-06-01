# Supplemental Figure 8f: Dot plot of level-1 cell type marker genes across annotated SALTIRE3 scRNA-seq atlas cell types
library(Seurat)
library(ggplot2)
library(scico)

source("../../utils.R")

# Load data
atlas_SALTIRE3_sc <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/merged_objects.rds")

SALTIRE3_level1_atlas_markers <- c("DCN", "COL1A2", "LUM", "PCOLCE", "PRELP", #VICs
                       "CD3G", "CD2", "CD3D", "IL7R", "CD52", #T
                       "NKG7", "KLRD1", "KLRB1", "PRF1", "CTSW", #NK
                       "CD163", "FCER1G", "LYZ", "FCGR2A", "C1QC", #Myeloid
                       "TPSAB1", "CPA3", "KIT", "TPSB2", "CLEC4OP", #Granulocytes
                       "VWF", "CLEC14A", "PLVAP", "CDH5", "RAMP2", #endo
                       "CD79A", "MS4A1", "BANK1", "RALGPS2", "LY9", #B
                       "MKI67", "TOP2A", "PCLAF", "RRM2", "CDK1") #Cycling
atlas_SALTIRE3_sc$annotations_level1_readable <- factor(atlas_SALTIRE3_sc$annotations_level1_readable, levels = c("Proliferating", "B cell", "Endothelial cell", "Granulocyte", "Myeloid", "NK cell", "T cell", "Valvular interstitial cell"))
plt <- DotPlot(atlas_SALTIRE3_sc, features = SALTIRE3_level1_atlas_markers, group.by = "annotations_level1_readable", col.min = 0, col.max = 1.5, dot.scale = 2) +
  scale_colour_gradientn(colours = scico(10, palette = "oslo", direction = -1)) +
  coord_equal() +
  scale_y_discrete(labels = c("Valvular interstitial cell" = "VICs", "T cell" = "T cells", "NK cell" = "NK cells", "Myeloid" = "Myeloid", "Granulocyte" = "Granulocytes", "Endothelial cell" = "Endothelial cells", "Proliferating" = "Proliferating", "B cell" = "B cells")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 1)),
        axis.text.y = element_text(size = 5, colour = "black", margin = margin(r = 1)),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"), legend.position = "none")

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_8f.svg"), width = 10, height = 3, units = "cm")
write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_8f.csv"), row.names = TRUE)
