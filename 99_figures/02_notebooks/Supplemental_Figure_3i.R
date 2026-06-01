# Supplemental Figure 3i: Dot plot of canonical marker gene expression across UBC snRNA-seq atlas cell types.
library(Seurat)
library(ggplot2)
library(scico)

source("../../utils.R")

atlas <- readRDS("../../02_UBC_snRNA_seq/03_output/11_UBC_snRNAseq_annotation_combining/atlas.rds")

atlas$annotations_level1_readable <- factor(atlas$annotations_level1_readable, levels = c("VICs", "Myeloid", "B cells", "T cells", "NK cells", "Endothelial cells", "Adipocytes", "Mast cells", "Proliferating"))

ubc_atlas_markers <- c("CDH19", "DCN", "NLGN4X", "PDGFRA", "FBLN5", #VICs
                       "CD163", "MRC1", "TGFBI", "HLA-DRA", "C1QC", #MPs
                       "BLK", "MZB1", "MS4A1", "JCHAIN", "IGKC", #B
                       "THEMIS", "CD247", "CD2", "IL7R", "CD3E", #T
                       "KLRC1", "KLRC3", "NKG7", "PRF1", "GZMA", #NK
                       "TEK", "PODXL", "FLT1", "CDH5", "VWF", #endo
                       "LPL", "PLIN1", "ADIPOQ", "FABP4", "PLIN4", #adipocytes
                       "CPA3", "KIT", "HDC", "TPSB2", "TPSAB1",
                       "POLQ", "MELK", "KNL1", "ANLN", "MKI67")

plt <- DotPlot(atlas, features = ubc_atlas_markers, group.by = "annotations_level1_readable", col.min = 0, col.max = 1.5, dot.scale = 2) +
  scale_colour_gradientn(colours = scico(10, palette = "oslo", direction = -1)) +
  coord_equal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 1)),
        axis.text.y = element_text(size = 5, colour = "black", margin = margin(r = 1)),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"), legend.position = "none")

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_3i.svg"), width = 11, height = 3, units = "cm")

write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_3i.csv"), row.names = TRUE)
