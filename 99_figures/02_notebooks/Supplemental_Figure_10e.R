# Supplemental Figure 10e: Dot plot of VIC subtype marker genes in Xenium atlas

library(Seurat)
library(ggplot2)
library(scico)

source("../../utils.R")

atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")
atlas_lim_nodular_tissue_vics <- subset(atlas_lim_nodular_tissue, annotations_level1 == "Valvular interstitial cell" & annotations_level2 != "cycling")

atlas_lim_nodular_tissue_vics$annotations_level2 <- factor(atlas_lim_nodular_tissue_vics$annotations_level2, levels = rev(c("VIC-FAP-osteogenic", "VIC-transitional", "VIC-contractile", "VIC-quiescent", "VIC-spongiosa", "VIC-neural-crest", "VIC-IFN")))

Xenium_vic_markers <- c("COMP", "COL11A1", "COL10A1", "FAP", "ENPP1", "IBSP", "CHAD", # FAP+ osteogenic
                    "NTM", "ELN", # Transitional
                    "ACTA2", "MYH11", "ACTB", "CNN1", # Contractile
                    "APOE", "DCN", "PDGFRA", # Quiescent
                    "CXCL14", "GFAP", "PLIN2",  # Spongiosa
                    "SOX10",  # Neural crest-like
                    "MX1", "IFIT3") # IFN-stimulated

plt <- DotPlot(atlas_lim_nodular_tissue_vics, features = Xenium_vic_markers, group.by = "annotations_level2", col.min = 0, col.max = 1.5, dot.scale = 2) +
  scale_colour_gradientn(colours = scico(10, palette = "oslo", direction = -1)) +
  coord_equal() +
  scale_y_discrete(labels = c("VIC-transitional" = "Transitional",
                              "VIC-quiescent" = "Quiescent",
                              "VIC-spongiosa" = "Spongiosa",
                              "VIC-FAP-osteogenic" = "FAP+ osteogenic",
                              "VIC-contractile" = "Contractile",
                              "VIC-neural-crest" = "Neural crest-like",
                              "VIC-IFN" = "IFN-stimulated")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 1)),
        axis.text.y = element_text(size = 5, colour = "black", margin = margin(r = 1)),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"),
        legend.position = "none")

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10e.svg"), width = 6, height = 3, units = "cm")
write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_10e.csv"), row.names = TRUE)
