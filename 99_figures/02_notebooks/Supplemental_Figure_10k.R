# Supplemental Figure 10k: Dot plot of niche marker genes across spatial niches in Xenium atlas

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(scico)
library(rstatix)

source("../../utils.R")

atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")
atlas_lim_nodular_tissue$niche <- factor(atlas_lim_nodular_tissue$niche, levels = c("Immune", "Peri-nodular", "Spongiosa", "Cap"))
Idents(atlas_lim_nodular_tissue) <- "niche"
FAM(atlas_lim_nodular_tissue, "Supplemental_Figure_10k_atlas_lim_nodular_tissue_niche", source_data_output_path)
dotplot_genes <- c(
  # Cap niche
  "MYH11", "ACTA2", "CNN1",
  # Spongiosa niche
  "GFAP", "SOX10", "CXCL14",
  # Osteogenic / mineralisation (peri-nodular)
  "IBSP", "SPP1", "COL11A1", "CHAD", "BGLAP",
  "ENPP1", "RUNX2",
  # Myeloid co-enrichment in peri-nodular
  "CD14", "CD163", "CSF1R", "CD68",
  # Immune niche
  "JCHAIN", "IGKC", "CD2", "CD3D", "VWF"
)
plt <- DotPlot(atlas_lim_nodular_tissue, features = dotplot_genes, group.by = "niche", col.min = 0, col.max = 1.5, dot.scale = 2) +
  scale_colour_gradientn(colours = scico(10, palette = "oslo", direction = -1)) +
  coord_equal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 1)),
        axis.text.y = element_text(size = 5, colour = "black", margin = margin(r = 1)),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"),
        legend.position = "none")

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10k.svg"), width = 7.5, height = 3, units = "cm")
write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_10k.csv"), row.names = TRUE)
