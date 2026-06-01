# Supplemental Figure 9c: Dot plot of VIC subtype marker genes across VIC subclusters from SALTIRE3 scRNA-seq
library(Seurat)
library(ggplot2)
library(scico)

source("../../utils.R")

# Load data
vics <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/03_S3_scRNAseq_vic_annotation/vics.rds")

vics$annotations_level3_readable <- factor(vics$annotations_level3_readable, levels = rev(c("VIC-quiescent", "VIC-transitional", "VIC-FAP-osteogenic", "VIC-contractile", "VIC-spongiosa", "VIC-neural-crest", "VIC-IFN")))

sc_vic_markers <- c("RARRES2", "APOE", "ADH1B", #Quiescent
                    "NTM", "MXRA5", "HMCN1", #Quiescent
                    "OMD", "COMP", "FAP", "BGLAP", "IBSP", "ENPP1", #FAP+ osteogenic
                    "MYH11", "MCAM", "CNN1", #Contractile
                    "CCN5", "GFRA1", "SFRP1", "GFAP", "PLA2G2A",  #Spongiosa
                    "TOX", "SOX10", "GRIK2", "CHL1", #Neural crest-like
                    "IFI44L", "IFIT3", "MX1") #IFN-stimulated

plt <- DotPlot(vics, features = sc_vic_markers, group.by = "annotations_level3_readable", col.min = 0, col.max = 1.5, dot.scale = 2) +
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
        panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"), legend.position = "none")

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_9c.svg"), width = 7.5, height = 3, units = "cm")
write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_9c.csv"), row.names = TRUE)
