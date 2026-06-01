# Supplemental Figure 9i: Dot plot of myeloid subtype marker genes across myeloid subclusters from SALTIRE3 scRNA-seq
library(Seurat)
library(ggplot2)
library(scico)

source("../../utils.R")

# Load data
mp_cells <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/myeloid_cells.rds")

mp_cells$annotations_level3_readable <- factor(mp_cells$annotations_level3_readable, levels = rev(c("Resident macrophage",
                                                                                                    "Resident activated macrophage",
                                                                                                    "Scar-associated macrophage",
                                                                                                    "Osteoclast",
                                                                                                    "IFN-stimulated macrophage",
                                                                                                    "Classical monocyte",
                                                                                                    "Non-classical monocyte",
                                                                                                    "cDC1",
                                                                                                    "cDC2",
                                                                                                    "moDC",
                                                                                                    "pDC")))

myeloid_markers <- c(
  # Resident macrophages
  "C1QA", "C1QB", "C1QC", "FOLR2", "LYVE1",

  # Activated macrophages
  "EGR1", "EGR2", "CCL3", "CCL4",

  # SPP1 macrophages
  "LPL", "SPP1", "FN1", "APOE", "CD109",

  # Osteoclast-like
  "ACP5", "CTSK", "MMP9",

  # IFN macrophages
  "IFIT1", "IFI44L",


  # Monocytes
  "S100A12", "SERPINB2", "EREG", "AREG",
  "FCN1", "CCR2",  "FCGR3A",
  "CLEC9A", "XCR1",
  "CD1C", "CD1B",
  "CCL19", "LAMP3",

  # DCs
  "GZMB", "JCHAIN", "CLEC4C"


)
plt <- DotPlot(mp_cells, features = myeloid_markers, group.by = "annotations_level3_readable", col.min = 0, col.max = 1.5, dot.scale = 2) +
  scale_colour_gradientn(colours = scico(10, palette = "oslo", direction = -1)) +
  coord_equal() +
  #scale_y_discrete(labels = c("MP-macro-res" = "Resident macrophage", "MP-macro-res-activated" = "Resident activated macrophage", "MP-macro-SA" = "Scar-associated macrophage","MP-macro-osteoclast" = "Osteoclast","MP-macro-IFN" = "IFN-stim macrophage","MP-mono-classical" = "Classical monocyte","MP-mono-nonclassical" = "Non-classical monocyte", "MP-DC-c1" = "cDC1", "MP-DC-c2" = "cDC2", "MP-DC-m" = "moDC", "MP-DC-p" = "pDC"))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 1)),
        axis.text.y = element_text(size = 5, colour = "black", margin = margin(r = 1)),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"), legend.position = "none")

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_9i.svg"), width = 12, height = 4, units = "cm")
write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_9i.csv"), row.names = TRUE)
