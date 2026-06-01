# Figure 2b: Dot plot of marker gene expression across VIC subpopulations (UBC snRNA-seq).
library(Seurat)
library(ggplot2)
library(scico)

source("../../utils.R")

# Load data
vics <- readRDS("../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/vics.rds")

vics$annotations_level2_readable <- factor(vics$annotations_level2_readable, levels = rev(c("Transitional", "Quiescent", "Spongiosa", "FAP+ osteogenic", "Contractile", "Neural crest-like", "IFN-stimulated")))
vicn_markers <- c("ELN", "COL1A1", "COL1A2", "COL3A1", "NTM", "CADM2", "KAZN", "GFAP", "FAP", "COMP", "ENPP1", "CRTAC1", "MYH11", "CNN1", "ACTA2", "SOX10", "MX1")
plt <- DotPlot(vics, features = vicn_markers, group.by = "annotations_level2_readable", col.min = 0, col.max = 1.5, dot.scale = 2) +
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

ggsave(plot = plt, filename = paste0(main_figure_output_path, "Figure_2b.svg"), width = 6, height = 6, units = "cm")

plt_size_legend <- DotPlot(
  vics,
  features = vicn_markers,
  group.by = "annotations_level2_readable",
  col.min = 0,
  col.max = 1.5,
  dot.scale = 2
) +
  scale_colour_gradientn(colours = scico(10, palette = "oslo", direction = -1), guide = "none") +
  guides(
    size = guide_legend(
      title = "Percent expressed",
      override.aes = list(colour = "black")
    )
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6, colour = "black"),
    legend.text = element_text(size = 5, colour = "black"),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm")
  )

# Extract just the legend
size_legend <- cowplot::get_legend(plt_size_legend)

# Save legend as standalone SVG
ggsave(
  filename = paste0(main_figure_output_path, "Figure_2b_size_legend.svg"),
  plot = size_legend,
  width = 2.2,
  height = 3,
  units = "cm"
)

write.csv(plt$data, file = paste0(source_data_output_path, "Figure_2b.csv"), row.names = TRUE)
