# Supplemental Figure 9g: Faceted dot plot of T/NK subtype marker genes across T/NK subclusters from SALTIRE3 scRNA-seq
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(scico)

source("../../utils.R")

# Load data
tnk_cells <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/tnk_cells.rds")

tnk_markers <- c(
  # CD3
  "CD3D", "CD3E", "CD3G",
  # TCR
  #"TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2",
  # CD4
  "CD4", "CD40LG",
  # Naive
  "CCR7", "SELL", #"TCF7", "LEF1", "MAL",
  # Activation / EM
  "HLA-DRA", "CD69", "IFNG", "TNF", "GZMK",
  # Th17
  "CCR6", "RORC", "IL23R",
  # Treg
  "FOXP3", "IL2RA", "CTLA4", #"TIGIT",
  # CD8
  "CD8A", "CD8B",
  # Cytotoxic
  "NKG7", "GNLY", "PRF1", "GZMB", "GZMH",
  # MAIT
  "SLC4A10", "TRAV1-2", "KLRB1",
  # γδ T
  "TRDV2", "TRGV9",
  # NK
  "KLRD1", "FCGR3A", "NCAM1"
)

marker_groups <- tibble(
  features.plot = tnk_markers,
  gene_group = c(
    rep("CD3", 3),
    #rep("TCR", 6),
    rep("CD4", 2),
    rep("Naive", 2),
    rep("Activation / EM", 5),
    rep("Th17", 3),
    rep("Treg", 3),
    rep("CD8", 2),
    rep("Cytotoxic", 5),
    rep("MAIT", 3),
    rep("γδ T", 2),
    rep("NK", 3)
  )
)

# build the Seurat dotplot first
p <- DotPlot(
  tnk_cells,
  features = tnk_markers,
  group.by = "annotations_level3_readable",
  col.min = 0,
  col.max = 1.5,
  dot.scale = 2
)

# add facet information to the plotting data
plot_data <- p$data %>%
  left_join(marker_groups, by = "features.plot")

# preserve your desired gene order
plot_data$features.plot <- factor(plot_data$features.plot, levels = tnk_markers)

# preserve facet order
plot_data$gene_group <- factor(
  plot_data$gene_group,
  levels = c("CD3",
             #"TCR",
             "CD4", "Naive", "Activation / EM", "Th17", "Treg", "CD8", "Cytotoxic", "MAIT", "γδ T", "NK")
)

plot_data$id <- factor(
  plot_data$id,
  levels = rev(c("T-CD4-naive", "T-CD4-EM", "T-CD4-cyto", "T-CD4-Th1", "T-CD4-Th17", "T-CD4-Treg", "T-CD8-cyto", "T-CD8-EM", "MAIT", "T-CD8-MAIT", "T-gd-cyto", "T-gd-V9", "NK-CD56-bright", "NK-CD56-dim"))
)

plt <- ggplot(
  plot_data,
  aes(x = features.plot, y = id, size = pct.exp, colour = avg.exp.scaled)
) +
  geom_point() +
  facet_grid(~ gene_group, scales = "free_x", space = "free_x", switch = "x") +
  scale_colour_gradientn(colours = scico(10, palette = "oslo", direction = -1)) +
  scale_size(range = c(0, 2)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 5, colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 1)),
    axis.text.y = element_text(size = 5, colour = "black", margin = margin(r = 1)),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"),
    legend.position = "none",
    strip.placement = "none",
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "mm")
  )
ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_9g.svg"), width = 12, height = 4, units = "cm")
write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_9g.csv"), row.names = TRUE)
