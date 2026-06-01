# Supplemental Figure 3g: Per-donor nuclei count and QC metrics violin plots for UBC snRNA-seq atlas.
library(Seurat)
library(ggplot2)
library(cowplot)

source("../../utils.R")

atlas <- readRDS("../../02_UBC_snRNA_seq/03_output/11_UBC_snRNAseq_annotation_combining/atlas.rds")
atlas_metadata <- atlas@meta.data

ubc_donor_colours_included <- UBC_donor_colours[sort(unique(sapply(strsplit(atlas_metadata$donor_id, "_"), `[`, 2)))]

atlas_metadata$pool_donor <- paste0(sapply(strsplit(atlas_metadata$pool_id, "_"), `[`, 2), "_", sapply(strsplit(atlas_metadata$donor_id, "_"), `[`, 2))
atlas_metadata$pool_donor <- factor(atlas_metadata$pool_donor, levels = sort(unique(atlas_metadata$pool_donor)))
atlas_metadata$donor <- sapply(strsplit(atlas_metadata$donor_id, "_"), `[`, 2)
atlas_metadata$donor <- factor(atlas_metadata$donor, levels = sort(unique(atlas_metadata$donor)))

plt1 <- ggplot(atlas_metadata, aes(x = donor)) +
  geom_bar(fill = ubc_donor_colours_included) +
  theme_minimal(base_size = 14) +
  scale_y_log10(expand = expansion(mult = c(0, 0.05)), breaks = c(10, 100, 1000, 10000)) +
  labs(
    y = "Number of Nuclei"
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(colour = "black", size = 5),
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text = element_text(colour = "black", size = 5),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        legend.position = "none")

plt2 <- ggplot(atlas_metadata[atlas_metadata$nCount_RNA < 20000, ], aes(x = donor, y = nCount_RNA, fill = donor)) +
  geom_violin(linewidth = 0.2, width = 1) +
  scale_fill_manual(values = ubc_donor_colours_included) +
  labs(
    x = "Sample",
    y = "Number of UMI"
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(colour = "black", size = 5),
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text = element_text(colour = "black", size = 5),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = 0))

plt3 <- ggplot(atlas_metadata, aes(x = donor, y = nFeature_RNA, fill = donor)) +
  geom_violin(linewidth = 0.2, width = 1) +
  scale_fill_manual(values = ubc_donor_colours_included) +
  labs(
    x = "Sample",
    y = "Number of genes"
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_text(colour = "black", size = 5),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 5),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = 0))

aligned <- align_plots(plt1, plt2, plt3, align = "v", axis = "l")

ggsave(plot = plot_grid(plotlist = aligned, ncol = 1, align = "v"), filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_3g.svg"), width = 6, height = 6, units = "cm")

write.csv(atlas_metadata[, c("nFeature_RNA", "nCount_RNA", "donor_id")], file = paste0(source_data_output_path, "Supplemental_Figure_3g.csv"), row.names = TRUE)
