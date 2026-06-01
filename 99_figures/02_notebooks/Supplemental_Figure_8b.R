# Supplemental Figure 8b: Per-donor cell counts, UMI counts, and gene counts for SALTIRE3 scRNA-seq atlas
library(Seurat)
library(ggplot2)
library(cowplot)

source("../../utils.R")

# Load data
atlas_SALTIRE3_sc <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/merged_objects.rds")
df_atlas_metadata <- atlas_SALTIRE3_sc@meta.data
SALTIRE3_donor_colours_included <- SALTIRE3_donor_colours[sort(unique(sapply(strsplit(df_atlas_metadata$donor_id, "_"), `[`, 2)))]

df_atlas_metadata$donor <- sapply(strsplit(df_atlas_metadata$donor_id, "_"), `[`, 2)
df_atlas_metadata$donor <- factor(df_atlas_metadata$donor, levels = sort(unique(df_atlas_metadata$donor)))

plt1 <- ggplot(df_atlas_metadata, aes(x = donor)) +
  geom_bar(fill = SALTIRE3_donor_colours_included) +
  theme_minimal(base_size = 14) +
  scale_y_log10(expand = expansion(mult = c(0, 0.05)), breaks = c(10, 100, 1000, 10000, 100000)) + labs(
    y = "Number of Cells"
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

# Plot: Number of UMIs per cell by sample
plt2 <- ggplot(df_atlas_metadata[df_atlas_metadata$nCount_RNA < 40000, ], aes(x = donor, y = nCount_RNA, fill = donor)) +
  geom_violin(linewidth = 0.2, width = 1) +
  scale_fill_manual(values = SALTIRE3_donor_colours_included) +
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

plt3 <- ggplot(df_atlas_metadata, aes(x = donor, y = nFeature_RNA, fill = donor)) +
  geom_violin(linewidth = 0.2, width = 1) +
  scale_fill_manual(values = SALTIRE3_donor_colours_included) +
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
ggsave(plot = plot_grid(plotlist = aligned, ncol = 1, align = "v"), filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_8b.svg"), width = 3, height = 6, units = "cm")
write.csv(df_atlas_metadata[, c("donor", "nCount_RNA", "nFeature_RNA")], file = paste0(source_data_output_path, "Supplemental_Figure_8b.csv"), row.names = TRUE)
