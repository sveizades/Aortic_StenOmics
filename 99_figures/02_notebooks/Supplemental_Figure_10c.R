# Supplemental Figure 10c: Cell counts, UMIs and genes per cell by sample after filtering (Xenium atlas QC violins)

library(Seurat)
library(ggplot2)
library(cowplot)

source("../../utils.R")

atlas_SALTIRE3_Xenium <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")
df_atlas_metadata <- atlas_SALTIRE3_Xenium@meta.data

SALTIRE3_donor_colours_included <- SALTIRE3_donor_colours[sort(unique(sapply(strsplit(df_atlas_metadata$donor_id, "_"), `[`, 2)))]

df_atlas_metadata$donor <- sapply(strsplit(df_atlas_metadata$donor_id, "_"), `[`, 2)
df_atlas_metadata$donor <- factor(df_atlas_metadata$donor, levels = sort(unique(df_atlas_metadata$donor)))

df_atlas_metadata$sample_id <- gsub("^CureAS_", "", df_atlas_metadata$sample_id) |> gsub("_", " ", x = _)
df_atlas_metadata$sample_id <- factor(df_atlas_metadata$sample_id, levels = sort(unique(df_atlas_metadata$sample_id)))

plt1 <- ggplot(df_atlas_metadata, aes(x = sample_id, fill = donor)) +
  geom_bar() +
  scale_fill_manual(values = SALTIRE3_donor_colours_included) +
  theme_minimal(base_size = 5) +
  scale_y_log10(expand = expansion(mult = c(0, 0.05)), breaks = c(10, 100, 1000, 10000, 100000)) +
  labs(
    y = "Number of Cells"
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", size = 5),
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text = element_text(colour = "black", size = 5),
        plot.margin = unit(c(0, 0, 0, 1), "mm"),
        legend.position = "none")

# Plot: Number of UMIs per cell by sample
plt2 <- ggplot(df_atlas_metadata,
               aes(x = sample_id, y = nCount_Xenium, fill = donor)) +
  geom_violin(linewidth = 0.2, width = 1, adjust = 3, scale = "width", trim = FALSE) +
  scale_fill_manual(values = SALTIRE3_donor_colours_included) +
  labs(x = "Sample", y = "Number of UMI") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", size = 5),
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text = element_text(colour = "black", size = 5),
        plot.margin = unit(c(0, 0, 0, 1), "mm"),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = 0))

plt3 <- ggplot(df_atlas_metadata,
               aes(x = sample_id, y = nFeature_Xenium, fill = donor)) +
  geom_violin(linewidth = 0.2, width = 1, adjust = 3, scale = "width", trim = FALSE) +
  scale_fill_manual(values = SALTIRE3_donor_colours_included) +
  labs(x = "Sample", y = "Number of genes") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", size = 5),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 5),
        plot.margin = unit(c(0, 0, 0, 1), "mm"),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = 0))

ggsave(plot = plot_grid(plt1, plt2, plt3, ncol = 1, align = "v", axis = "lr", rel_heights = c(1.1, 1.1, 2)), filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10c.svg"), width = 3, height = 6, units = "cm")
write.csv(df_atlas_metadata[, c("donor", "nCount_Xenium", "nFeature_Xenium")], file = paste0(source_data_output_path, "Supplemental_Figure_10c.csv"), row.names = TRUE)
