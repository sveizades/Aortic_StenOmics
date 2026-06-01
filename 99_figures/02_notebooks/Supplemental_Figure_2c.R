# Supplemental Figure 2c: GO pathway enrichment bubble plot across disease conditions (UBC bulk RNA-seq).
library(ggplot2)
library(dplyr)
library(scico)
library(tidyr)
library(stringr)
library(cowplot)

source("../../utils.R")

bulk_pathways_ctrl <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/control_pathways.csv", row.names = 1)
bulk_pathways_mm <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/mildmoderate_pathways.csv", row.names = 1)
bulk_pathways_sev <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/severe_pathways.csv", row.names = 1)

# Pathways to plot
selected_pathways <- c("GO:0060348", "GO:0001649", "GO:0045933", "GO:0032963", "GO:0030199", "GO:0030198", "GO:0043062", "GO:0032963", "GO:0001501", "GO:0001503", "GO:0032964", "GO:0051216", "GO:0002062", "GO:0060350", "GO:0002063", "GO:0030595", "GO:0042730", "GO:0050727")

# Add condition labels
bulk_pathways_ctrl$condition <- "control"
bulk_pathways_mm$condition <- "mildmoderate"
bulk_pathways_sev$condition <- "severe"

# Combine all
pathway_all <- bind_rows(
  bulk_pathways_ctrl,
  bulk_pathways_mm,
  bulk_pathways_sev
)

# Keep only selected IDs and their Descriptions
selected_df <- pathway_all %>%
  filter(ID %in% selected_pathways) %>%
  distinct(ID, Description)

# Join selected Descriptions back in full data (avoids filtering too early)
pathway_all <- pathway_all %>%
  inner_join(selected_df, by = c("ID", "Description")) %>%
  complete(Description, condition, fill = list(
    pvalue = NA,
    FoldEnrichment = NA,
    p.adjust = NA
  )) %>%
  mutate(
    p_for_plot = ifelse(is.na(p.adjust), 1, ifelse(p.adjust < 1e-300, 1e-300, p.adjust)),
    logp = -log10(p_for_plot),
    is_significant = p_for_plot < 0.05,
    size_val = ifelse(is_significant, logp, 0.01),
    color_val = ifelse(is_significant, FoldEnrichment, NA)
  )

pathway_all$condition <- plyr::mapvalues(pathway_all$condition, from = c("control", "mildmoderate", "severe"), to = c("Control", "Mild/moderate", "Severe"))
pathway_all$Description <- str_to_sentence(pathway_all$Description)

plt <- ggplot(pathway_all, aes(x = condition, y = Description)) +
  geom_point(aes(size = size_val, fill = color_val), stroke = 0, shape = 21) +
  scale_fill_gradientn(colours = scico(10, palette = "acton", direction = -1)[1:5], na.value = "black") +
  scale_size(range = c(0.01, 3), breaks = c(1, 3, 5), name = expression(-log[10](p))) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.05),
        axis.text.y = element_text(colour = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 0.7, linetype = "solid", colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.position = "none") +
  coord_equal()

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_2c.svg"), width = 6, height = 6, units = "cm")

write.csv(
  pathway_all,
  file = paste0(source_data_output_path, "Supplemental_Figure_2c.csv"),
  row.names = FALSE
)

# Save scale bar
df_legend <- data.frame(x = 1, y = 1, fill = 1:5)
# Create a plot with just the color scale
legend_plot <- ggplot(df_legend, aes(x = x, y = y, fill = fill)) +
  geom_tile() +
  scale_fill_gradientn(colours = scico(10, palette = "acton", direction = -1)[1:5]) +
  theme_void() +
  theme(legend.position = "right")
# Save only the legend as a plot
ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_2c_color_scale.svg"), plot = legend_plot, width = 3, height = 5)

plt_size_legend <- plt +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 5, colour = "black"),
    legend.title = element_text(size = 5, colour = "black"),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm")
  ) +
  guides(
    fill = "none",
    size = guide_legend(
      override.aes = list(fill = "white", colour = "black", shape = 21, stroke = 0.2)
    )
  )

size_legend <- cowplot::get_legend(plt_size_legend)

ggsave(
  filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_2c_dot_scale.svg"),
  plot = size_legend,
  width = 2.2,
  height = 3,
  units = "cm",
  bg = "transparent"
)
