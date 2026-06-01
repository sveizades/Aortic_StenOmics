# Supplemental Figure 10d: Total transcripts per gene for Xenium genes vs negative control probes and codewords

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

source("../../utils.R")

atlas_SALTIRE3_Xenium <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")

get_gene_totals <- function(assay_obj, category) {
  counts <- GetAssayData(assay_obj, layer = "counts")
  data.frame(
    total_transcripts = rowSums(counts) + 1,
    category = category
  )
}

df_plot <- bind_rows(
  get_gene_totals(atlas_SALTIRE3_Xenium@assays$Xenium,          "Xenium genes"),
  get_gene_totals(atlas_SALTIRE3_Xenium@assays$ControlCodeword,  "Neg Control\nCodeword"),
  get_gene_totals(atlas_SALTIRE3_Xenium@assays$ControlProbe,     "Neg Control\nProbe"),
  get_gene_totals(atlas_SALTIRE3_Xenium@assays$BlankCodeword,    "Blank\nCodeword")
)

# Order categories
df_plot$category <- factor(df_plot$category,
                           levels = c("Xenium genes", "Neg Control\nCodeword", "Neg Control\nProbe", "Blank\nCodeword"))

category_colours <- c(
  "Xenium genes"          = "#4C9BE8",
  "Neg Control\nCodeword" = "#E8754C",
  "Neg Control\nProbe"    = "#7DBE6A",
  "Blank\nCodeword"       = "#B87DDB"
)

plt <- ggplot(df_plot, aes(x = category, y = total_transcripts, fill = category)) +
  geom_violin(trim = FALSE, alpha = 0.7, colour = NA) +
  geom_jitter(width = 0.1, size = 0.2, shape = 16) +
  scale_fill_manual(values = category_colours) +
  scale_y_log10(labels = scales::comma) +
  labs(y = "Total transcripts per gene", x = NULL) +
  theme_classic(base_size = 5) +
  theme(
    axis.title.y = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5, colour = "black"),
    axis.text.y  = element_text(size = 5),
    axis.line  = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(0.5, "mm"),
    plot.margin = margin(1, 1, 1, 1, "mm"),
    legend.position = "none"
  )

ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_10d.svg"),
       plot = plt,
       width = 40, height = 40, units = "mm")
write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_10d.csv"), row.names = TRUE)
