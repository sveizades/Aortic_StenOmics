# Supplemental Figure 10a: Histogram of genes per cell before filtering (Xenium pre-filtering QC)

library(Seurat)
library(ggplot2)

source("../../utils.R")

Xenium_pre_filtering <- readRDS("../../04_SALTIRE3_Xenium/03_output/01_SALTIRE3_xenium_preprocessing/Xenium_prefiltering.rds")
df <- Xenium_pre_filtering@meta.data

# Histogram with threshold lines
plt <- ggplot(df, aes(x = nFeature_Xenium)) +
  geom_histogram(binwidth = 5, fill = "#001861") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "red", linewidth = 0.2) +
  labs(x = "Number of genes", y = "Cell count") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 5),
        axis.title = element_text(colour = "black", size = 5),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = 0))

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10a.svg"), width = 3.5, height = 3, units = "cm")
write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_10a.csv"), row.names = TRUE)
