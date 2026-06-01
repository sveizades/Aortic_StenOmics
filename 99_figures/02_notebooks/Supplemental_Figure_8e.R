# Supplemental Figure 8e: Histogram of mitochondrial transcript percentage per cell with filtering threshold from unprocessed SALTIRE3 scRNA-seq data
library(Seurat)
library(ggplot2)

source("../../utils.R")

# Load data
atlas_SALTIRE3_sc_pre_filtering <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/01_S3_scRNAseq_preprocessing/SALTIRE3cellatlas_unprocessed.rds")
df <- atlas_SALTIRE3_sc_pre_filtering@meta.data

# Histogram with threshold lines
plt <- ggplot(df, aes(x = percent_mito)) +
  geom_histogram(binwidth = 0.25, fill = "#001861") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "red", linewidth = 0.2) +
  labs(x = "% mitochondrial transcripts", y = "Cell count") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 5),
        axis.title.x = element_text(colour = "black", size = 5),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = 0))
ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_8e.svg"), width = 3, height = 3, units = "cm")
