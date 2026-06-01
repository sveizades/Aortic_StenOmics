# Supplemental Figure 8c: Histogram of number of genes per cell with filtering threshold from unprocessed SALTIRE3 scRNA-seq data
library(Seurat)
library(ggplot2)

source("../../utils.R")

# Load data
atlas_SALTIRE3_sc_pre_filtering <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/01_S3_scRNAseq_preprocessing/SALTIRE3cellatlas_unprocessed.rds")
df <- atlas_SALTIRE3_sc_pre_filtering@meta.data

# Histogram with threshold lines
plt <- ggplot(df[df$nFeature_RNA < 10000, ], aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill = "#001861") +
  geom_vline(xintercept = 350, linetype = "dashed", color = "red", linewidth = 0.2) +
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
  scale_y_continuous(expand = expansion(mult = 0)) +
  scale_x_continuous(breaks = c(0, 5000, 10000))
ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_8c.svg"), width = 3.5, height = 3, units = "cm")
write.csv(df[, c("nCount_RNA", "nFeature_RNA", "percent_mito")], file = paste0(source_data_output_path, "Supplemental_Figure_8c-e.csv"), row.names = TRUE)
