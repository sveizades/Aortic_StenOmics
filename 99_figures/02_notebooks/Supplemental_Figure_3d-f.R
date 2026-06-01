# Supplemental Figure 3d-f: Pre-filtering QC histograms for nFeature_RNA, nCount_RNA, and percent mitochondrial transcripts (UBC snRNA-seq).
library(Seurat)
library(ggplot2)

source("../../utils.R")

alas_pre_filtering <- readRDS("../../02_UBC_snRNA_seq/03_output/01_UBC_snRNAseq_preprocessing/UBCnucatlas_unprocessed.rds")

df <- alas_pre_filtering@meta.data

# Histogram with threshold lines
plt <- ggplot(df[df$nFeature_RNA < 10000, ], aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill = "#001861") +
  geom_vline(xintercept = 384.5321, linetype = "dashed", color = "red", size = 0.2) +
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

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_3d.svg"), width = 3.5, height = 3, units = "cm")

plt <- ggplot(df[df$nCount_RNA < 30000, ], aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 300, fill = "#001861") +
  geom_vline(xintercept = 360.1015, linetype = "dashed", color = "red", size = 0.2) +
  labs(x = "Number of UMIs", y = "Cell count") +
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
  scale_y_continuous(expand = expansion(mult = 0)) +
  scale_x_continuous(breaks = c(0, 15000, 30000))

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_3e.svg"), width = 3, height = 3, units = "cm")

plt <- ggplot(df[df$percent_mito < 15, ], aes(x = percent_mito)) +
  geom_histogram(binwidth = 0.25, fill = "#001861") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 0.2) +
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

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_3f.svg"), width = 3, height = 3, units = "cm")

write.csv(df[, c("nFeature_RNA", "nCount_RNA", "percent_mito")], file = paste0(source_data_output_path, "Supplemental_Figure_3d-f.csv"), row.names = TRUE)
