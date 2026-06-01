# Supplemental Figure 3b: Donor similarity matrix heatmap for UBC snRNA-seq demultiplexing quality control.
library(pheatmap)
library(ggplot2)
library(scico)

source("../../utils.R")

SQ_similaritymatrix <- read.csv("../../02_UBC_snRNA_seq/03_output/02_UBC_snRNAseq_dehashing/CureAS_SQ_donor_similarity_matrix.csv", row.names = 1)
colnames(SQ_similaritymatrix) <- factor(sapply(strsplit(colnames(SQ_similaritymatrix), "_"), `[`, 2), levels = rev(c("INV", "UKG", "O2N", "PUM")))
rownames(SQ_similaritymatrix) <- factor(paste0("SQ ", rownames(SQ_similaritymatrix)), levels = c("SQ 2", "SQ 3", "SQ 1", "SQ 0"))

heatmap <- pheatmap(
  SQ_similaritymatrix,
  cluster_rows = F,
  cluster_cols = F,
  fontsize_number = 5,
  fontsize_row = 5,
  fontsize_col = 5,
  color = scico(100, palette = "vik"),
  border_color = NA,
  treeheight_row = 0,
  treeheight_col = 0,
  legend = FALSE, cutree_rows = 1, cutree_cols = 1,
  angle_col = 45
)

write.csv(SQ_similaritymatrix, file = paste0(source_data_output_path, "Supplemental_Figure_3b.csv"), row.names = TRUE)

ggsave(plot = heatmap, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_3b.svg"), width = 3, height = 3, units = "cm")
