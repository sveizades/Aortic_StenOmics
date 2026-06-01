# Figure 1f: Heatmap of differentially expressed genes across disease severity groups (UBC bulk RNA-seq).
library(ggplot2)
library(pheatmap)
library(gprofiler2)
library(DESeq2)
library(openxlsx)

source("../../utils.R")

# Load data
bulk_DE_AS_v_noAS <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_AS_v_noAS.csv", row.names = 1)
bulk_DE_mm_v_ctrl <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_mm_v_ctrl.csv", row.names = 1)
bulk_DE_sev_v_ctrl <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_sev_v_ctrl.csv", row.names = 1)
bulk_DE_sev_v_mm <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_sev_v_mm.csv", row.names = 1)
bulk_vsd <- readRDS("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/bulk_vsd.rds")

anyDEG <- unique(c(bulk_DE_AS_v_noAS[abs(bulk_DE_AS_v_noAS$log2FoldChange) > 1 & bulk_DE_AS_v_noAS$padj < 0.05 & !is.na(bulk_DE_AS_v_noAS$padj), "ensembl_gene_id"],
                   bulk_DE_mm_v_ctrl[abs(bulk_DE_mm_v_ctrl$log2FoldChange) > 1 & bulk_DE_mm_v_ctrl$padj < 0.05 & !is.na(bulk_DE_mm_v_ctrl$padj), "ensembl_gene_id"],
                   bulk_DE_sev_v_ctrl[abs(bulk_DE_sev_v_ctrl$log2FoldChange) > 1 & bulk_DE_sev_v_ctrl$padj < 0.05 & !is.na(bulk_DE_sev_v_ctrl$padj), "ensembl_gene_id"],
                   bulk_DE_sev_v_mm[abs(bulk_DE_sev_v_mm$log2FoldChange) > 1 & bulk_DE_sev_v_mm$padj < 0.05 & !is.na(bulk_DE_sev_v_mm$padj), "ensembl_gene_id"]))
vsd_scaled <- t(scale(t(assay(bulk_vsd)[anyDEG, ])))
vsd_scaled[vsd_scaled < -2.5] <- -2.5
vsd_scaled[vsd_scaled > 2.5] <- 2.5
vsd_scaled <- vsd_scaled[, colnames(bulk_vsd)[order(as.character(bulk_vsd$clinical_classification))]]

conv <- gconvert(
  rownames(vsd_scaled),
  organism = "hsapiens",
  target = "ENTREZGENE",
  filter_na = FALSE
)
conv <- conv[match(rownames(vsd_scaled), conv$input), ]
gene.symbols <- conv$target
gene.symbols[is.na(gene.symbols)] <- rownames(vsd_scaled)[is.na(gene.symbols)]
rownames(vsd_scaled) <- gene.symbols

annotation_col <- data.frame(row.names = colnames(bulk_vsd), clinical_classification = as.character(bulk_vsd$clinical_classification))
plt <- pheatmap(vsd_scaled,
                annotation_col = annotation_col, color = Seurat::PurpleAndYellow(k = 10),
                show_rownames = TRUE,
                show_colnames = FALSE, cluster_cols = FALSE,
                clustering_method = "ward.D2", annotation_colors = list(clinical_classification = UBC_disease_colours), treeheight_row = 0, legend = FALSE, annotation_legend = FALSE, fontsize = 5)
plt$gtable$grobs[[4]]$label <- ""
plt_flags <- add.flag(plt, kept.labels = c("FAP", "FN1", "COL1A1", "COL1A2", "COL3A1", "THBS2", "THY1", "CCL11", "FCGR3A", "COL10A1",
                                           "APOE", "LPL", "CNTFR", "NRXN1", "NFKBIA", "PDK6", "IL6R", "CXCL1", "CXCL5", "COMP", "CDKN2A", "CTHRC1", "ASPN", "SPARC", "PDK4"), repel.degree = 0)
ggsave(plot = plt_flags, filename = paste0(main_figure_output_path, "Figure_1f.svg"), width = 6, height = 6, units = "cm")

source_heatmap_matrix <- as.data.frame(vsd_scaled)
source_heatmap_matrix <- tibble::rownames_to_column(source_heatmap_matrix, "gene_id_or_symbol")
source_sample_annotation <- annotation_col |>
  tibble::rownames_to_column("sample_id")
wb <- createWorkbook()
addWorksheet(wb, "heatmap_matrix")
writeData(wb, "heatmap_matrix", source_heatmap_matrix)
addWorksheet(wb, "sample_annotation")
writeData(wb, "sample_annotation", source_sample_annotation)
saveWorkbook(
  wb,
  file = file.path(source_data_output_path, "Figure_1f.xlsx"),
  overwrite = TRUE
)

# Save scale bar
svg(paste0(main_figure_output_path, "Figure_1f_color_scale.svg"), width = 8, height = 2)
par(mar = c(0, 0, 0, 0))
image(seq(0, 1000), c(0, 1), matrix(1:1000, nrow = 1000, ncol = 1), col = Seurat::PurpleAndYellow(1000), axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
dev.off()
