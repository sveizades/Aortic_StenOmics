# Figure 2h: Pseudotime heatmap of gene expression dynamics along VIC lineages (UBC snRNA-seq).
library(pheatmap)
library(viridis)
library(grid)
library(Seurat)
library(stats)
library(slingshot)
library(plyr)
library(cowplot)
library(dplyr)

source("../../utils.R")

# Load data
vics <- readRDS("../../02_UBC_snRNA_seq/03_output/12_UBC_snRNAseq_vic_velocity/vics_seurat_slingshot.rds")
vic_sce <- readRDS("../../02_UBC_snRNA_seq/03_output/12_UBC_snRNAseq_vic_velocity/vics_sce_slingshot.rds")
vics$slingPseudotime_1_scaled <- vic_sce$slingPseudotime_1_scaled
vics$slingPseudotime_2_scaled <- vic_sce$slingPseudotime_2_scaled

# Let's say you have this from slingshot:
cellWeights <- slingCurveWeights(vic_sce)

# Assign each cell to its most likely trajectory
vics$lineage <- apply(cellWeights, 1, function(x) {
  if (is.na(x[1])) return(NA)
  paste0("Lineage", which.max(x))
})

order <- c(Cells(vics)[vics$lineage == "Lineage1"][order(vics$slingPseudotime_1_scaled[vics$lineage == "Lineage1"], decreasing = TRUE)], Cells(vics)[vics$lineage == "Lineage2"][order(vics$slingPseudotime_2_scaled[vics$lineage == "Lineage2"], decreasing = FALSE)])

vics$orderN <- as.numeric(mapvalues(Cells(vics), from = order, to = 1:length(order)))

groups <- cut(vics$orderN, breaks = 500)
Idents(vics) <- groups

cols <- c("RYR3", "BRINP1", "WARS2-IT1", "NTRK3", "EFEMP1", "DCLK1", "PDE7A", "FGF7", "TARID", "ADH4", "SELENOP", "PXDNL", "BCAS1", "EYS", "LINC01697", "TCF21", "EYA2", "ADH1B", "CCDC85A", "GLDN", "GRM7", "KAZN", "C7", "CHRM2", "CNTN6", "NPY1R", "AGTR1", "PTPRQ", "CDH20", "LINC00968", "RFTN1", "OSR1", "LINC03051", "RGS7BP", "GABRE", "ISM1", "NAV3", "COL4A3", "ABCA9-AS1", "SVEP1", "LINC01554", "MAFF", "LRP1B", "PODN", "NID2", "DENND2D", "EMID1", "SPON1", "MID1", "ELL2", "EPAS1", "SLC16A12", "FAM110B", "ABCA9", "ANOS1", "MYH10", "FBLN1", "CFH", "NRG1", "SLC9A9",
          "SLC6A20", "VIT", "ITGB4", "GRIK2", "COL9A1", "CLDN10", "SLC9C2", "CXCL2", "ELOVL2", "RBFOX1", "IL1RAPL1", "ANKRD45", "SLC6A1", "NRCAM", "FAM135B", "TNFRSF21", "LINC03021.1", "LINC00092",
          "CADM2", "SLC7A14-AS1", "CDH2", "BRINP3", "STK32A", "FRMD5", "ST3GAL6", "PRICKLE1", "PDZD2", "TM4SF1", "LTC4S", "PDE10A", "TRIM9", "ITGA6", "COBL", "TENM2", "SAMD12", "SPOCK3", "MYRIP", "CNBD1", "FRAS1", "STXBP5L", "PIP5K1B", "POSTN", "PRAG1", "FN1", "COMP", "DGKI", "PDGFC", "SLC20A1", "IVNS1ABP", "BHLHE40", "CHD7", "ADGRG2", "PODXL", "RGCC", "FAP", "CD55", "ENPP1", "RUNX1", "NGEF", "SIPA1L2", "UGP2", "GRIP1", "ENOX1", "KIAA1217", "CRTAC1", "HTRA1", "PMEPA1", "CD109", "ATP10A", "COL13A1", "SGCD", "GALNT5", "ITGA10", "SULF1", "PRELP", "ABI3BP", "COL15A1", "SERPINE2", "CDON", "UNC5B", "CREB3L2", "RETREG1", "TMEM164", "PAPSS2", "CHST3", "DISC1", "VAT1L", "XKR6", "ACVR1", "PRKCA", "ALCAM",
          "TPM2", "CPXM2", "GALNT18", "PLS3", "CCDC3", "CSRP1", "ABCC9", "MEF2C", "HDAC9", "DGKG", "SYTL2", "MYOM1", "CNNM2", "SYNE2", "LIMS2", "INHBA-AS1", "FGF1", "TBC1D1", "MCAM", "ADGRL3", "MYH11", "NFASC", "A2M", "ITGA7", "SLC8A1", "ADGRE5", "DENND3", "TBL1X", "CTNNA3", "MERTK", "ASPN", "SGCD", "FAM13C", "PAG1", "PRUNE2", "RIN3", "MTSS1", "PCED1B", "INPP4B", "YWHAH", "ACTA2", "NR2F2", "FLNA", "MYH9", "LAMA5", "CACNA1C", "CRIM1", "HIP1", "RASSF3", "MYLK", "PELI2", "CLMN", "ACTA2", "CNN1", "MYH11", "CADM2", "KAZN", "C7", "CD55", "ENPP1", "COMP", "CRTAC1", "POSTN", "FAP", "RUNX1")
cols <- unique(cols)
groups <- cut(vics$orderN, breaks = 100)
Idents(vics) <- groups
group_cts <- paste(groups, vics$annotations_level2, sep = "&")
group_cts <- table(group_cts)
names_vector <- names(group_cts)
intervals <- sub("&.*", "", names_vector)
cell_types <- sub(".*&", "", names_vector)
count_df <- data.frame(Interval = intervals, CellType = cell_types, Count = unname(group_cts))
result <- count_df %>%
  dplyr::group_by(Interval, CellType) %>%
  dplyr::summarize(TotalCount = sum(Count.Freq)) %>%
  group_by(Interval) %>%
  dplyr::filter(TotalCount == max(TotalCount)) %>%
  ungroup()
result <- result %>%
  arrange(Interval) %>%
  mutate(
    PrevInterval = lag(Interval),
    NextInterval = lead(Interval)
  ) %>%
  dplyr::filter(!is.na(PrevInterval) | !is.na(NextInterval)) %>%
  group_by(Interval) %>%
  dplyr::summarise(
    MostCommonCellType = CellType[which.max(TotalCount)],
    TotalCount = max(TotalCount)
  ) %>%
  ungroup()
col_anno <- as.data.frame(result$MostCommonCellType)
colnames(col_anno) <- "Group"
# Gather gene expression values
avg_expr <- AverageExpression(vics, assays = "RNA", slot = "data")
avg_expr <- avg_expr[["RNA"]]
avg_expr <- avg_expr[rownames(avg_expr) %in% cols, ]
# Smoothing
for (row in 1:nrow(avg_expr)) {
  x <- 1:ncol(avg_expr)
  y <- avg_expr[row, ]
  smoothed <- loess(y ~ x, span = 0.2)$fitted
  avg_expr[row, ] <- smoothed
}
# Row annotations
row_anno <- rownames(avg_expr)
row_anno[!row_anno %in% c("ACTA2", "CNN1", "MYH11", "CADM2", "KAZN", "C7", "CD55", "ENPP1", "COMP", "CRTAC1", "POSTN", "FAP", "RUNX1")] <- ""
rownames(col_anno) <- colnames(avg_expr)

mat2 <- t(apply(avg_expr, 1, scale))
mat2[mat2 > 2.5] <- 2.5
mat2[mat2 < -2.5] <- -2.5
mat2 <- mat2[cols, ]
row_anno <- rownames(mat2)

mat2_lin1 <- mat2[, 1:60]
heat_lin1 <- pheatmap(mat2_lin1, cluster_cols = FALSE, cluster_rows = FALSE,
                      labels_col = rep("", ncol(avg_expr)),
                      show_rownames = FALSE,
                      color = c(viridis(1000, option = "C", begin = 0, end = 0.2),
                              viridis(300, option = "C", begin = 0.2, end = 0.8),
                              viridis(1000, option = "C", begin = 0.8, end = 1)), legend = FALSE)
mat2_lin2 <- mat2[, 61:100]
heat_lin2 <- pheatmap(mat2_lin2, cluster_cols = FALSE, cluster_rows = FALSE,
                      labels_col = rep("", ncol(avg_expr)),
                      labels_row = row_anno,
                      color = c(viridis(1000, option = "C", begin = 0, end = 0.2),
                              viridis(300, option = "C", begin = 0.2, end = 0.8),
                              viridis(1000, option = "C", begin = 0.8, end = 1)), legend = TRUE, fontsize = 5)

p7 <- plot_grid(heat_lin1[[4]], add.flag(heat_lin2,
                                         kept.labels = c("ACTA2", "CNN1", "MYH11", "CADM2", "KAZN", "C7", "CD55", "ENPP1", "COMP", "CRTAC1", "POSTN", "FAP", "RUNX1"),
                                         repel.degree = 0), ncol = 2)
heat_lin2 <- pheatmap(mat2_lin2, cluster_cols = FALSE, cluster_rows = FALSE,
                      labels_col = rep("", ncol(avg_expr)),
                      show_rownames = FALSE,
                      color = c(viridis(1000, option = "C", begin = 0, end = 0.2),
                              viridis(300, option = "C", begin = 0.2, end = 0.8),
                              viridis(1000, option = "C", begin = 0.8, end = 1)), legend = FALSE)
p8 <- plot_grid(heat_lin1[[4]], heat_lin2[[4]], ncol = 2)
save_plot(filename = paste0(main_figure_output_path, "Figure_2h-legends.svg"), plot = p7, base_height = 1.3)
save_plot(filename = paste0(main_figure_output_path, "Figure_2h.svg"), plot = p8, base_height = 1.3)
write.csv(mat2, file = paste0(source_data_output_path, "Figure_2h.csv"), row.names = TRUE)
