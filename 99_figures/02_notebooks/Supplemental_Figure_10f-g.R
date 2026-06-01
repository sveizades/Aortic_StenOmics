# Supplemental Figure 10f-g: Cross-modality VIC marker correlation heatmaps (f) and FAP+ osteogenic DEG bubble plot across scRNA-seq, snRNA-seq, and Xenium (g)

library(Seurat)
library(ggplot2)
library(scico)
library(dplyr)
library(tidyr)
library(openxlsx)

source("../../utils.R")

wb <- createWorkbook()

atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")
atlas_lim_nodular_tissue_vics <- subset(atlas_lim_nodular_tissue, annotations_level1 == "Valvular interstitial cell")

atlas_lim_nodular_tissue_vics$annotations_level2 <- plyr::mapvalues(atlas_lim_nodular_tissue_vics$annotations_level2, from = names(VIC_celltype_names), to = as.character(VIC_celltype_names))
atlas_lim_nodular_tissue_vics$annotations_level2 <- factor(atlas_lim_nodular_tissue_vics$annotations_level2, levels = sort(as.character(VIC_celltype_names)))

sn_vics <- readRDS("../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/vics.rds")

sn_vics$annotations_level2_readable <- factor(sn_vics$annotations_level2_readable, levels = sort(as.character(VIC_solo_celltype_levels)))

sc_vics <- readRDS("../../03_SALTIRE3_scRNA_seq/03_output/03_S3_scRNAseq_vic_annotation/vics.rds")

sc_vics$annotations_level2_readable <- as.character(sc_vics$annotations_level2_readable)
sc_vics$annotations_level2_readable <- plyr::mapvalues(sc_vics$annotations_level2_readable, from = names(VIC_celltype_names), to = as.character(VIC_celltype_names))
sc_vics$annotations_level2_readable <- factor(sc_vics$annotations_level2_readable, levels = sort(as.character(VIC_celltype_names)))

Idents(atlas_lim_nodular_tissue_vics) <- "annotations_level2"
Idents(sn_vics) <- "annotations_level2_readable"
Idents(sc_vics) <- "annotations_level2_readable"
deg_X <- FindAllMarkers(atlas_lim_nodular_tissue_vics, densify = T, max.cells.per.ident = 2000)
deg_sc <- FindAllMarkers(sc_vics, densify = T, max.cells.per.ident = 2000)
deg_sn <- FindAllMarkers(sn_vics, densify = T, max.cells.per.ident = 2000)

marker_genes <- c(
  "COL1A2",
  "COL3A1",
  "COL1A1",
  "LUM",
  "FN1",
  "CDH11",
  "RUNX2",
  "ANPEP",
  "ASPN",
  "FAP",
  "POSTN",
  "CRTAC1",
  "COMP"
)

osteogenic_cluster <- "FAP+ osteogenic"

deg_marker_df <- bind_rows(
  deg_sc %>%
    filter(cluster == osteogenic_cluster) %>%
    mutate(dataset = "scRNA-seq"),
  deg_sn %>%
    filter(cluster == osteogenic_cluster) %>%
    mutate(dataset = "snRNA-seq"),
  deg_X %>%
    filter(cluster == osteogenic_cluster) %>%
    mutate(dataset = "Xenium")
) %>%
  filter(gene %in% marker_genes) %>%
  mutate(
    dataset = factor(dataset, levels = c("scRNA-seq", "snRNA-seq", "Xenium")),
    gene = factor(gene, levels = rev(marker_genes)),
    neg_log10_p = -log10(p_val_adj),
    avg_log2FC_plot = pmax(pmin(avg_log2FC, 2), -2)
  )

plt <- ggplot(
  deg_marker_df,
  aes(
    x = dataset,
    y = gene,
    size = neg_log10_p,
    fill = avg_log2FC_plot
  )
) +
  geom_point(shape = 21, colour = "black", stroke = 0.3) +
  scale_size_continuous(
    name = expression(-log[10](p)),
    range = c(1.5, 5),
    breaks = c(5, 10, 20)
  ) +
  scale_fill_gradient2(
    name = "logFC",
    low = "#053061",
    mid = "white",
    high = "#7F2704",
    midpoint = 0,
    limits = c(-2, 2),
    breaks = c(-2, 0, 2)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )

ggsave(
  plot = plt,
  filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10g.svg"),
  width = 26.106,
  height = 59.972,
  units = "mm"
)
write.csv(deg_marker_df, file = paste0(source_data_output_path, "Supplemental_Figure_10g.csv"), row.names = TRUE)

deg_sc2 <- deg_sc %>%
  filter(p_val_adj < 0.05) %>%
  select(cluster, gene, avg_log2FC)

deg_sn2 <- deg_sn %>%
  filter(p_val_adj < 0.05) %>%
  select(cluster, gene, avg_log2FC)

deg_X2 <- deg_X %>%
  filter(p_val_adj < 0.05) %>%
  select(cluster, gene, avg_log2FC)

mat_sc <- deg_sc2 %>%
  pivot_wider(names_from = cluster,
              values_from = avg_log2FC)

mat_sn <- deg_sn2 %>%
  pivot_wider(names_from = cluster,
              values_from = avg_log2FC)

mat_X <- deg_X2 %>%
  pivot_wider(names_from = cluster,
              values_from = avg_log2FC)

mat_sc <- as.data.frame(mat_sc)
rownames(mat_sc) <- mat_sc$gene
mat_sc$gene <- NULL

mat_sn <- as.data.frame(mat_sn)
rownames(mat_sn) <- mat_sn$gene
mat_sn$gene <- NULL

mat_X <- as.data.frame(mat_X)
rownames(mat_X) <- mat_X$gene
mat_X$gene <- NULL

common_genes <- intersect(rownames(mat_sc), rownames(mat_sn))
mat_sc_sub <- mat_sc[common_genes, , drop = FALSE]
mat_sn_sub <- mat_sn[common_genes, , drop = FALSE]
mat_sc_sub[is.na(mat_sc_sub)] <- 0
mat_sn_sub[is.na(mat_sn_sub)] <- 0
cor_mat <- matrix(NA,
                  nrow = ncol(mat_sn_sub),
                  ncol = ncol(mat_sc_sub))
rownames(cor_mat) <- colnames(mat_sn_sub)
colnames(cor_mat) <- colnames(mat_sc_sub)
for (i in 1:ncol(mat_sn_sub)) {
  for (j in 1:ncol(mat_sc_sub)) {
    cor_mat[i, j] <- cor(mat_sn_sub[, i],
                         mat_sc_sub[, j],
                         method = "pearson")
  }
}
sn_levels <- sort(rownames(cor_mat))
sc_levels <- sort(colnames(cor_mat))
plot_df <- as.data.frame(as.table(cor_mat)) |>
  rename(sn = Var1, sc = Var2, r = Freq) |>
  mutate(
    sn = factor(sn, levels = sn_levels),
    sc = factor(sc, levels = sc_levels)
  )
plot_df$r[plot_df$r < -0.5] <- -0.5
plot_df$r[plot_df$r > 0.5] <- 0.5
plt <- ggplot(plot_df, aes(sc, sn, fill = r)) +
  geom_tile(color = "white", linewidth = 0.2) +
  coord_fixed() +
  scico::scale_fill_scico(palette = "vik", midpoint = 0) +
  labs(
    x = "scRNA-seq", y = "snRNA-seq"
  ) +
  theme_classic(base_size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.text = element_text(size = 5, colour = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 6, colour = "black")
  )
ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10f-scsn.svg"), width = 4.5, height = 4.5, units = "cm")
addWorksheet(wb, "scRNA-seq vs snRNA-seq")
writeData(wb, "scRNA-seq vs snRNA-seq", plot_df)
addWorksheet(wb, "scRNA-seq vs Xenium")

common_genes <- intersect(rownames(mat_sc), rownames(mat_X))
mat_sc_sub <- mat_sc[common_genes, , drop = FALSE]
mat_X_sub <- mat_X[common_genes, , drop = FALSE]
mat_sc_sub[is.na(mat_sc_sub)] <- 0
mat_X_sub[is.na(mat_X_sub)] <- 0
cor_mat <- matrix(NA,
                  nrow = ncol(mat_X_sub),
                  ncol = ncol(mat_sc_sub))
rownames(cor_mat) <- colnames(mat_X_sub)
colnames(cor_mat) <- colnames(mat_sc_sub)
for (i in 1:ncol(mat_X_sub)) {
  for (j in 1:ncol(mat_sc_sub)) {
    cor_mat[i, j] <- cor(mat_X_sub[, i],
                         mat_sc_sub[, j],
                         method = "pearson")
  }
}
X_levels <- sort(rownames(cor_mat))
sc_levels <- sort(colnames(cor_mat))
plot_df <- as.data.frame(as.table(cor_mat)) |>
  rename(X = Var1, sc = Var2, r = Freq) |>
  mutate(
    X = factor(X, levels = X_levels),
    sc = factor(sc, levels = sc_levels)
  )
plot_df$r[plot_df$r > 0.5] <- 0.5
plt <- ggplot(plot_df, aes(sc, X, fill = r)) +
  geom_tile(color = "white", linewidth = 0.2) +
  coord_fixed() +
  scico::scale_fill_scico(palette = "vik", midpoint = 0) +
  labs(
    x = "scRNA-seq", y = "Xenium"
  ) +
  theme_classic(base_size = 5) +
  theme(
    axis.text = element_text(size = 5, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 6, colour = "black")
  )
ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10f-scXenium.svg"), width = 4.5, height = 4.5, units = "cm")
writeData(wb, "scRNA-seq vs Xenium", plot_df)
addWorksheet(wb, "snRNA-seq vs Xenium")

common_genes <- intersect(rownames(mat_sn), rownames(mat_X))
mat_sn_sub <- mat_sn[common_genes, , drop = FALSE]
mat_X_sub <- mat_X[common_genes, , drop = FALSE]
mat_sn_sub[is.na(mat_sn_sub)] <- 0
mat_X_sub[is.na(mat_X_sub)] <- 0
cor_mat <- matrix(NA,
                  nrow = ncol(mat_X_sub),
                  ncol = ncol(mat_sn_sub))
rownames(cor_mat) <- colnames(mat_X_sub)
colnames(cor_mat) <- colnames(mat_sn_sub)
for (i in 1:ncol(mat_X_sub)) {
  for (j in 1:ncol(mat_sn_sub)) {
    cor_mat[i, j] <- cor(mat_X_sub[, i],
                         mat_sn_sub[, j],
                         method = "pearson")
  }
}
X_levels <- sort(rownames(cor_mat))
sn_levels <- sort(colnames(cor_mat))
plot_df <- as.data.frame(as.table(cor_mat)) |>
  rename(X = Var1, sn = Var2, r = Freq) |>
  mutate(
    X = factor(X, levels = X_levels),
    sn = factor(sn, levels = sn_levels)
  )
plot_df$r[plot_df$r > 0.5] <- 0.5
plt <- ggplot(plot_df, aes(X, sn, fill = r)) +
  geom_tile(color = "white", linewidth = 0.2) +
  coord_fixed() +
  scico::scale_fill_scico(palette = "vik", midpoint = 0) +
  labs(
    x = "snRNA-seq", y = "Xenium"
  ) +
  theme_classic(base_size = 5) +
  theme(
    axis.text = element_text(size = 5, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 6, colour = "black")
  )
ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10f-snXenium.svg"), width = 4.5, height = 4.5, units = "cm")
writeData(wb, "snRNA-seq vs Xenium", plot_df)
saveWorkbook(wb, paste0(source_data_output_path, "Supplemental_Figure_10f.xlsx"), overwrite = TRUE)
