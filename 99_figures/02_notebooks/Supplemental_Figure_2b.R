# Supplemental Figure 2b: Bulk RNA-seq gene expression boxplots for selected inflammatory and metabolic genes (UBC cohort).
library(dplyr)
library(readxl)
library(DESeq2)
library(ggplot2)
library(grid)
library(Seurat)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(gprofiler2)

source("../../utils.R")

bulk_dds <- readRDS("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DESeq2_condition_obj.rds")
bulk_DE_mm_v_ctrl <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_mm_v_ctrl.csv", row.names = 1)
bulk_DE_sev_v_ctrl <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_sev_v_ctrl.csv", row.names = 1)
bulk_DE_sev_v_mm <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_sev_v_mm.csv", row.names = 1)

# Genes to plot
selected_genes <- c("APOE", "CNTFR", "NRXN1", "NFKBIA", "PDK4", "IL6R")
selected_genes_ensembl <- gconvert(selected_genes, organism = "hsapiens", target = "ENSG", filter_na = F)$target

# Normalise counts to log (counts/counts + 1)
normalised_counts <- Seurat::NormalizeData(assay(bulk_dds))[selected_genes_ensembl, ]
rownames(normalised_counts) <- selected_genes

# Convert to long format and join with metadata
vst_long <- normalised_counts %>%
  as.data.frame() %>%
  rownames_to_column("hgnc_symbol") %>%
  pivot_longer(-hgnc_symbol, names_to = "donor_id", values_to = "vst_expr") %>%
  left_join(colData(bulk_dds) %>% as.data.frame() %>% rownames_to_column("donor_id"), by = "donor_id")

# Extract relevant DE results for selected genes & significant p-values only (padj < 0.05)
df_ctrl_mm <- bulk_DE_mm_v_ctrl %>%
  filter(hgnc_symbol %in% selected_genes, padj < 0.05) %>%
  mutate(group1 = "control", group2 = "mildmoderate", comparison = "Ctrl vs Mild")

df_ctrl_sev <- bulk_DE_sev_v_ctrl %>%
  filter(hgnc_symbol %in% selected_genes, padj < 0.05) %>%
  mutate(group1 = "control", group2 = "severe", comparison = "Ctrl vs Severe")

df_mm_sev <- bulk_DE_sev_v_mm %>%
  filter(hgnc_symbol %in% selected_genes, padj < 0.05) %>%
  mutate(group1 = "mildmoderate", group2 = "severe", comparison = "Mild vs Severe")

# Combine all
pvals_df <- bind_rows(df_ctrl_mm, df_ctrl_sev, df_mm_sev) %>%
  dplyr::select(hgnc_symbol, group1, group2, p.value = padj, comparison)

offset_multipliers <- tibble(
  comparison = c("Ctrl vs Mild", "Ctrl vs Severe", "Mild vs Severe"),
  offset_factor = c(0.2, 0.05, 0.35)
)

max_expr_df <- vst_long %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::summarise(max_expr = max(vst_expr), .groups = "drop")

pvals_df <- pvals_df %>%
  left_join(max_expr_df, by = "hgnc_symbol") %>%
  group_by(hgnc_symbol) %>%
  arrange(hgnc_symbol, comparison) %>%
  dplyr::mutate(offset_rank = row_number()) %>%
  ungroup() %>%
  dplyr::mutate(y.position = max_expr * (1.05 + 0.15 * offset_rank))

pvals_df$pvalue_label <- ifelse(
  pvals_df$p.value < 0.001,
  "p<0.001",
  paste0("p=", round(pvals_df$p.value, 3))
)

plt <- ggplot(vst_long, aes(x = clinical_classification, y = vst_expr)) +
  stat_boxplot(aes(colour = clinical_classification), geom = "errorbar", linetype = 1, width = 0.5, alpha = 1) +
  geom_boxplot(outlier.shape = NA, aes(fill = clinical_classification, colour = clinical_classification)) +
  geom_jitter(aes(fill = clinical_classification, colour = clinical_classification), shape = 16,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              color = "black",
              alpha = 1,
              size = 0.7) +
  facet_wrap(~ hgnc_symbol, scales = "free_y", nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_cartesian(ylim = c(0, NA), clip = "on") +
  scale_color_manual(values = UBC_disease_colours) +
  scale_fill_manual(values = UBC_disease_colours_lighter) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 5),
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),         # removes grey box
        strip.text = element_text(color = "black", size = 7),
        plot.margin = unit(c(0, 0, 1, 0), "mm"),
        legend.position = "none") +
  stat_pvalue_manual(pvals_df,
                     label = "pvalue_label",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = "y.position",
                     tip.length = 0.01,
                     step.increase = 0.0, size = 1.75)

ggsave(plot = plt, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_2b.svg"), width = 2 * length(selected_genes), height = 4, units = "cm")

source_data <- vst_long %>%
  dplyr::select(
    hgnc_symbol,
    donor_id,
    clinical_classification,
    vst_expr
  )

write.csv(
  source_data,
  file = paste0(source_data_output_path, "Supplemental_Figure_2b.csv"),
  row.names = FALSE
)
