# Supplemental Figure 4c: CytoTRACE2 potency score boxplots across VIC subtypes (UBC snRNA-seq).
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(stringr)
library(ggrastr)

source("../../utils.R")

# Read CytoTRACE2 results
ct2 <- read_tsv("../../02_UBC_snRNA_seq/03_output/12_UBC_snRNAseq_vic_velocity/cytotrace2_results/cytotrace2_results.txt")
vics <- readRDS("../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/vics.rds")

# If the first column came in unnamed, rename it
colnames(ct2)[1] <- "cell_id_ct2"

# Make CytoTRACE2 barcodes match Seurat colnames
# turns CureAS_JU_TTATGCTAGTCCAGGA.1 -> CureAS_JU_TTATGCTAGTCCAGGA-1
ct2 <- ct2 %>%
  mutate(cell_id_seurat = str_replace(cell_id_ct2, "\\.([0-9]+)$", "-\\1"))

# Keep only cells present in the Seurat object
ct2 <- ct2 %>%
  filter(cell_id_seurat %in% colnames(vics))

# Add CytoTRACE2 metadata to Seurat object
meta_to_add <- ct2 %>%
  select(
    cell_id_seurat,
    CytoTRACE2_Score,
    CytoTRACE2_Potency,
    CytoTRACE2_Relative,
    preKNN_CytoTRACE2_Score,
    preKNN_CytoTRACE2_Potency
  ) %>%
  tibble::column_to_rownames("cell_id_seurat")

vics <- AddMetaData(vics, metadata = meta_to_add)

# Build plotting dataframe
plot_df <- vics@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(
    cell_id,
    phenotype = annotations_level2_readable,
    CytoTRACE2_Score,
    CytoTRACE2_Potency
  ) %>%
  filter(!is.na(phenotype), !is.na(CytoTRACE2_Score))

plot_df$phenotype <- factor(
  plot_df$phenotype,
  levels = c("Spongiosa", "IFN-stimulated", "Transitional", "Quiescent", "FAP+ osteogenic", "Neural crest-like", "Contractile")
)

# Optional: fill colours similar to your figure
vic_cols <- c(
  "Neural crest-like" = "#5DA0A3",
  "Quiescent" = "#5DA0A3",
  "FAP+ osteogenic" = "#67ABA2",
  "Transitional" = "#67ABA2",
  "Contractile" = "#5A4A99",
  "IFN-stimulated" = "#67ABA2",
  "Spongiosa" = "#73B6A1"
)

potency_labels <- data.frame(
  y = c(0.92, 0.75, 0.58, 0.41, 0.25, 0.08),
  label = c("Totipotent", "Pluripotent", "Multipotent",
            "Oligopotent", "Unipotent", "Differentiated")
)

p <- ggplot(plot_df, aes(x = phenotype, y = CytoTRACE2_Score, fill = phenotype)) +
  geom_boxplot(
    width = 0.8,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  ggrastr::rasterise(
    geom_jitter(
      width = 0.2,
      size = 0.01,
      alpha = 0.1,
      colour = "black"
    ),
    dpi = 300
  ) +
  scale_fill_manual(values = vic_cols, drop = FALSE) +
  labs(
    y = "Potency score"
  ) +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 5, colour = "black"),
    axis.text.x = element_text(size = 5, colour = "black", angle = 45, hjust = 1, vjust = 1.05),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 5, colour = "black"),
    plot.title = element_text(size = 5, colour = "black", hjust = 0.5),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(0, 10, 0, 0), "mm")
  ) +
  annotate(
    "text",
    x = length(levels(plot_df$phenotype)),
    y = c(0.92, 0.75, 0.58, 0.41, 0.25, 0.08),
    label = c("Totipotent", "Pluripotent", "Multipotent",
              "Oligopotent", "Unipotent", "Differentiated"),
    hjust = 0,
    size = 5 / .pt
  )

ggsave(
  paste0(supplemental_figure_output_path, "Supplemental_Figure_4c.svg"),
  p,
  width = 8,
  height = 6,
  units = "cm"
)

write.csv(plot_df, file = paste0(source_data_output_path, "Supplemental_Figure_4c.csv"), row.names = TRUE)
