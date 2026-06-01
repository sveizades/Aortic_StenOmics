# Supplemental Figure 10l: Osteogenic module score by spatial niche in VICs (Xenium atlas)

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(scico)
library(rstatix)

source("../../utils.R")

atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")
Idents(atlas_lim_nodular_tissue) <- "niche"

osteogenic_genes <- c(
  # Master TFs
  "RUNX2",     # osteoblast master regulator
  "SOX9",      # chondrogenic master regulator

  # Bone matrix / osteoblast markers
  "BGLAP",     # osteocalcin
  "IBSP",      # bone sialoprotein
  "SPP1",      # osteopontin
  "ALPL",      # alkaline phosphatase (mineralisation)

  # Cartilage / osteochondral matrix
  "ACAN",      # aggrecan
  "CHAD",      # chondroadherin
  "COL10A1",   # hypertrophic chondrocyte (endochondral ossification)
  "COL2A1",    # type II collagen (cartilage)
  "COL11A1",   # cartilage-associated collagen
  "COMP",      # cartilage oligomeric matrix protein

  # Mineralisation / bone remodelling
  "ENPP1",     # pyrophosphate metabolism
  "BMP2",      # osteogenic morphogen
  "TNFSF11"    # RANKL
)

# Score all VICs by niche
niche_col <- "niche"
atlas_lim_nodular_tissue$niche <- factor(atlas_lim_nodular_tissue$niche, levels = rev(c("Cap", "Spongiosa", "Peri-nodular", "Immune")))

atlas_lim_nodular_tissue <- AddModuleScore(
  atlas_lim_nodular_tissue,
  features = list(osteogenic = osteogenic_genes),
  name = "osteogenic_score",
  ctrl = 20,
  nbin = 12
)

# subset to VICs only
vic_cells <- atlas_lim_nodular_tissue@meta.data %>%
  filter(startsWith(annotations_level2, "VIC"))

# Violin/box plot of osteogenic score per niche
niche_order <- c("Peri-nodular", "Cap", "Spongiosa", "Immune")

comparisons <- lapply(niche_order[-1], function(x) c("Peri-nodular", x))

stat_df <- vic_cells %>%
  wilcox_test(osteogenic_score1 ~ niche, ref.group = "Peri-nodular") %>%
  adjust_pvalue(method = "BH") %>%
  mutate(
    p_label = ifelse(p.adj < 0.001, "p<0.001", paste0("p=", round(p.adj, 3))),
    y.position = 2 * c(1.15, 1.3, 1.45)
  )

p_score <- ggplot(
  vic_cells,
  aes(x = niche, y = osteogenic_score1, fill = niche)
) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  stat_pvalue_manual(stat_df, label = "p_label", size = 1.5, bracket.size = 0.3) +
  scale_fill_manual(values = spatial_niche_cols) +
  labs(y = "Osteogenic module score") +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5, colour = "black"),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(0.1, "cm")
  )

ggsave(
  plot = p_score,
  filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_10l.svg"),
  width = 4, height = 4, units = "cm"
)
write.csv(vic_cells, file = paste0(source_data_output_path, "Supplemental_Figure_10l.csv"), row.names = TRUE)
