# Figure 4h: Xenium cell-state enrichment across spatial niches

library(dplyr)
library(tidyr)
library(scico)
library(ggplot2)

source("../../utils.R")

atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")
atlas_lim_nodular_tissue$niche <- factor(atlas_lim_nodular_tissue$niche, levels = c("Cap", "Spongiosa", "Peri-nodular", "Immune"))

meta <- atlas_lim_nodular_tissue@meta.data 

# Proportion of each cell type within each niche, per donor
df_niche_prop <- meta %>%
  dplyr::count(donor_id, niche, annotations_level2, name = "n") %>%
  dplyr::group_by(donor_id, niche) %>%
  dplyr::mutate(
    prop_within_niche = n / sum(n)
  ) %>%
  dplyr::ungroup()

# Enrichment-style summary:
# observed = mean proportion of a cell type in a niche across donors
# expected = mean proportion of that same cell type in all other niches
df_enrich <- df_niche_prop %>%
  dplyr::group_by(annotations_level2, niche) %>%
  dplyr::summarise(
    observed = mean(prop_within_niche, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    expected = mean(
      df_niche_prop$prop_within_niche[
        df_niche_prop$annotations_level2 == annotations_level2 &
          df_niche_prop$niche != niche
      ],
      na.rm = TRUE
    ),
    log2_fc = log2((observed + 1e-6) / (expected + 1e-6))
  ) %>%
  dplyr::ungroup()

# Read p-values
p_long <- read.csv("../../04_SALTIRE3_Xenium/01_data/scNiche_annotations_level2_pval.csv", row.names = 1) %>%
  tibble::rownames_to_column("annotations_level2") %>%
  tidyr::pivot_longer(
    cols = -annotations_level2,
    names_to = "niche",
    values_to = "p"
  )


p_long$niche <- plyr::mapvalues(p_long$niche, from = paste0("Niche", 0:3), c("Cap", "Spongiosa", "Peri-nodular", "Immune"))
p_long$niche <- factor(p_long$niche, levels = c("Cap", "Spongiosa", "Peri-nodular", "Immune"))

# Join FC and p-values into one long plotting dataframe
df_plot <- df_enrich %>%
  dplyr::left_join(p_long, by = c("annotations_level2", "niche")) %>%
  dplyr::mutate(
    niche = factor(niche, levels = c("Cap", "Spongiosa", "Peri-nodular", "Immune"))
  )
# Cap FC values for plotting if desired
df_plot <- df_plot %>%
  dplyr::mutate(
    log2_fc_plot = pmax(pmin(log2_fc, 1), -1)
  )


## Build matrix only for row ordering
mat <- df_plot %>%
  dplyr::select(annotations_level2, niche, log2_fc_plot) %>%
  tidyr::pivot_wider(
    names_from = niche,
    values_from = log2_fc_plot,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("annotations_level2") %>%
  as.matrix()

row_order <- c(
  "VIC-transitional",
  "VIC-contractile",
  "VIC-spongiosa",
  "VIC-neural-crest",
  "VIC-quiescent",
  "VIC-IFN",
  "VIC-FAP-osteogenic",
  "Macrophage",
  "Monocyte",
  "VEC-endocardial",
  "VEC-vascular",
  "Dendritic cell",
  "Neutrophil",
  "Mast cell",
  "NK",
  "T-CD4",
  "T-CD8",
  "T-gd",
  "B cell",
  "Plasma cell",
  "Proliferating"
)


df_plot <- df_plot %>%
  dplyr::mutate(
    annotations_level2 = factor(annotations_level2, levels = rev(row_order))
  )

# Significant tiles
df_sig <- df_plot %>%
  dplyr::filter(!is.na(p), p < 0.05)

plt <- ggplot(df_plot, aes(x = niche, y = annotations_level2, fill = log2_fc_plot)) +
  geom_tile(color = "white", linewidth = 0.1) +
  geom_tile(
    data = df_sig,
    aes(x = niche, y = annotations_level2),
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_scico(
    palette = "roma",
    direction = -1,
    limits = c(-1, 1)
  ) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      colour = "black",
      size = 5
    ),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(colour = "black", size = 5)
  ) +
  coord_equal()


ggsave(plot = plt, 
       filename = paste0(main_figure_output_path, "Figure_4h.svg"), 
       width = 3.5, 
       height = 5.5, 
       units = "cm")
write.csv(df_plot, file = paste0(source_data_output_path, "Figure_4h.csv"), row.names = TRUE)
