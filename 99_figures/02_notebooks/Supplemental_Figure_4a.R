# Supplemental Figure 4a: VIC subtype proportion boxplots across disease conditions (UBC snRNA-seq).
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

source("../../utils.R")

# Read in the saved statistical results and proportions data
stats_results <- read.csv(
  "../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/VIC_proportions_statistics.csv"
)

proportions_data <- read.csv(
  "../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/VIC_proportions_data.csv"
)

# VIC subtypes to exclude because they are shown separately in Figure 2d
excluded_celltypes <- c("FAP+ osteogenic", "Contractile")

# Prepare data for all remaining VIC subtypes
props_long <- proportions_data %>%
  mutate(
    clinical_classification = factor(
      clinical_classification,
      levels = c("control", "mildmoderate", "severe")
    ),
    celltype = factor(
      annotations_level2_readable,
      levels = VIC_solo_celltype_levels
    )
  ) %>%
  filter(!celltype %in% excluded_celltypes) %>%
  droplevels()

# Get p-values for the same VIC subtypes
stats_for_plot <- stats_results %>%
  filter(celltype %in% VIC_solo_celltype_levels) %>%
  mutate(
    celltype = factor(
      celltype,
      levels = VIC_solo_celltype_levels
    )
  ) %>%
  filter(!celltype %in% excluded_celltypes) %>%
  droplevels()

# Calculate max proportion across plotted celltypes for consistent p-value positioning
max_prop_overall <- max(props_long$proportion, na.rm = TRUE)

# Create p-value labels only for significant results
pval_labels <- stats_for_plot %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(
    pvalue_label = ifelse(
      adj.P.Val < 0.001,
      "p<0.001",
      paste0("p=", round(adj.P.Val, 3))
    ),
    y_position = max_prop_overall * 1.15
  )

# Create the plot
plt <- ggplot(
  props_long,
  aes(x = clinical_classification, y = proportion)
) +
  stat_boxplot(
    aes(colour = clinical_classification),
    geom = "errorbar",
    linetype = 1,
    width = 0.5,
    alpha = 1
  ) +
  geom_boxplot(
    outlier.shape = NA,
    aes(
      fill = clinical_classification,
      colour = clinical_classification
    ),
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(
      fill = clinical_classification,
      colour = clinical_classification
    ),
    shape = 16,
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.75
    ),
    colour = "black",
    alpha = 1,
    size = 0.7
  ) +
  geom_text(
    data = pval_labels,
    aes(
      x = 2,
      y = y_position,
      label = pvalue_label
    ),
    size = 1.75,
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ celltype,
    scales = "free_x",
    nrow = 1,
    strip.position = "bottom"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.1)),
    labels = scales::label_percent(scale = 1)
  ) +
  scale_color_manual(values = UBC_disease_colours) +
  scale_fill_manual(values = UBC_disease_colours_lighter) +
  theme(
    panel.grid.major = element_blank(),
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
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 7),
    strip.placement = "outside",
    strip.clip = "off",
    plot.title = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    legend.position = "none"
  )

ggsave(
  plot = plt,
  filename = paste0(
    supplemental_figure_output_path,
    "Supplemental_Figure_4a.svg"
  ),
  width = 2 * nlevels(props_long$celltype),
  height = 4,
  units = "cm"
)
