# Supplemental Figure 4b: Dumbbell plot comparing VIC subtype t-statistics between progressive and severe-specific disease contrasts.
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

source("../../utils.R")

# Read in the saved statistical results and proportions data
topTable_results_all <- read.csv(
  "../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/VIC_proportions_statistics.csv"
)

severe_results <- read.csv(
  "../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/VIC_proportions_severe_statistics.csv"
)

paired_df <- topTable_results_all %>%
  dplyr::select(celltype, t_progressive = t) %>%
  left_join(
    severe_results %>%
      dplyr::select(celltype, t_severe = t),
    by = "celltype"
  ) %>%
  mutate(delta = t_severe - t_progressive
  ) %>%
  mutate(celltype = factor(celltype,
                                 levels = celltype[order(t_progressive)]))

# Add significance to long format
sig_lookup <- bind_rows(
  topTable_results_all %>%
    dplyr::select(celltype, adj.P.Val) %>%
    mutate(contrast = "Progressive"),
  severe_results %>%
    dplyr::select(celltype, adj.P.Val) %>%
    mutate(contrast = "Severe-specific")
) %>%
  mutate(neg_log10_p = -log10(adj.P.Val))

paired_long <- paired_df %>%
  tidyr::pivot_longer(cols = c(t_progressive, t_severe),
                      names_to = "contrast",
                      values_to = "t") %>%
  mutate(contrast = factor(contrast,
                           levels = c("t_progressive", "t_severe"),
                           labels = c("Progressive", "Severe-specific"))) %>%
  left_join(sig_lookup, by = c("celltype", "contrast"))

paired_long <- paired_long %>%
  mutate(point_colour = case_when(
    contrast == "Progressive"     & adj.P.Val < 0.05 ~ "Progressive",
    contrast == "Severe-specific" & adj.P.Val < 0.05 ~ "Severe-specific",
    TRUE ~ "NS"
  ))

plt <- ggplot(paired_df, aes(y = celltype)) +
  geom_vline(
    xintercept = 0,
    colour = "grey40",
    linewidth = 0.2
  ) +
  geom_segment(
    aes(
      x = t_progressive,
      xend = t_severe,
      y = celltype,
      yend = celltype
    ),
    colour = "black",
    linewidth = 0.3
  ) +
  geom_point(
  data = paired_long,
  aes(
    x = t,
    colour = point_colour,
    shape = contrast,
    size = ifelse(point_colour == "NS", neg_log10_p * 0.6, neg_log10_p)
  ),
  stroke = 0.1
) +
  scale_colour_manual(
    values = c(
      "Progressive"     = "#4393c3",
      "Severe-specific" = "#d6604d",
      "NS"              = "grey80"
    ),
    breaks = c("Progressive", "Severe-specific"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "Progressive"     = 16,
      "Severe-specific" = 17
    ),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(1, 4),
    name = "-log10(adj. p)",
    breaks = c(0.5, 1, 1.5, 2)
  ) +
  scale_x_continuous(
    limits = c(-4, 5),
    breaks = seq(-4, 4, 2),
    expand = c(0, 0)
  ) +
  labs(
    x = "t-statistic",
    y = NULL
  ) +
  ggpubr::theme_pubr() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    legend.position = "none"
  )

ggsave(
  plot = plt,
  filename = paste0(supplemental_figure_output_path, "S305.svg"),
  width = 8,
  height = 3,
  units = "cm"
)

legend_cs_df <- data.frame(
  x = 1,
  y = c("Progressive", "Severe-specific"),
  point_colour = c("Progressive", "Severe-specific"),
  contrast = c("Progressive", "Severe-specific")
)

legend_cs <- ggplot(
  legend_cs_df,
  aes(
    x = x,
    y = y,
    colour = point_colour,
    shape = contrast
  )
) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = c(
      "Progressive"     = "#4393c3",
      "Severe-specific" = "#d6604d"
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "Progressive"     = 16,
      "Severe-specific" = 17
    ),
    name = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.text.x = element_blank()
  )

ggsave(
  plot = legend_cs,
  filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_4b-color-legend.svg"),
  width = 3,
  height = 2,
  units = "cm"
)

size_vals <- c(0.5, 1, 1.5, 2)

legend_size_df <- data.frame(
  x = 1,
  y = factor(size_vals, levels = rev(size_vals)),
  neg_log10_p = size_vals
)

legend_size <- ggplot(
  legend_size_df,
  aes(
    x = x,
    y = y,
    size = neg_log10_p
  )
) +
  geom_point(shape = 16, colour = "black") +
  scale_size_continuous(
    range = c(1, 4),
    breaks = size_vals,
    name = "-log10(adj. p)"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.text.x = element_blank()
  )

ggsave(
  plot = legend_size,
  filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_4b-size-legend.svg"),
  width = 3,
  height = 3,
  units = "cm"
)

write.csv(paired_df, file = paste0(source_data_output_path, "Supplemental_Figure_4b.csv"), row.names = TRUE)
