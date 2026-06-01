# Figure 2k: Paired dot-line plots of FAP and CNN1 mIF staining by valve region (UBC cohort).
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(openxlsx)

source("../../utils.R")

# Load data
df_long <- read.csv("../../05_statistics/03_output/01_vic_mIF_statistics/mIF_data.csv")
ihc_stats_table <- read.csv("../../05_statistics/03_output/01_vic_mIF_statistics/mIF_region_LMM_pairwise_with_means.csv")

marker_to_plot <- "FAP"

# Prepare data for plotting
plot_df <- df_long %>%
  dplyr::filter(antibody == marker_to_plot) %>%
  dplyr::mutate(
    region   = factor(region, levels = c("Inner Nodule", "Peri-Nodular Margin", "Other Tissue")),
    slide_id = factor(slide_id)
  )

# Build p-value annotation df (ONLY significant)
ymax <- max(plot_df$pct_stained, na.rm = TRUE)
pvals_plot <- ihc_stats_table %>%
  filter(marker == marker_to_plot,
         p_adj < 0.05) %>%   # <-- only significant
  transmute(
    group1 = region_1,
    group2 = region_2,
    p_adj
  ) %>%
  arrange(p_adj) %>%
  mutate(
    # Format labels
    p_label = ifelse(
      p_adj < 0.001,
      "p<0.001",
      paste0("p=", sprintf("%.3f", round(p_adj, 3)))
    ),
    # Manual y positions
    y.position = ymax + (seq_len(n()) * 0.15 * ymax)
  )

# Colours
region_cols <- c(
  "Inner Nodule"          = "#9a9945",
  "Peri-Nodular Margin"   = "#c76674",
  "Other Tissue"          = "#8961b3"
)

# Plot FAP
p_fap <- ggplot(plot_df, aes(x = region, y = pct_stained, group = slide_id)) +
  geom_line(colour = "black", linewidth = 0.2) +
  geom_point(aes(colour = region), size = 0.8) +
  scale_colour_manual(values = region_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(y = paste0(marker_to_plot, "+ percent"), x = NULL) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(
    limits = c(0, 12.5),
    breaks = seq(0, 12, by = 4),
    expand = c(0, 0),
    labels = scales::label_percent(scale = 1)
  ) +
  theme_classic(base_size = 7) +
  theme(
    axis.title.y = element_text(colour = "black", size = 7),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.text.y  = element_text(colour = "black", size = 5),
    axis.line    = element_line(colour = "black"),
    plot.margin  = grid::unit(c(0, 0, 1, 0), "mm"),
    legend.position = "none"
  ) +
  ggpubr::stat_pvalue_manual(
    pvals_plot,
    label = "p_label",
    y.position = "y.position",
    tip.length = 0.01,
    size = 1.75
  )

wb <- createWorkbook()
addWorksheet(wb, "CNN1")
writeData(wb, "CNN1", plot_df)

marker_to_plot <- "CNN1"

# Prepare data for plotting
plot_df <- df_long %>%
  dplyr::filter(antibody == marker_to_plot) %>%
  dplyr::mutate(
    region   = factor(region, levels = c("Inner Nodule", "Peri-Nodular Margin", "Other Tissue")),
    slide_id = factor(slide_id)
  )

# Build p-value annotation df (ONLY significant)
ymax <- max(plot_df$pct_stained, na.rm = TRUE)
pvals_plot <- ihc_stats_table %>%
  filter(marker == marker_to_plot,
         p_adj < 0.05) %>%   # <-- only significant
  transmute(
    group1 = region_1,
    group2 = region_2,
    p_adj
  ) %>%
  arrange(p_adj) %>%
  mutate(
    # Format labels
    p_label = ifelse(
      p_adj < 0.001,
      "p<0.001",
      paste0("p=", sprintf("%.3f", round(p_adj, 3)))
    ),
    # Manual y positions
    y.position = ymax + (seq_len(n()) * 0.15 * ymax)
  )

# Colours
region_cols <- c(
  "Inner Nodule"          = "#9a9945",
  "Peri-Nodular Margin"   = "#c76674",
  "Other Tissue"          = "#8961b3"
)

# Plot CNN1
p_cnn1 <- ggplot(plot_df, aes(x = region, y = pct_stained, group = slide_id)) +
  geom_line(colour = "black", linewidth = 0.2) +
  geom_point(aes(colour = region), size = 0.8) +
  scale_colour_manual(values = region_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(y = paste0(marker_to_plot, "+ percent"), x = NULL) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(
    limits = c(0, 2.5),
    breaks = seq(0, 2, by = 1),
    expand = c(0, 0),
    labels = scales::label_percent(scale = 1)
  ) +
  theme_classic(base_size = 7) +
  theme(
    axis.title.y = element_text(colour = "black", size = 7),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.text.y  = element_text(colour = "black", size = 5),
    axis.line    = element_line(colour = "black"),
    plot.margin  = grid::unit(c(1, 0, 0, 0), "mm"),
    legend.position = "none"
  ) +
  ggpubr::stat_pvalue_manual(
    pvals_plot,
    label = "p_label",
    y.position = "y.position",
    tip.length = 0.01,
    size = 1.75
  )

# Stack and align
out_p <- (p_fap / p_cnn1)

# Save
ggsave(
  filename = paste0(main_figure_output_path, "Figure_2k.svg"),
  plot = out_p,
  width = 3.2, height = 4, units = "cm"
)

addWorksheet(wb, "FAP")
writeData(wb, "FAP", plot_df)
saveWorkbook(
  wb,
  file.path(source_data_output_path, "Figure_2k.xlsx"),
  overwrite = TRUE
)
