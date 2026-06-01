# Figure 1d: Stacked bar chart of pathology grades by disease severity (UBC cohort).
library(readxl)
library(dplyr)
library(ggplot2)
library(grid)

source("../../utils.R")

# Load data
df <- read_excel("../01_data/UBC_cohort_pathology_grading.xlsx", sheet = 1)

# Prepare data
df_plot <- df %>%
  mutate(
    grade_plot = factor(Grade, levels = rev(c("1", "2", "3", "4"))),
    Clinical_Classification = factor(
      Clinical_Classification,
      levels = c("control", "mildmoderate", "severe"),
      labels = c("Control", "Mild/moderate", "Severe")
    )
  ) %>%
  select(Donor_ID, Clinical_Classification, grade_plot)

# Count observations
count_df <- df_plot %>%
  count(Clinical_Classification, grade_plot, name = "n")

# Define palette
pal <- c(
  "1" = "#6A4C93",
  "2" = "#1982C4",
  "3" = "#2CA58D",
  "4" = "#FF6B35"
)

# Plot
plt_grade <- ggplot(
  count_df,
  aes(x = Clinical_Classification,
      y = n,
      fill = grade_plot)
) +
  geom_col(width = 0.55, colour = NA) +
  scale_fill_manual(values = pal) +

  scale_y_continuous(
    limits = c(0, NA),
    expand = c(0, 0)
  ) +

  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black"),

    axis.text.x = element_text(
      colour = "black",
      size = 5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),

    axis.text.y  = element_text(colour = "black", size = 5),
    axis.line    = element_line(colour = "black"),
    plot.margin  = unit(c(1, 1, 1, 1), "mm"),
    legend.position = "none"
  )

# Export
write.csv(
  df_plot,
  file = paste0(source_data_output_path, "Figure_1d.csv"),
  row.names = FALSE
)
ggsave(
  paste0(main_figure_output_path, "Figure_1d.svg"),
  plot = plt_grade,
  width = 2,
  height = 4,
  units = "cm",
  bg = "transparent"
)
