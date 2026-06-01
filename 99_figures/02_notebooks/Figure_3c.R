# Figure 3c: Boxplot of FAPI TBRmax by aortic stenosis disease severity (SALTIRE3 cohort).
library(ggplot2)
library(dplyr)
library(scales)

source("../../utils.R")

# Load data
FAPI_data <- read.csv("../../05_statistics/01_data/SALTIRE3_raw_data_formatted.csv")
FAPI_stats <- read.csv("../../05_statistics/03_output/02_FAPI_statistics/Figure_3c_ASdiagnosis_FAPI_trend_stats.csv")

label_p <- function(p) {
  if (is.na(p)) return("p=NA")
  if (p < 0.001) return("p<0.001")
  paste0("p=", round(p, 3))
}

df <- FAPI_data %>%
  mutate(
    Disease_Severity = factor(
      Disease_Severity,
      levels = c("Control", "Sclerosis", "Mild", "Moderate", "Severe"),
      ordered = TRUE
    ),
    Sex = factor(Sex)
  )

# Use your already-computed trend p-value from fit_trend
p_trend <- FAPI_stats[FAPI_stats$term == "severity_score", "p.value"]
p_lab <- label_p(p_trend)

max_y <- max(df$TBRmax, na.rm = TRUE)

plt_F302 <- ggplot(df, aes(x = Disease_Severity, y = TBRmax)) +
  stat_boxplot(aes(colour = Disease_Severity),
               geom = "errorbar",
               linetype = 1,
               width = 0.5,
               alpha = 1,
               linewidth = 0.4) +
  geom_boxplot(
    outlier.shape = NA,
    aes(fill = Disease_Severity, colour = Disease_Severity),
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(fill = Disease_Severity, colour = Disease_Severity),
    shape = 16,
    position = position_jitter(width = 0.18, height = 0),
    color = "black",
    alpha = 1,
    size = 0.5
  ) +
  annotate(
    "text",
    x = 3,                         # centre-ish (Mild)
    y = max_y,
    label = p_lab,
    size = 1.6
  ) +
  scale_y_continuous(
    #limits = c(0, NA),
    expand = expansion(mult = c(0.10, 0.10)), breaks = c(1.2, 1.4, 1.6)
  ) +
  scale_color_manual(values = SALTIRE_disease_colours) +
  scale_fill_manual(values = SALTIRE_disease_colours_lighter) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),

    axis.title.y = element_blank(),
    axis.title.x = element_blank(),

    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.25),

    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 5),

    axis.line = element_line(colour = "black", linewidth = 0.25),

    strip.background = element_blank(),
    strip.text = element_blank(),

    plot.title = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    legend.position = "none"
  )

ggsave(plot = plt_F302,
       filename = paste0(main_figure_output_path, "Figure_3c.svg"),
       width = 4,
       height = 4,
       units = "cm")
write.csv(df, file = paste0(source_data_output_path, "Figure_3c.csv"), row.names = TRUE)
