# Figure 3g: Scatter plot of FAPI TBRmax versus calcific valve volume with linear regression (SALTIRE3 cohort).
library(ggplot2)
library(dplyr)
library(scales)

source("../../utils.R")

# Load data
FAPI_data <- read.csv("../../05_statistics/01_data/SALTIRE3_raw_data_formatted.csv")
FAPI_data <- FAPI_data[!is.na(FAPI_data$Calcific_Volume), ]
FAPI_stats <- read.csv("../../05_statistics/03_output/02_FAPI_statistics/Figure_3g_CalcificThickening_TBRmax_correlation.csv")

label_p <- function(p) {
  if (p < 0.001) "p<0.001"
  else paste0("p=", round(p, 3))
}
# Use your already-computed trend p-value from fit_trend
p_vmax <- FAPI_stats$p_value
r <- FAPI_stats$Correlation
#p_lab <- label_p(p_trend)

#lab_vmax <- paste0("R²=", round(beta_vmax, 3), "\n", label_p(p_vmax))
label_r2 <- paste0("r=", round(r, 3))
label_p  <- ifelse(p_vmax < 0.001, "p<0.001", paste0("p=", round(p_vmax, 3)))

label_text <- paste(label_r2, label_p, sep = "\n")

#label_text <- bquote(R^2 == .(round(FAPI_stats_R$x,3)) ~ "," ~
#                       .(if(p_vmax < 0.001) "p<0.001" else paste0("p=", round(p_vmax,3))))
x_breaks <- c(0, 10, 28, 100, 280, 1000, 2800)
plt_F306 <- ggplot(FAPI_data, aes(x = Calcific_Volume, y = TBRmax)) +
  geom_point(size = 0.3, color = "black", alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.4, color = "blue", fill = "blue", alpha = 0.1) +
  annotate("text",
           x = min(FAPI_data$Calcific_Volume, na.rm = TRUE),
           y = max(FAPI_data$TBRmax, na.rm = TRUE),
           label = label_text,
           hjust = 0, vjust = 1,
           size = 1.6, color = "black") +
  scale_x_continuous(
    trans = "log1p",
    breaks = x_breaks,
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0.05, 0.03))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.10))) +
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
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    legend.position = "none"
  ) +
  theme(axis.text.x = element_text(colour = "black", size = 5),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.25))

ggsave(
  plot = plt_F306,
  filename = paste0(main_figure_output_path, "Figure_3g.svg"),
  width = 5, height = 4.3, units = "cm"
)
write.csv(FAPI_data[, c("Donor_ID", "Calcific_Volume", "TBRmax")], file = paste0(source_data_output_path, "Figure_3g.csv"), row.names = TRUE)
