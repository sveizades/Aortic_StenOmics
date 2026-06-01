# Figure 5b: Baseline [68Ga]FAPI-46 PET TBRmax versus CT calcium progression

library(ggplot2)
library(dplyr)
library(scales)
library(grid)

source("../../utils.R")
follow_up <- read.csv('../../05_statistics/01_data/AorticStenomics_FollowUp_Data.csv')

FAPI_stats <- read.csv("../../05_statistics/03_output/03_FAPI_follow_up/TBRmax_CTcalc_correlation.csv")

label_p <- function(p) {
  if (p < 0.001) "p<0.001"
  else paste0("p=", round(p, 3))
}

p_vmax <- FAPI_stats$p.value
r <- FAPI_stats$estimate
#p_lab <- label_p(p_trend)

#lab_vmax <- paste0("R²=", round(beta_vmax, 3), "\n", label_p(p_vmax))
label_r2 <- paste0("r=", round(r, 3))
label_p  <- ifelse(p_vmax < 0.001, "p<0.001", paste0("p=", round(p_vmax, 3)))

label_text <- paste(label_r2, label_p, sep = "\n")

#label_text <- bquote(R^2 == .(round(FAPI_stats_R$x,3)) ~ "," ~
#                       .(if(p_vmax < 0.001) "p<0.001" else paste0("p=", round(p_vmax,3))))
x_breaks <- c(0, 10, 30, 100, 300, 1000, 3000)
plt_F502 <- ggplot(follow_up, aes(x = Baseline_TBRmax, y = CTCalciumScore_Annualised_difference)) +
  geom_point(size = 0.3, color = "black", alpha = 1) +
  geom_smooth(method = "lm", se = T, linewidth = 0.4, color = "blue", fill = "blue", alpha = 0.1) +
  annotate("text",
           x = min(follow_up$Baseline_TBRmax, na.rm = TRUE),
           y = max(follow_up$CTCalciumScore_Annualised_difference, na.rm = TRUE)*6,
           label = label_text,
           hjust = 0, vjust = 1,
           size = 1.6, color = "black") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10, 100, 1000, 3000),
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
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
  )+
  theme(axis.text.x = element_text(colour = "black", size = 5),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.25))

ggsave(
  plot = plt_F502,
  filename = paste0(main_figure_output_path, "Figure_5b.svg"),
  width = 4, height = 3.8, units = "cm"
)
write.csv(follow_up[, c("Donor_ID", "Baseline_TBRmax", "CTCalciumScore_Annualised_difference")], file = paste0(source_data_output_path, "Figure_5b.csv"), row.names = TRUE)
