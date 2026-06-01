# Figure 2d-1: Boxplot of FAP+ osteogenic VIC proportion by disease severity (UBC snRNA-seq).
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(openxlsx)

source("../../utils.R")

# Load data
stats_results <- read.csv("../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/VIC_proportions_statistics.csv")
proportions_data <- read.csv("../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/VIC_proportions_data.csv")

# Select which celltype to plot
selected_celltype <- "FAP+ osteogenic"  # Replace with actual VIC subtype

# Filter data for selected celltype
props_long <- proportions_data %>%
  filter(annotations_level2_readable == selected_celltype) %>%
  mutate(clinical_classification = factor(clinical_classification,
                                          levels = c("control", "mildmoderate", "severe")))

# Get p-value for selected celltype
celltype_pval <- stats_results %>%
  filter(celltype == selected_celltype)

# Calculate max proportion for p-value positioning
max_prop <- max(props_long$proportion)

# Format p-value label
pvalue_label <- ifelse(
  celltype_pval$adj.P.Val < 0.001,
  "p<0.001",
  paste0("p=", round(celltype_pval$adj.P.Val, 3))
)

# Create the plot
plt <- ggplot(props_long, aes(x = clinical_classification, y = proportion)) +
  stat_boxplot(aes(colour = clinical_classification), geom = "errorbar",
               linetype = 1, width = 0.5, alpha = 1) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = clinical_classification, colour = clinical_classification), linewidth = 0.5) +
  geom_jitter(aes(fill = clinical_classification, colour = clinical_classification),
              shape = 16,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              color = "black",
              alpha = 1,
              size = 0.5) +
  # Add p-value as text annotation (no brackets)
  annotate("text",
           x = 2,
           y = max_prop * 1.15,
           label = pvalue_label,
           size = 1.75) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)),
                     labels = scales::label_percent(scale = 1)) +
  scale_color_manual(values = UBC_disease_colours) +
  scale_fill_manual(values = UBC_disease_colours_lighter) +
  theme(panel.grid.major = element_blank(),
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
        strip.text = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "mm"),
        legend.position = "none")

ggsave(plot = plt,
       filename = paste0(main_figure_output_path, "Figure_2d-1.svg"),
       width = 1.5,
       height = 2.5,
       units = "cm")

wb <- createWorkbook()
addWorksheet(wb, "stats_results")
writeData(wb, "stats_results", stats_results)
addWorksheet(wb, "proportions_data")
writeData(wb, "proportions_data", proportions_data)
saveWorkbook(
  wb,
  file.path(source_data_output_path, "Figure_2d_Supplemental_Figure_4a.xlsx"),
  overwrite = TRUE
)
