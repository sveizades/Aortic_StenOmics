# Figure 1e: PCA plot of bulk RNA-seq samples coloured by clinical classification (UBC cohort).
library(ggplot2)

source("../../utils.R")

# Load data
bulk_pca_data <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/pca_data.csv")

write.csv(bulk_pca_data[, c("name", "PC1", "PC2", "clinical_classification")], file = paste0(source_data_output_path, "Figure_1e.csv"), row.names = FALSE)

plt <- ggplot(bulk_pca_data, aes(x = PC1, y = PC2, color = clinical_classification)) +
  geom_point(size = 1, shape = 16) +
  stat_ellipse(geom = "polygon",
               aes(fill = clinical_classification),
               alpha = 0.15, linewidth = 0.1) +
  scale_color_manual(values = UBC_disease_colours) +
  scale_fill_manual(values = UBC_disease_colours) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 5),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        legend.position = "none") +
  scale_x_continuous(breaks = c(-10, 0, 10)) +
  scale_y_continuous(breaks = c(-10, 0, 10)) +
  coord_equal()

ggsave(plot = plt, filename = paste0(main_figure_output_path, "Figure_1e.svg"), width = 4.5, height = 4.5, units = "cm")
