# Supplemental Figure 3h: Stacked bar chart of cell type proportions per donor coloured by disease group (UBC snRNA-seq atlas).
library(Seurat)
library(ggplot2)
library(dplyr)

source("../../utils.R")

atlas <- readRDS("../../02_UBC_snRNA_seq/03_output/11_UBC_snRNAseq_annotation_combining/atlas.rds")

atlas$donor <- sapply(strsplit(atlas$donor_id, "_"), `[`, 2)

aggregated_metadata <- atlas@meta.data %>%
  group_by(donor) %>%  # Group by sample_id
  summarize(
    disease = paste(unique(clinical_classification), collapse = ", ")
  )
names(aggregated_metadata)[1] <- "Sample"

sample_proprotions <- as.data.frame(table(atlas$donor, atlas$annotations_level1_readable)); colnames(sample_proprotions) <- c("Sample", "Cluster", "Proportion")
sample_proprotions <- merge(sample_proprotions, aggregated_metadata, by = "Sample")

sample_proprotions <- sample_proprotions %>%
  mutate(Group = paste(disease, sep = " | "))

sample_proprotions <- sample_proprotions %>%
  arrange(Group)

sample_proprotions$Group <- factor(sample_proprotions$Group, levels = c("control", "mildmoderate", "severe"))

df <- sample_proprotions %>%
  group_by(Group) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>%
  ungroup()

p <- ggplot(df, aes(fill = Cluster, y = Proportion, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values = c(UBC_atlas_colours, UBC_disease_colours)) +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "line"),
        plot.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, margin = margin(t = -2)),
        axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

header_bars <- df %>%
  group_by(Group) %>%
  summarize(
    xmin = min(as.numeric(Sample)) - 0.5,
    xmax = max(as.numeric(Sample)) + 0.5
  )

stacked_bar <- p + geom_rect(data = header_bars, aes(xmin = xmin, xmax = xmax, ymin = 1.01, ymax = 1.06, fill = Group), inherit.aes = FALSE)
stacked_bar <- stacked_bar + NoLegend()

ggsave(plot = stacked_bar, filename = paste0(supplemental_figure_output_path, "Supplemental_Figure_3h.svg"), width = 6, height = 3, units = "cm")

write.csv(df, file = paste0(source_data_output_path, "Supplemental_Figure_3h.csv"), row.names = TRUE)
