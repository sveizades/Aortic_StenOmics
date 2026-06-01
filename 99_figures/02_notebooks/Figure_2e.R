# Figure 2e: Bubble plot of GO pathway enrichment across VIC subpopulations (UBC snRNA-seq).
library(openxlsx)
library(ggplot2)
library(scico)
library(tidyr)
library(dplyr)
library(stringr)

source("../../utils.R")

# Load data
vic_cluster_pathways <- read.xlsx(xlsxFile = "../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/pathways_vics_annotations_level2_readable.xlsx", sheet = "Transitional")
vic_cluster_pathways$cluster <- "Transitional"
for (i in c("IFN-stimulated", "Neural crest-like", "Contractile", "FAP+ osteogenic", "Spongiosa", "Quiescent")) {
  tmp_pathway_df <- read.xlsx(xlsxFile = "../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/pathways_vics_annotations_level2_readable.xlsx", sheet = i)
  tmp_pathway_df$cluster <- i
  vic_cluster_pathways <- rbind(vic_cluster_pathways, tmp_pathway_df)
}

# Pathways to plot
selected_pathways_vics <- c("GO:0030199", #collagen fibril organization
                            "GO:0061035", #cartilage development
                            "GO:0002062", #chondrocyte differentiation
                            "GO:0030282", #bone mineralization
                            "GO:0032970", #regulation of actin filament-based process
                            "GO:0030239", #myofibril assembly
                            "GO:0051017", #actin filament bundle assembly
                            "GO:0061572") #actin filament bundle organisation

pathway_all <- vic_cluster_pathways
# Keep only selected IDs and their Descriptions
selected_df <- pathway_all %>%
  filter(ID %in% selected_pathways_vics) %>%
  distinct(ID, Description)

# Join selected Descriptions back in full data (avoids filtering too early)
pathway_all <- pathway_all %>%
  inner_join(selected_df, by = c("ID", "Description")) %>%
  complete(Description, cluster, fill = list(
    pvalue = NA,
    FoldEnrichment = NA,
    p.adjust = NA
  )) %>%
  mutate(
    p_for_plot = ifelse(is.na(p.adjust), 1, ifelse(p.adjust < 1e-300, 1e-300, p.adjust)),
    logp = -log10(p_for_plot),
    is_significant = p_for_plot < 0.05,
    size_val = ifelse(is_significant, logp, 0.01),
    color_val = ifelse(is_significant, FoldEnrichment, NA)
  )
#pathway_all$color_val[pathway_all$color_val>10] <- 10
pathway_all$Description <- str_to_sentence(pathway_all$Description)
pathway_all$Description <- factor(pathway_all$Description, levels = rev(str_to_sentence(selected_df$Description[match(selected_pathways_vics, selected_df$ID)])))
pathway_all$cluster <- factor(pathway_all$cluster, levels = rev(c("IFN-stimulated", "Neural crest-like", "Contractile", "FAP+ osteogenic", "Spongiosa", "Quiescent", "Transitional")))
plt <- ggplot(pathway_all, aes(x = cluster, y = Description)) +
  geom_point(aes(size = size_val, fill = color_val), shape = 21, colour = "white", stroke = 0) +
  scale_fill_gradientn(colours = scico(10, palette = "acton", direction = -1)[1:5], na.value = "black") +
  scale_size(range = c(0.01, 3), breaks = c(1, 2, 3, 4, 5), name = expression(-log[10](p))) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.05),
        axis.text.y = element_text(colour = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 0.7, linetype = "solid", colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.position = "none") +
  coord_equal()
ggsave(plot = plt, filename = paste0(main_figure_output_path, "Figure_2e.svg"), width = 9, height = 3.5, units = "cm")

plt_size_legend <- plt +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 5, colour = "black"),
    legend.title = element_text(size = 5, colour = "black"),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm")
  ) +
  guides(
    fill = "none",
    size = guide_legend(
      override.aes = list(fill = "white", colour = "black", shape = 21, stroke = 0.2)
    )
  )

size_legend <- cowplot::get_legend(plt_size_legend)

ggsave(
  filename = paste0(main_figure_output_path, "Figure_2e_size_legend.svg"),
  plot = size_legend,
  width = 2.2,
  height = 3,
  units = "cm",
  bg = "transparent"
)
