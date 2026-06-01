# Figure 2g: Force-directed layout of VIC trajectories with Slingshot pseudotime curves, alternative view (UBC snRNA-seq).
library(SingleCellExperiment)
library(slingshot)

source("../../utils.R")

# Load data
vic_sce <- readRDS("../../02_UBC_snRNA_seq/03_output/12_UBC_snRNAseq_vic_velocity/vics_sce_slingshot.rds")

# Extract embedding and color info
embedding <- reducedDims(vic_sce)$FA
plot_df <- as.data.frame(embedding)
plot_df$color <- UBC_atlas_colours_vics[vic_sce$annotations_level2_readable]

# Extract Slingshot curves
sds <- SlingshotDataSet(vic_sce)
curves <- slingCurves(sds)

# Convert curves to a dataframe
curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
  cdat <- as.data.frame(curves[[i]]$s)
  cdat$curve <- paste0("curve", i)
  cdat$order <- seq_len(nrow(cdat))
  return(cdat)
}))

# Plot
colnames(plot_df)[1:2] <- c("x", "y")
colnames(curve_df)[1:2] <- c("x", "y")
plt <- ggplot(plot_df, aes(x = x, y = y)) +
  ggrastr::rasterise(geom_point(aes(color = color), size = 0.6, stroke = 0, shape = 16), dpi = 600) +
  geom_path(data = curve_df, aes(x = x, y = y, group = curve), color = "black", size = 0.3) +
  scale_color_identity() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        legend.position = "none", panel.border = element_blank()) +
  coord_equal()
ggsave(plot = plt, filename = paste0(main_figure_output_path, "Figure_2g.svg"), width = 4, height = 4, units = "cm")
