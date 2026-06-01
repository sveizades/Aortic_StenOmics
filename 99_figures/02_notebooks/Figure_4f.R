# Figure 4f: Xenium VIC sub-state proportions by distance to nodule edge

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

source("../../utils.R")
atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")

breaks <- c(seq(0, 1000, by = 50), Inf)

atlas_lim_nodular_tissue$dist_bin_1000 <- cut(
  atlas_lim_nodular_tissue$dist_to_nodule_edge,
  breaks = breaks,
  include.lowest = TRUE,
  right = FALSE,
  labels = c(
    sprintf("bin%04d", seq(0, 950, by = 50)),
    "bin1000"
  )
)

atlas_lim_nodular_tissue$dist_bin_1000 <- factor(atlas_lim_nodular_tissue$dist_bin_1000, levels = sort(unique(atlas_lim_nodular_tissue$dist_bin_1000)))

vic_levels <- rev(c(
  "VIC-FAP-osteogenic", "VIC-transitional", "VIC-contractile",
  "VIC-spongiosa", "VIC-quiescent", "VIC-neural-crest", "VIC-IFN"
))

# 1) within-sample proportions (sample_id × dist_bin_1000 × VIC subtype)
df_props_sample <- atlas_lim_nodular_tissue@meta.data %>%
  dplyr::filter(
    startsWith(annotations_level2, "VIC"),
    !is.na(dist_bin_1000),
    !is.na(annotations_level2),
    !is.na(sample_id)
  ) %>%
  dplyr::count(sample_id, dist_bin_1000, annotations_level2, name = "n") %>%
  dplyr::group_by(sample_id, dist_bin_1000) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

# 2) average proportions across samples (equal weight per sample)
df_props_avg <- df_props_sample %>%
  dplyr::group_by(dist_bin_1000, annotations_level2) %>%
  dplyr::summarise(
    prop = mean(prop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # "bin0250" -> 250
    dist_mid = as.numeric(stringr::str_extract(as.character(dist_bin_1000), "\\d+")),
    annotations_level2 = factor(annotations_level2, levels = vic_levels)
  )

# 3) smooth within subtype (on the averaged curve)
df_smooth <- df_props_avg %>%
  dplyr::group_by(annotations_level2) %>%
  dplyr::arrange(dist_mid) %>%
  dplyr::mutate(
    prop_smooth = stats::loess(prop ~ dist_mid, span = 0.6)$fitted
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prop_smooth = pmax(prop_smooth, 0)) %>%
  # 4) renormalise so stacks sum to 1 at each distance
  dplyr::group_by(dist_mid) %>%
  dplyr::mutate(prop_smooth = prop_smooth / sum(prop_smooth)) %>%
  dplyr::ungroup()

# 5) plot
plt <- ggplot(df_smooth,
       aes(x = dist_mid, y = prop_smooth, fill = annotations_level2)) +
  geom_area(position = "stack", alpha = 0.95) +
  scale_fill_manual(values = vic_cols_black_background) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1, suffix = "%"),
    expand = c(0, 0)
  )+
  scale_x_continuous(labels = comma, expand = c(0, 0), breaks = c(0, 500, 1000)) +
  ggpubr::theme_pubr() +                              
  theme(
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),axis.ticks.length = unit(.1, "cm"),
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5))
ggsave(plot = plt, 
       filename = paste0(main_figure_output_path, "Figure_4f.svg"), 
       width = 2.75, 
       height = 2.5,
       units = "cm")
write.csv(df_smooth, file = paste0(source_data_output_path, "Figure_4f.csv"), row.names = TRUE)
