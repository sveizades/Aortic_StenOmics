# Figure 4i: Spatial niche density by distance to nodule edge

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

source("../../utils.R")
atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")

df <- data.frame(
  cell = colnames(atlas_lim_nodular_tissue),
  dist = atlas_lim_nodular_tissue$dist_to_nodule_edge,
  niche = atlas_lim_nodular_tissue$niche,
  donor = atlas_lim_nodular_tissue$slide_id
)

make_density_summary <- function(df,
                                 dist_col   = dist,
                                 donor_col  = donor_id,
                                 group_col  = celltype,
                                 n_grid     = 512,
                                 from       = 0,
                                 to         = NULL,
                                 bw         = "nrd0",
                                 min_cells_per_donor = 50) {
  
  df2 <- df %>%
    filter(!is.na({{ dist_col }})) %>%
    filter({{ dist_col }} >= from)
  
  if (!is.null(to)) df2 <- df2 %>% filter({{ dist_col }} <= to)
  
  # common grid for all KDEs
  x_grid <- seq(from,
                ifelse(is.null(to), max(df2 %>% pull({{ dist_col }}), na.rm = TRUE), to),
                length.out = n_grid)
  
  # compute KDE per donor within each group, then interpolate onto common grid
  dens_long <- df2 %>%
    dplyr::group_by({{ group_col }}, {{ donor_col }}) %>%
    dplyr::filter(n() >= min_cells_per_donor) %>%   # <-- important for stability
    group_modify(~{
      x <- .x %>% pull({{ dist_col }})
      # density() needs >=2 points
      if (length(x) < 2) return(tibble())
      d <- density(x, n = n_grid, from = min(x_grid), to = max(x_grid), bw = bw)
      tibble(dist = x_grid, density = approx(d$x, d$y, xout = x_grid, rule = 2)$y)
    }) %>%
    ungroup()
  
  # summarize across donors: mean ± SE at each distance
  dens_sum <- dens_long %>%
    dplyr::group_by({{ group_col }}, dist) %>%
    dplyr::summarise(
      mean = mean(density, na.rm = TRUE),
      se   = sd(density, na.rm = TRUE) / sqrt(sum(!is.na(density))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ymin = pmax(mean - se, 0),
      ymax = mean + se
    )
  
  list(dens_long = dens_long, dens_sum = dens_sum)
}
res <- make_density_summary(
  df,
  dist_col  = dist,
  donor_col = donor,
  group_col = niche,      # or lineage
  n_grid    = 512,
  from      = 0,
  to        = 1000,          # optional: match the paper x-range vibe
  bw        = "nrd0",
  min_cells_per_donor = 50
)

plt <- ggplot(res$dens_sum, aes(x = dist, y = mean, colour = niche, fill = niche)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.5) +
  scale_x_continuous(labels = comma, expand = c(0, 0)) +
  scale_y_continuous(labels = scales::label_number(scale = 1000, accuracy = 1), expand = c(0, 0)) +
  ggpubr::theme_pubr() +
  theme(
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    plot.margin = unit(c(1, 3, 1, 1), "mm"))+
  scale_colour_manual(values = spatial_niche_cols) +
  scale_fill_manual(values = spatial_niche_cols)
ggsave(plot = plt, 
       filename = paste0(main_figure_output_path, "Figure_4i.svg"), 
       width = 2.75, 
       height = 3.5, 
       units = "cm")
write.csv(res$dens_sum, file = paste0(source_data_output_path, "Figure_4i.csv"), row.names = TRUE)
