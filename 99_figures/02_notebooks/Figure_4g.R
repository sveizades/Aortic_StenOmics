# Figure 4g: Distance-to-nodule association models for Xenium cell states

library(lme4)
library(purrr)
library(broom.mixed)
library(Seurat)
library(ggplot2)
library(tibble)

source("../../utils.R")
atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")
df_glm <- atlas_lim_nodular_tissue@meta.data %>%
  dplyr::filter(
    !is.na(dist_to_nodule_edge),
    !is.na(annotations_level2),
    !is.na(donor_id)
  ) %>%
  dplyr::mutate(
    dist_scaled = as.numeric(scale(dist_to_nodule_edge))
  )


celltypes <- unique(df_glm$annotations_level2)

assoc_tbl <- map_dfr(celltypes, function(ct) {
  
  df_tmp <- df_glm %>%
    dplyr::mutate(is_ct = annotations_level2 == ct)
  fit <- try(
    glmer(
      is_ct ~ dist_scaled + (1 | donor_id),
      data = df_tmp,
      family = binomial(),
      control = glmerControl(optimizer = "bobyqa")
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) return(NULL)
  
  coef_row <- broom.mixed::tidy(fit, effects = "fixed") %>%
    dplyr::filter(term == "dist_scaled")
  
  tibble(
    annotations_level2 = ct,
    beta = coef_row$estimate,
    OR   = exp(coef_row$estimate),
    p    = coef_row$p.value
  )
}) %>%
  dplyr::mutate(
    p_adj = p.adjust(p, method = "BH")
  )
lvl_map <- atlas_lim_nodular_tissue@meta.data
lvl_map <-lvl_map %>%
  dplyr::filter(!is.na(annotations_level1), !is.na(annotations_level2)) %>%
  dplyr::distinct(annotations_level2, annotations_level1)


assoc_tbl2 <- assoc_tbl %>%
  dplyr::left_join(lvl_map, by = "annotations_level2") %>%
  dplyr::mutate(
    sig = dplyr::case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )
plot_df <- assoc_tbl2 %>%
  dplyr::arrange(beta) %>%
  dplyr::mutate(
    annotations_level2 = factor(annotations_level2, levels = sort(unique(annotations_level2)))
  )

plot_df <- plot_df %>%
  dplyr::mutate(
    outline_linetype = dplyr::case_when(
      p_adj < 0.001 ~ "solid",
      p_adj < 0.01  ~ "dashed",
      p_adj < 0.05  ~ "dotted",
      TRUE          ~ "blank"
    ),
    outline_linewidth = dplyr::case_when(
      p_adj < 0.001 ~ 0.5,
      p_adj < 0.01  ~ 0.5,
      p_adj < 0.05  ~ 0.5,
      TRUE          ~ 0
    )
  )
plt <- ggplot(
  plot_df,
  aes(
    y = annotations_level2,
    x = beta,
    fill = annotations_level1
  )
) +
  scale_fill_manual(values = A1_celltype_cols) +
  geom_col(
    aes(
      linetype = outline_linetype,
      linewidth = 0.2
    ),
    colour = "black",
    width = 0.8
  ) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.2) +
  scale_linetype_identity() +
  scale_linewidth_identity() +
  ggpubr::theme_pubr() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 5),
        axis.text.y = element_text(colour = "black", size = 5),
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "mm"),
        legend.position = 'none', axis.ticks = element_blank())
ggsave(plot = plt, 
       filename = paste0(main_figure_output_path, "Figure_4g.svg"), 
       width = 4, 
       height = 4, 
       units = "cm")
fdr_legend <- data.frame(
  beta = 0,
  annotations_level2 = c("p_adj < 0.001", "p_adj < 0.01", "p_adj < 0.05"),
  annotations_level1 = NA,
  outline_linetype = c("dotted", "dashed", "solid"),
  outline_linewidth = c(0.5, 0.5, 0.5)
)

legend <- ggplot() +
  geom_col(
    data = fdr_legend,
    aes(
      x = beta,
      y = annotations_level2,
      linetype = outline_linetype,
      linewidth = outline_linewidth
    ),
    fill = "white",
    colour = "black",
    width = 0.8
  ) +
  scale_linetype_identity(name = "FDR") +
  scale_linewidth_identity() +
  theme_void()
ggsave(plot = legend, 
       filename = paste0(main_figure_output_path, "Figure_4g-legend.svg"), 
       width = 4, 
       height = 6, 
       units = "cm")
write.csv(plot_df, file = paste0(source_data_output_path, "Figure_4g.csv"), row.names = TRUE)
