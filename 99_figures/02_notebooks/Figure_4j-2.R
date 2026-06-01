# Figure 4j-2: Spatial niche and VIC subtype maps for BXNFB_D1

library(Seurat)
library(ggplot2)
library(dplyr)
library(sf)
pixel_size <- 0.2125
grey <- "grey40"
A2_celltype_cols_lim <- c(
  "VIC-spongiosa"  = grey,
  "VIC-FAP-osteogenic"    = "#E7298A",
  "VIC-quiescent"     = grey,
  "VIC-contractile"           = grey,
  "VIC-transitional"  = grey,
  "VIC-neural-crest"  =  grey,
  "VIC-IFN"           = grey
)

spatial_niche_cols_lim <- c(
  "Cap" = grey,
  "Spongiosa" = grey,
  "Peri-nodular" = "#66CC55",
  "Immune" = grey
)

affine_xy <- function(xy, M) {
  # xy: n x 2
  out <- cbind(xy, 1) %*% t(M)
  out[, 1:2, drop = FALSE]
}

# 3) Apply affine to POLYGON/MULTIPOLYGON while preserving rings
affine_geom <- function(g, M, pixel_size = NULL) {
  g <- sf::st_cast(g, "MULTIPOLYGON")
  
  xy <- sf::st_coordinates(g)
  xy2 <- affine_xy(xy[, c("X","Y"), drop = FALSE], M)
  
  if (!is.null(pixel_size)) {
    xy2 <- xy2 * pixel_size
  }
  
  # rebuild MULTIPOLYGON structure using L1 (polygon id) and L2 (ring id)
  df <- data.frame(
    X = xy2[,1], Y = xy2[,2],
    L1 = xy[, "L1"], L2 = xy[, "L2"]
  )
  
  polys <- lapply(split(df, df$L1), function(p1) {
    rings <- lapply(split(p1, p1$L2), function(r) {
      as.matrix(r[, c("X","Y"), drop = FALSE])
    })
    rings
  })
  
  sf::st_multipolygon(polys)
}
source("../../utils.R")
atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")
annotations <- rbind(read_sf("../../04_SALTIRE3_Xenium/01_data/E251106C_BI5XO-A10_PH2N_Scan1_all_annotations.geojson"), read_sf("../../04_SALTIRE3_Xenium/01_data/E251106C_BXNFB-A10_PH2N_Scan1.geojson"))
annotations <- annotations[!is.na(annotations$name),]
df_BXNFB_D1 <- data.frame(x=GetTissueCoordinates(atlas_lim_nodular_tissue, image = "BXNFB_D1")$x, y=GetTissueCoordinates(atlas_lim_nodular_tissue, image = "BXNFB_D1")$y, annotations_level1 = atlas_lim_nodular_tissue$annotations_level1[atlas_lim_nodular_tissue$slide_id == "BXNFB_D1"], annotations_level2 = atlas_lim_nodular_tissue$annotations_level2[atlas_lim_nodular_tissue$slide_id == "BXNFB_D1"], annotations_level3 = atlas_lim_nodular_tissue$annotations_level3[atlas_lim_nodular_tissue$slide_id == "BXNFB_D1"], niche = atlas_lim_nodular_tissue$niche[atlas_lim_nodular_tissue$slide_id == "BXNFB_D1"], dist_to_nodule_edge = atlas_lim_nodular_tissue$dist_to_nodule_edge[atlas_lim_nodular_tissue$slide_id == "BXNFB_D1"])
M_BXNFB_D1 <- matrix(scan(paste0("../../04_SALTIRE3_Xenium/01_data/Post_Xenium_IF_alignment_files_BXNFB_D1/matrix.csv"), sep=","), nrow=3, byrow=TRUE)

#extract nodular tissue object
sf_obj_BXNFB_D1 <- annotations[annotations$name == paste0("BXNFB_D1 Nodular Tissue"),]
sf::st_crs(sf_obj_BXNFB_D1) <- sf::st_crs(NA)

#transform nodular tissue to Xenium coordinate
sf_obj_xenium_BXNFB_D1 <- sf_obj_BXNFB_D1
sf::st_geometry(sf_obj_xenium_BXNFB_D1) <- sf::st_sfc(
  lapply(sf::st_geometry(sf_obj_BXNFB_D1), affine_geom, M = M_BXNFB_D1, pixel_size = pixel_size),
  crs = sf::st_crs(NA)  # keep planar coords
)

xmin <- 2498.9608517339625
xmax <- 3567.292264626808
ymin <- 1228.1014371087194
ymax <- 2375.604174149951

plt <- ggplot() +
  # cells
  ggrastr::geom_point_rast(
    data = df_BXNFB_D1,
    aes(x, y, colour = niche),
    size = 0.1,
    raster.dpi = 600, stroke = 0
  ) +
  scale_colour_manual(values = spatial_niche_cols_lim)+
  # annotation polygon(s)
  geom_sf(
    data = sf_obj_xenium_BXNFB_D1,
    fill = NA,
    color = "white",
    linewidth = 0.1, linetype = "22"
  ) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = NA,
    colour = "red",
    linewidth = 0.3
  ) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", colour = NA),
    plot.background  = element_rect(fill = "black", colour = NA),
    legend.position = "none"
  )
ggsave(plot = plt, 
       filename = paste0(main_figure_output_path, "Figure_4j-FAPIbrightniches.svg"), 
       width = 3, 
       height = 3, 
       units = "cm")


df_subset <- df_BXNFB_D1 |>
  dplyr::filter(
    x >= xmin,
    x <= xmax,
    y >= ymin,
    y <= ymax
  )

df_subset <- df_subset[startsWith(df_subset$annotations_level2, "VIC"),]
plt <- ggplot() +
  ggrastr::geom_point_rast(
    data = df_subset,
    aes(x, y, colour = annotations_level2),
    size = 0.5, stroke = 0,
    raster.dpi = 600
  ) +
  scale_colour_manual(values = A2_celltype_cols_lim)+
  # annotation polygon(s)
  geom_sf(
    data = sf_obj_xenium_BXNFB_D1,
    fill = NA,
    color = "white",
    linewidth = 0.5, linetype = "dashed"
  )+
  coord_sf(
    ylim   = c(1228.1014371087194, 2375.604174149951),
    xlim   = c(2498.9608517339625, 3567.292264626808),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", colour = NA),
    plot.background  = element_rect(fill = "black", colour = NA),
    legend.position = "none"
  )

ggsave(plot = plt, 
       filename = paste0(main_figure_output_path, "Figure_4j-FAPIbright-niches-zoom.svg"), 
       width = 3, 
       height = 3, 
       units = "cm")
write.csv(df_BXNFB_D1, file = paste0(source_data_output_path, "Figure_4j-2.csv"), row.names = TRUE)


# -------------------------------------------------------------------------
# Full-size non-rasterised FAP+ osteogenic VIC plot
# -------------------------------------------------------------------------

df_BXNFB_D1_vics <- df_BXNFB_D1 |>
  dplyr::filter(startsWith(annotations_level2, "VIC"))

plt <- ggplot() +
  geom_point(
    data = df_BXNFB_D1_vics,
    aes(x, y, colour = annotations_level2),
    size = 0.1,
    stroke = 0
  ) +
  scale_colour_manual(values = A2_celltype_cols_lim) +
  geom_sf(
    data = sf_obj_xenium_BXNFB_D1,
    fill = NA,
    color = "white",
    linewidth = 0.1,
    linetype = "22"
  ) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", colour = NA),
    plot.background  = element_rect(fill = "black", colour = NA),
    legend.position = "none"
  )

ggsave(
  plot = plt,
  filename = paste0(main_figure_output_path, "Figure_4j-FAPIbright-FAP_osteogenic_full_vector.svg"),
  width = 3,
  height = 3,
  units = "cm"
)

