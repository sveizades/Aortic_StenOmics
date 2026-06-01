# Figure 4c: Xenium spatial map of VIC sub-states in BI5XO_A1

library(Seurat)
library(ggplot2)
library(sf)
source("../../utils.R")

pixel_size <- 0.2125

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

atlas_lim_nodular_tissue <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")

M_BI5XO_A1 <- matrix(scan(paste0("../../04_SALTIRE3_Xenium/01_data/Post_Xenium_IF_alignment_files_BI5XO_A1/matrix.csv"), sep=","), nrow=3, byrow=TRUE)

annotations <- rbind(read_sf("../../04_SALTIRE3_Xenium/01_data/E251106C_BI5XO-A10_PH2N_Scan1_all_annotations.geojson"), read_sf("../../04_SALTIRE3_Xenium/01_data/E251106C_BXNFB-A10_PH2N_Scan1.geojson"))
annotations <- annotations[!is.na(annotations$name),]
annotations$name <- gsub(x = annotations$name, pattern = "B15XO", replacement = "BI5XO")

#extract nodular tissue object
sf_obj_BI5XO_A1 <- annotations[annotations$name == paste0("BI5XO_A1 Nodular Tissue"),]
sf::st_crs(sf_obj_BI5XO_A1) <- sf::st_crs(NA)

#transform nodular tissue to Xenium coordinate
sf_obj_xenium_BI5XO_A1 <- sf_obj_BI5XO_A1
sf::st_geometry(sf_obj_xenium_BI5XO_A1) <- sf::st_sfc(
  lapply(sf::st_geometry(sf_obj_BI5XO_A1), affine_geom, M = M_BI5XO_A1, pixel_size = pixel_size),
  crs = sf::st_crs(NA)  # keep planar coords
)

df_BI5XO_A1 <- data.frame(x=GetTissueCoordinates(atlas_lim_nodular_tissue, image = "BI5XO_A1")$x, y=GetTissueCoordinates(atlas_lim_nodular_tissue, image = "BI5XO_A1")$y, annotations_level1 = atlas_lim_nodular_tissue$annotations_level1[atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"], annotations_level2 = atlas_lim_nodular_tissue$annotations_level2[atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"], annotations_level3 = atlas_lim_nodular_tissue$annotations_level3[atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"], niche = atlas_lim_nodular_tissue$niche[atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"], dist_to_nodule_edge = atlas_lim_nodular_tissue$dist_to_nodule_edge[atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"])


df_BI5XO_A1_vics <- df_BI5XO_A1
df_BI5XO_A1_vics$annotations_level2[!startsWith(df_BI5XO_A1$annotations_level2, "VIC")] <- "other"

df_BI5XO_A1_vics$annotations_level2 <- factor(
  df_BI5XO_A1_vics$annotations_level2,
  levels = names(vic_cols_black_background)
)

plot <- ggplot() +
  # rasterized cells
  ggrastr::geom_point_rast(
    data = df_BI5XO_A1_vics,
    aes(x, y, colour = annotations_level2),
    size = 0.2, stroke = 0, shape = 16, raster.dpi = 900
  ) +
  scale_colour_manual(values = vic_cols_black_background, drop = FALSE) +
  # vector annotation polygon(s)
  geom_sf(
    data = sf_obj_xenium_BI5XO_A1,
    fill = NA,
    color = "white",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", colour = NA),
    plot.background  = element_rect(fill = "black", colour = NA),
    legend.background = element_rect(fill = "black"),
    legend.text = element_text(colour = "white"),
    legend.position = 'none'
  )

ggsave(
  filename = file.path(main_figure_output_path, "Figure_4c.svg"),
  plot = plot,
  width = 9, height = 5, units = "cm"
)
write.csv(df_BI5XO_A1_vics, file = paste0(source_data_output_path, "Figure_4c.csv"), row.names = TRUE)
