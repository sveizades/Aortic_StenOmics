# Supplemental Figure 11: Whole-slide and zoom spatial niche/cell-type maps for BI5XO_A1 Xenium slide (panels 11a and 11b)

library(Seurat)
library(ggplot2)
library(sf)
library(dplyr)
library(ggrastr)
library(openxlsx)

source("../../utils.R")

pixel_size <- 0.2125
grey <- "#2D3132"

# Functions

affine_xy <- function(xy, M) {
  # xy: n x 2
  out <- cbind(xy, 1) %*% t(M)
  out[, 1:2, drop = FALSE]
}

# Apply affine to POLYGON/MULTIPOLYGON while preserving rings
affine_geom <- function(g, M, pixel_size = NULL) {
  g <- sf::st_cast(g, "MULTIPOLYGON")

  xy <- sf::st_coordinates(g)
  xy2 <- affine_xy(xy[, c("X", "Y"), drop = FALSE], M)

  if (!is.null(pixel_size)) {
    xy2 <- xy2 * pixel_size
  }

  # Rebuild MULTIPOLYGON structure using L1 (polygon id) and L2 (ring id)
  df <- data.frame(
    X = xy2[, 1],
    Y = xy2[, 2],
    L1 = xy[, "L1"],
    L2 = xy[, "L2"]
  )

  polys <- lapply(split(df, df$L1), function(p1) {
    rings <- lapply(split(p1, p1$L2), function(r) {
      as.matrix(r[, c("X", "Y"), drop = FALSE])
    })
    rings
  })

  sf::st_multipolygon(polys)
}

make_square_sf <- function(xmin, xmax, ymin, ymax) {
  square_coords <- matrix(
    c(
      xmin, ymin,
      xmax, ymin,
      xmax, ymax,
      xmin, ymax,
      xmin, ymin
    ),
    ncol = 2,
    byrow = TRUE
  )

  square_poly <- st_polygon(list(square_coords))
  st_sfc(square_poly)
}

get_zoom_subset <- function(df, xmin, xmax, ymin, ymax) {
  df |>
    dplyr::filter(
      x >= xmin,
      x <= xmax,
      y >= ymin,
      y <= ymax
    )
}

make_threshold_colour_vector <- function(
    df_subset,
    colour_vector,
    grey_colour = grey,
    min_cells = 30
) {
  cell_counts <- df_subset |>
    dplyr::count(annotations_level2, name = "n_cells")

  low_groups <- cell_counts |>
    dplyr::filter(n_cells < min_cells) |>
    dplyr::pull(annotations_level2)

  colour_vector_full <- colour_vector
  colour_vector_full[low_groups] <- grey_colour

  colour_vector_full
}

make_selected_colour_vector <- function(
    df_subset,
    colour_vector,
    grey_colour = grey
) {
  unlabelled_groups <- setdiff(
    unique(as.character(df_subset$annotations_level2)),
    names(colour_vector)
  )

  c(
    colour_vector,
    setNames(rep(grey_colour, length(unlabelled_groups)), unlabelled_groups)
  )
}

plot_niche_zoom <- function(df_subset) {
  ggplot() +
    ggrastr::geom_point_rast(
      data = df_subset,
      aes(x, y, colour = niche),
      size = 0.05,
      raster.dpi = 600
    ) +
    scale_colour_manual(values = spatial_niche_cols) +
    coord_sf(expand = FALSE) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "black", colour = NA),
      plot.background  = element_rect(fill = "black", colour = NA),
      legend.position = "none"
    )
}

plot_celltype_zoom <- function(df_subset, colour_vector) {
  ggplot() +
    ggrastr::geom_point_rast(
      data = df_subset,
      aes(x, y, colour = annotations_level2),
      size = 0.05,
      raster.dpi = 600
    ) +
    scale_colour_manual(values = colour_vector) +
    coord_sf(expand = FALSE) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "black", colour = NA),
      plot.background  = element_rect(fill = "black", colour = NA),
      legend.position = "none"
    )
}

add_source_sheet <- function(wb, sheet_name, df) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = df)
}

# For Spongiosa, Cap, and Perinodular:
# cell types with >= 30 cells are shown in colour; those with < 30 cells are greyed.
get_threshold_colour_table <- function(
    df_subset,
    colour_vector,
    grey_colour = grey,
    min_cells = 20
) {
  cell_counts <- df_subset |>
    dplyr::count(annotations_level2, name = "n_cells")

  cell_counts |>
    dplyr::mutate(
      labelled_in_colour = n_cells >= min_cells,
      plotted_colour = dplyr::if_else(
        labelled_in_colour,
        unname(colour_vector[as.character(annotations_level2)]),
        grey_colour
      )
    ) |>
    dplyr::arrange(
      dplyr::desc(labelled_in_colour),
      dplyr::desc(n_cells),
      annotations_level2
    )
}

# For Immune:
# only cell types explicitly included in A2_celltype_cols_lim are shown in colour;
# all others are greyed.
get_selected_colour_table <- function(
    df_subset,
    colour_vector,
    grey_colour = grey
) {
  cell_counts <- df_subset |>
    dplyr::count(annotations_level2, name = "n_cells")

  cell_counts |>
    dplyr::mutate(
      labelled_in_colour = annotations_level2 %in% names(colour_vector),
      plotted_colour = dplyr::if_else(
        labelled_in_colour,
        unname(colour_vector[as.character(annotations_level2)]),
        grey_colour
      )
    ) |>
    dplyr::arrange(
      dplyr::desc(labelled_in_colour),
      dplyr::desc(n_cells),
      annotations_level2
    )
}

# Load data

atlas_lim_nodular_tissue <- readRDS(
  "../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds"
)

df_BI5XO_A1 <- data.frame(
  x = GetTissueCoordinates(atlas_lim_nodular_tissue, image = "BI5XO_A1")$x,
  y = GetTissueCoordinates(atlas_lim_nodular_tissue, image = "BI5XO_A1")$y,
  annotations_level1 = atlas_lim_nodular_tissue$annotations_level1[
    atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"
  ],
  annotations_level2 = atlas_lim_nodular_tissue$annotations_level2[
    atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"
  ],
  annotations_level3 = atlas_lim_nodular_tissue$annotations_level3[
    atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"
  ],
  niche = atlas_lim_nodular_tissue$niche[
    atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"
  ],
  dist_to_nodule_edge = atlas_lim_nodular_tissue$dist_to_nodule_edge[
    atlas_lim_nodular_tissue$slide_id == "BI5XO_A1"
  ]
)

# Load and transform nodular tissue annotation

M_BI5XO_A1 <- matrix(
  scan(
    "../../04_SALTIRE3_Xenium/01_data/Post_Xenium_IF_alignment_files_BI5XO_A1/matrix.csv",
    sep = ","
  ),
  nrow = 3,
  byrow = TRUE
)

annotations <- rbind(
  read_sf(
    "../../04_SALTIRE3_Xenium/01_data/E251106C_BI5XO-A10_PH2N_Scan1_all_annotations.geojson"
  ),
  read_sf(
    "../../04_SALTIRE3_Xenium/01_data/E251106C_BXNFB-A10_PH2N_Scan1.geojson"
  )
)

annotations <- annotations[!is.na(annotations$name), ]
annotations$name <- gsub(
  x = annotations$name,
  pattern = "B15XO",
  replacement = "BI5XO"
)

# Extract nodular tissue object
sf_obj_BI5XO_A1 <- annotations[
  annotations$name == "BI5XO_A1 Nodular Tissue",
]

sf::st_crs(sf_obj_BI5XO_A1) <- sf::st_crs(NA)

# Transform nodular tissue into Xenium coordinates
sf_obj_xenium_BI5XO_A1 <- sf_obj_BI5XO_A1

sf::st_geometry(sf_obj_xenium_BI5XO_A1) <- sf::st_sfc(
  lapply(
    sf::st_geometry(sf_obj_BI5XO_A1),
    affine_geom,
    M = M_BI5XO_A1,
    pixel_size = pixel_size
  ),
  crs = sf::st_crs(NA)
)

# Define subpanels and zoom regions

zoom_regions <- tibble::tribble(
  ~panel,   ~region,        ~xmin,               ~xmax,               ~ymin,               ~ymax,               ~colour_rule,
  "11ai",   "Spongiosa",    1252.131990565447,   2154.8678264887396,  3424.253946550486,  4326.9897826257775,  "threshold",
  "11aii",  "Cap",          2026.631073201302,   2565.6931921120527,  4329.6643693121105,  4868.726488222859,  "threshold",
  "11aiii", "Immune",       2918.544,            3280.915,            3717.849,            4080.22,             "selected",
  "11aiv",  "Perinodular",  8155.96726869808,    8746.456835413503,  3588.5104322433743,  4178.999998958796,   "threshold"
)

# White inset boxes shown on Supplemental Figure 11b

square_sf_spongiosa <- make_square_sf(
  xmin = 1252.131990565447,
  xmax = 2154.8678264887396,
  ymin = 3424.253946550486,
  ymax = 4326.9897826257775
)

square_sf_cap <- make_square_sf(
  xmin = 2026.631073201302,
  xmax = 2565.6931921120527,
  ymin = 4329.6643693121105,
  ymax = 4868.726488222859
)

square_sf_immune <- make_square_sf(
  xmin = 2918.544,
  xmax = 3280.915,
  ymin = 3717.849,
  ymax = 4080.22
)

square_sf_perinodular <- make_square_sf(
  xmin = 8155.96726869808,
  xmax = 8746.456835413503,
  ymin = 3588.5104322433743,
  ymax = 4178.999998958796
)

# Supplemental Figure 11b: whole-slide niche plot

plt <- ggplot() +
  ggrastr::geom_point_rast(
    data = df_BI5XO_A1,
    aes(x, y, colour = niche),
    size = 0.3,
    stroke = 0,
    raster.dpi = 600
  ) +
  scale_colour_manual(values = spatial_niche_cols) +
  geom_sf(
    data = square_sf_spongiosa,
    fill = NA,
    colour = "white",
    linewidth = 1
  ) +
  geom_sf(
    data = square_sf_cap,
    fill = NA,
    colour = "white",
    linewidth = 1
  ) +
  geom_sf(
    data = square_sf_immune,
    fill = NA,
    colour = "white",
    linewidth = 1
  ) +
  geom_sf(
    data = square_sf_perinodular,
    fill = NA,
    colour = "white",
    linewidth = 1
  ) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", colour = NA),
    plot.background  = element_rect(fill = "black", colour = NA),
    legend.position = "none"
  )

ggsave(
  filename = file.path(supplemental_figure_output_path, "Supplemental_Figure_11b.svg"),
  plot = plt,
  width = 12,
  height = 6,
  units = "cm"
)

# Colours for immune cell-type plot

A2_celltype_cols_lim <- c(
  "B cell"                = "#F94144",
  "Plasma cell"           = "#1434A4",
  "T-CD4"            = "#43AA8B",
  "T-CD8"            = "#DFFF00",
  "VEC-vascular"     = "#F9844A",
  "VIC-spongiosa" = "#E0B0FF",
  "VIC-quiescent"    = "#90DBF4",
  "VIC-contractile"          = "#800080"
)

# Supplemental Figure 11a zoom panels
# Each subpanel has two plots:
#   top row    = niches
#   bottom row = cell types

for (i in seq_len(nrow(zoom_regions))) {

  region_i <- zoom_regions[i, ]

  df_subset_i <- get_zoom_subset(
    df = df_BI5XO_A1,
    xmin = region_i$xmin,
    xmax = region_i$xmax,
    ymin = region_i$ymin,
    ymax = region_i$ymax
  )

  # Make colour vector for the cell-type plot
  celltype_cols_i <- if (region_i$colour_rule == "selected") {
    make_selected_colour_vector(
      df_subset = df_subset_i,
      colour_vector = A2_celltype_cols_lim,
      grey_colour = grey
    )
  } else {
    make_threshold_colour_vector(
      df_subset = df_subset_i,
      colour_vector = A2_celltype_cols_dark_background,
      grey_colour = grey,
      min_cells = 20
    )
  }

  # Niche plot
  plt_niches <- plot_niche_zoom(df_subset_i)

  ggsave(
    plot = plt_niches,
    filename = file.path(
      supplemental_figure_output_path,
      paste0("Supplemental_Figure_", region_i$panel, "_niches.svg")
    ),
    width = 4,
    height = 4,
    units = "cm"
  )

  # Cell-type plot
  plt_celltypes <- plot_celltype_zoom(
    df_subset = df_subset_i,
    colour_vector = celltype_cols_i
  )

  ggsave(
    plot = plt_celltypes,
    filename = file.path(
      supplemental_figure_output_path,
      paste0("Supplemental_Figure_", region_i$panel, "_celltypes.svg")
    ),
    width = 4,
    height = 4,
    units = "cm"
  )
}

# Source data workbook

wb <- createWorkbook()

# Whole-slide panel source data: Supplemental Figure 11b
source_11b <- df_BI5XO_A1 |>
  dplyr::select(
    x,
    y,
    niche,
    annotations_level2
  )

add_source_sheet(
  wb,
  "11b_whole_slide",
  source_11b
)

# Coordinates of the boxes shown on 11b
source_11b_boxes <- zoom_regions |>
  dplyr::select(
    panel,
    region,
    xmin,
    xmax,
    ymin,
    ymax
  )

add_source_sheet(
  wb,
  "11b_inset_boxes",
  source_11b_boxes
)

# Source data for zoom panels and summary of displayed colours
colour_summary_list <- list()

for (i in seq_len(nrow(zoom_regions))) {

  region_i <- zoom_regions[i, ]

  df_subset_i <- get_zoom_subset(
    df = df_BI5XO_A1,
    xmin = region_i$xmin,
    xmax = region_i$xmax,
    ymin = region_i$ymin,
    ymax = region_i$ymax
  )

  # Same cells underlie both plots for each subpanel,
  # but save separate sheets for the niche and cell-type views.
  niche_source_i <- df_subset_i |>
    dplyr::select(
      x,
      y,
      niche
    )

  celltype_source_i <- df_subset_i |>
    dplyr::select(
      x,
      y,
      annotations_level2
    )

  add_source_sheet(
    wb,
    paste0(region_i$panel, "_niches"),
    niche_source_i
  )

  add_source_sheet(
    wb,
    paste0(region_i$panel, "_celltypes"),
    celltype_source_i
  )

  # Record which cell types were shown in colour vs grey
  colour_summary_i <- if (region_i$colour_rule == "selected") {
    get_selected_colour_table(
      df_subset = df_subset_i,
      colour_vector = A2_celltype_cols_lim,
      grey_colour = grey
    )
  } else {
    get_threshold_colour_table(
      df_subset = df_subset_i,
      colour_vector = A2_celltype_cols_dark_background,
      grey_colour = grey,
      min_cells = 20
    )
  }

  colour_summary_i <- colour_summary_i |>
    dplyr::mutate(
      panel = region_i$panel,
      region = region_i$region,
      colour_rule = region_i$colour_rule,
      .before = 1
    )

  colour_summary_list[[i]] <- colour_summary_i
}

colour_summary <- dplyr::bind_rows(colour_summary_list)

add_source_sheet(
  wb,
  "celltype_colours",
  colour_summary
)

saveWorkbook(
  wb,
  file = file.path(
    source_data_output_path,
    "Supplemental_Figure_11_source_data.xlsx"
  ),
  overwrite = TRUE
)

# Text file listing the cell types shown in colour

labelled_colours_txt <- colour_summary |>
  dplyr::filter(labelled_in_colour) |>
  dplyr::group_by(panel, region) |>
  dplyr::summarise(
    labelled_celltypes = paste0(
      annotations_level2,
      " (",
      plotted_colour,
      ")",
      collapse = "; "
    ),
    .groups = "drop"
  )

txt_lines <- c(
  "Cell types shown in colour rather than grey in Supplemental Figure 11a cell-type plots",
  "",
  "Colouring rules:",
  "- 11ai, 11aii, 11aiv: cell types with >= 20 cells within the plotted region shown in colour; cell types with < 20 cells shown in grey.",
  "- 11aiii: selected immune-panel cell types from A2_celltype_cols_lim shown in colour; all others shown in grey.",
  "",
  apply(labelled_colours_txt, 1, function(x) {
    paste0(
      x[["panel"]],
      " [",
      x[["region"]],
      "]: ",
      x[["labelled_celltypes"]]
    )
  })
)

writeLines(
  txt_lines,
  con = file.path(
    source_data_output_path,
    "Supplemental_Figure_11_labelled_colours.txt"
  )
)
