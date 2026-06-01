# Supplemental Figure 9a: Dot plot of level-2 cell type marker genes across annotated SALTIRE3 scRNA-seq atlas cell types
library(Seurat)
library(ggplot2)
library(scico)
library(openxlsx)
library(dplyr)

source("../../utils.R")

# Load data
atlas_SALTIRE3_sc <- readRDS(
  "../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/merged_objects.rds"
)

# Order of cell types in the plot, using the labels that are actually in the object
celltype_levels <- rev(c(
  "VIC-transitional",    # NTM
  "VIC-contractile",     # MYH11
  "VIC-quiescent",       # RARRES2
  "VIC-FAP-osteogenic",  # BGLAP
  "VIC-neural-crest",    # SOX10
  "VIC-spongiosa",       # GFAP
  "VIC-IFN",             # IFIT1
  "VEC-vascular",        # ACKR1
  "VEC-endocardial",     # BMPER
  "T-gd",                # TRDV2
  "T-CD8",               # CD8B
  "T-CD4",               # IL7R
  "Plasma cell",         # MZB1
  "NK",                  # KLRF1
  "Neutrophil",          # FCGR3B
  "Monocyte",            # FCN1
  "Macrophage",          # C1QC
  "Dendritic cell",      # CD1C
  "Mast cell",           # TPSAB1
  "Proliferating",       # MKI67
  "B cell"               ))

# Labels to display on the plot
ct_labels <- c(
  "Proliferating"        = "Proliferating",
  "Monocyte"             = "Monocytes",
  "Macrophage"           = "Macrophages",
  "Mast cell"            = "Mast cells",
  "VIC-FAP-osteogenic"   = "VIC-FAP-osteogenic",
  "Neutrophil"           = "Neutrophils",
  "Dendritic cell"       = "Dendritic cells",
  "VEC-vascular"         = "EC-vascular",
  "NK"                   = "NK cells",
  "T-CD4"                = "CD4+ T cells",
  "B cell"               = "B cells",
  "Plasma cell"          = "Plasma cells",
  "T-gd"                 = "γδ T cells",
  "T-CD8"                = "CD8+ T cells",
  "VIC-contractile"      = "VIC-contractile",
  "VEC-endocardial"      = "EC-endocardial",
  "VIC-transitional"     = "VIC-transitional",
  "VIC-quiescent"        = "VIC-quiescent",
  "VIC-spongiosa"        = "VIC-spongiosa",
  "VIC-IFN"              = "VIC-IFN",
  "VIC-neural-crest"     = "VIC-neural-crest"
)

# Set plotting order
atlas_SALTIRE3_sc$annotations_level2_readable <- factor(
  atlas_SALTIRE3_sc$annotations_level2_readable,
  levels = celltype_levels
)

markers <- c(
  # VIC-transitional
  "NTM",
  # VIC-contractile
  "MYH11",
  # VIC-quiescent
  "RARRES2",
  # VIC-FAP-osteogenic
  "BGLAP",
  # VIC-neural-crest
  "SOX10",
  # VIC-spongiosa
  "GFAP",
  # VIC-IFN
  "IFIT1",
  # VEC-vascular
  "ACKR1",
  # VEC-endocardial
  "BMPER",
  # T-gd
  "TRDV2",
  # T-CD8
  "CD8B",
  # T-CD4
  "IL7R",
  # Plasma cell
  "MZB1",
  # NK
  "KLRF1",
  # Neutrophil
  "FCGR3B",
  # Monocyte
  "FCN1",
  # Macrophage
  "C1QC",
  # Dendritic cell
  "CD1C",
  # Mast cell
  "TPSAB1",
  # Proliferating
  "MKI67",
  # B cell
  "CD79A"
)

# Build once so the same object can be used for the plot and source data
dotplot_obj <- DotPlot(
  atlas_SALTIRE3_sc,
  features = markers,
  group.by = "annotations_level2_readable",
  col.min = 0,
  col.max = 1.5,
  dot.scale = 2
)

plt <- dotplot_obj +
  scale_colour_gradientn(
    colours = scico(10, palette = "oslo", direction = -1)
  ) +
  coord_equal() +
  scale_y_discrete(
    limits = celltype_levels,
    labels = ct_labels
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 5, colour = "black"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      margin = margin(t = 1)
    ),
    axis.text.y = element_text(
      size = 5,
      colour = "black",
      margin = margin(r = 1)
    ),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    panel.border = element_rect(
      linewidth = 0.7,
      linetype = "solid",
      colour = "black"
    ),
    legend.position = "none"
  )

ggsave(
  plot = plt,
  filename = paste0(
    supplemental_figure_output_path,
    "Supplemental_Figure_9a.svg"
  ),
  width = 6,
  height = 6,
  units = "cm"
)

write.csv(plt$data, file = paste0(source_data_output_path, "Supplemental_Figure_9a.csv"), row.names = TRUE)
