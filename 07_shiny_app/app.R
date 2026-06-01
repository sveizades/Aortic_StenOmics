library(shiny)
library(ggplot2)
library(dplyr)
library(Matrix)
library(viridis)
library(cowplot)
library(scico)

UBC_atlas_colours <- c("VICs"= '#332288', "Myeloid"='#CC6677', 'B cells'='#DDCC77', 'T cells'='#117733', 'NK cells'='#88CCEE', 'Endothelial cells'='#882255', 'Adipocytes'='#44AA99', 'Mast cells'='#999933', "Proliferating" = "#AA4499")
UBC_atlas_colours_vics <- c(Transitional= '#1B9E77', Quiescent='#D95F02', Spongiosa='#7570B3', "FAP+ osteogenic"='#E7298A', Contractile='#0F52BA', "Neural crest-like"='#E6AB02', "IFN-stimulated"='black')
SALTIRE3_donor_colours <- c(ESJ = "#5c3989", I1Z = "#56ae6a", MIR = "#6973d8", OJC = "#ac9c3d", UCC = "#c169b9", UVN = "#ba543d", WQJ = "#6b93d9", XRU = "#b84873")
A1_celltype_cols <- c(
  "Valvular interstitial cell"        = "#332288",
  "Endothelial cell"        = "#882255",
  "Granulocyte"     = "#44AA99",
  "Myeloid"      = "#CC6677",
  "NK cell"   = "#88CCEE",
  "Proliferating"     = "#AA4499",
  "T cell"              = "#117733",
  "B cell"      = "#DDCC77"
)
vic_cols_black_background <- c(
  "VIC-IFN" = "#7A6E2F",
  "VIC-spongiosa"  = "#958BC2",
  "VIC-FAP-osteogenic" = "#E7298A",
  "VIC-quiescent" = "#F18D2D",
  "VIC-contractile" = "#0D52BA",
  "VIC-transitional" = "#1C9E76",
  "VIC-neural-crest"  = "#8A5A44",
  "other" = "#2D3132")
spatial_niche_cols <- c(
  "Cap" = "#4477AA",
  "Spongiosa" = "#CCBB44",
  "Peri-nodular" = "#66CC55",
  "Immune" = "#AA3377"
)
A2_celltype_cols <- c(
  # VICs
  "VIC-spongiosa"        = "#958BC2",
  "VIC-quiescent"        = "#F18D2D",
  "VIC-transitional"     = "#1C9E76",
  "VIC-contractile"      = "#0D52BA",
  "VIC-FAP-osteogenic"   = "#E7298A",
  "VIC-neural-crest"     = "#8A5A44",
  "VIC-IFN"              = "#111111",
  
  # Endothelium
  "VEC-endocardial"      = "#FF7F00",
  "VEC-vascular"         = "#00A6D6",
  
  # T / NK
  "T-CD4"                = "#0072B2",
  "T-CD8"                = "#E69F00",
  "T-gd"                 = "#56B4E9",
  "NK"                   = "#CC79A7",
  
  # Cycling
  "Proliferating"        = "#AA4499",
  
  # Myeloid
  "Macrophage"           = "#66C2A5",
  "Dendritic cell"       = "#984EA3",
  "Neutrophil"           = "#A6761D",
  "Monocyte"             = "#E41A1C",
  
  # B / plasma / mast
  "B cell"               = "#FB8072",
  "Plasma cell"          = "#80B1D3",
  "Mast cell"            = "#FDB462"
)
A3_celltype_cols <- c(
  
  # ---------------- VICs ----------------
  "VIC-spongiosa"                  = "#958BC2",
  "VIC-quiescent"                  = "#F18D2D",
  "VIC-transitional"               = "#1C9E76",
  "VIC-contractile"                = "#0D52BA",
  "VIC-FAP-osteogenic"             = "#E7298A",
  "VIC-neural-crest"               = "#8A5A44",
  "VIC-IFN"                        = "#111111",
  "VIC-cycling"                    = "#AA4499",
  
  # ---------------- VEC ----------------
  "VEC-arterial"                   = "#1F78B4",
  "VEC-capillary"                  = "#A6CEE3",
  "VEC-venous"                     = "#33A02C",
  "VEC-endocardial"                = "#FF7F00",
  "VEC-lymphatic"                  = "#6A3D9A",
  "VEC-activated"                  = "#B15928",
  
  # ---------------- T cells ----------------
  "T-CD4-naive"                    = "#66C2A5",
  "T-CD4-EM"                       = "#FC8D62",
  "T-CD4-Th1"                      = "#8DA0CB",
  "T-CD4-Th17"                     = "#E78AC3",
  "T-CD4-Treg"                     = "#A6D854",
  "T-CD4-cyto"                     = "#FFD92F",
  "T-CD8-EM"                       = "#E5C494",
  "T-CD8-cyto"                     = "#B3B3B3",
  "T-CD8-MAIT"                     = "#1B9E77",
  "T-gd-V9"                        = "#D95F02",
  "T-gd-cyto"                      = "#7570B3",
  "T-cycling"                      = "#F781BF",
  
  # ---------------- NK ----------------
  "NK-CD56-dim"                    = "#00BFC4",
  "NK-CD56-bright"                 = "#C77CFF",
  "NK-cycling"                     = "#FF61C3",
  
  # ---------------- Myeloid ----------------
  "Classical monocyte"             = "#E41A1C",
  "Non-classical monocyte"         = "#377EB8",
  "Resident macrophage"            = "#4DAF4A",
  "Resident activated macrophage"  = "#984EA3",
  "IFN-stimulated macrophage"      = "#FF7F00",
  "Scar-associated macrophage"     = "#FFFF33",
  "Osteoclast"                     = "#A65628",
  "cDC1"                           = "#F781BF",
  "cDC2"                           = "#999999",
  "moDC"                           = "#66A61E",
  "pDC"                            = "#E6AB02",
  "Neutrophil"                     = "#A6761D",
  "Myeloid-cycling"                = "#1F78B4",
  
  # ---------------- B / plasma / mast ----------------
  "B-naive"                        = "#8DD3C7",
  "B-memory"                       = "#FB8072",
  "Plasma cell"                    = "#80B1D3",
  "Mast cell"                      = "#FDB462"
)

A3_celltype_cols <- c(
  
  # ---------------- VICs ----------------
  "VIC-spongiosa"                  = "#958BC2",
  "VIC-quiescent"                  = "#F18D2D",
  "VIC-transitional"               = "#1C9E76",
  "VIC-contractile"                = "#0D52BA",
  "VIC-FAP-osteogenic"             = "#E7298A",
  "VIC-neural-crest"               = "#8A5A44",
  "VIC-IFN"                        = "#111111",
  "VIC-cycling"                    = "#AA4499",
  
  # ---------------- VEC ----------------
  "VEC-arterial"                   = "#1F78B4",
  "VEC-capillary"                  = "#A6CEE3",
  "VEC-venous"                     = "#33A02C",
  "VEC-endocardial"                = "#FF7F00",
  "VEC-lymphatic"                  = "#6A3D9A",
  "VEC-activated"                  = "#B15928",
  
  # ---------------- T cells ----------------
  "T-CD4-naive"                    = "#66C2A5",
  "T-CD4-EM"                       = "#FC8D62",
  "T-CD4-Th1"                      = "#8DA0CB",
  "T-CD4-Th17"                     = "#E78AC3",
  "T-CD4-Treg"                     = "#A6D854",
  "T-CD4-cyto"                     = "#FFD92F",
  "T-CD8-EM"                       = "#E5C494",
  "T-CD8-cyto"                     = "#B3B3B3",
  "T-CD8-MAIT"                     = "#1B9E77",
  "T-gd-V9"                        = "#D95F02",
  "T-gd-cyto"                      = "#7570B3",
  "T-cycling"                      = "#F781BF",
  
  # ---------------- NK ----------------
  "NK-CD56-dim"                    = "#00BFC4",
  "NK-CD56-bright"                 = "#C77CFF",
  "NK-cycling"                     = "#FF61C3",
  
  # ---------------- Myeloid ----------------
  "Classical monocyte"             = "#E41A1C",
  "Non-classical monocyte"         = "#377EB8",
  "Resident macrophage"            = "#4DAF4A",
  "Resident activated macrophage"  = "#984EA3",
  "IFN-stimulated macrophage"      = "#FF7F00",
  "Scar-associated macrophage"     = "#FFFF33",
  "Osteoclast"                     = "#A65628",
  "cDC1"                           = "#F781BF",
  "cDC2"                           = "#999999",
  "moDC"                           = "#66A61E",
  "pDC"                            = "#E6AB02",
  "Neutrophil"                     = "#A6761D",
  "Myeloid-cycling"                = "#1F78B4",
  
  # ---------------- B / plasma / mast ----------------
  "B-naive"                        = "#8DD3C7",
  "B-memory"                       = "#FB8072",
  "Plasma cell"                    = "#80B1D3",
  "Mast cell"                      = "#FDB462"
)

A2_celltype_cols_dark_background <- c(
  # VICs
  "VIC-spongiosa"        = "#958BC2",
  "VIC-quiescent"        = "#F18D2D",
  "VIC-transitional"     = "#1C9E76",
  "VIC-contractile"      = "#0D52BA",
  "VIC-FAP-osteogenic"   = "#E7298A",
  "VIC-neural-crest"     = "#8A5A44",
  "VIC-IFN"              = "#7A6E2F",
  
  # Endothelium
  "VEC-endocardial"      = "#FF7F00",
  "VEC-vascular"         = "#00A6D6",
  
  # T / NK
  "T-CD4"                = "#0072B2",
  "T-CD8"                = "#E69F00",
  "T-gd"                 = "#56B4E9",
  "NK"                   = "#CC79A7",
  
  # Cycling
  "Proliferating"        = "#AA4499",
  
  # Myeloid
  "Macrophage"           = "#ffd700",
  "Dendritic cell"       = "#984EA3",
  "Neutrophil"           = "#A6761D",
  "Monocyte"             = "#E41A1C",
  
  # B / plasma / mast
  "B cell"               = "#FB8072",
  "Plasma cell"          = "#80B1D3",
  "Mast cell"            = "#FDB462"
)
VIC_celltype_names <- c("VIC-transitional" = "Transitional", 
                        "VIC-quiescent" = "Quiescent", 
                        "VIC-spongiosa" = "Spongiosa", 
                        "VIC-FAP-osteogenic" = "FAP+ osteogenic", 
                        "VIC-contractile" = "Contractile", 
                        "VIC-neural-crest" = "Neural crest-like", 
                        "VIC-IFN" = "IFN-stimulated")
VIC_celltype_levels <- c("VIC-transitional", 
                         "VIC-quiescent", 
                         "VIC-spongiosa", 
                         "VIC-FAP-osteogenic", 
                         "VIC-contractile", 
                         "VIC-neural-crest", 
                         "VIC-IFN-stim")
VIC_solo_celltype_levels <- c("Transitional", 
                              "Quiescent", 
                              "Spongiosa", 
                              "FAP+ osteogenic", 
                              "Contractile", 
                              "Neural crest-like", 
                              "IFN-stimulated")

# ---- Load pre-processed data ----
# Run prepare_data.R once before launching the app.

snrna   <- readRDS("01_data/snrna_atlas.rds")
vics    <- readRDS("01_data/snrna_vics.rds")
scrna   <- readRDS("01_data/scrna_atlas.rds")
spatial <- readRDS("01_data/spatial.rds")

# ---- Colour palettes ----

disease_cols <- c("Control" = "#4D8FCB", "Mild/moderate" = "#E7BF47", "Severe" = "#D95E5B")

# Level-2 palette for snRNA-seq atlas: VIC subtypes + broad non-VIC types
snrna_level2_cols <- c(
  A2_celltype_cols[startsWith(names(A2_celltype_cols), "VIC")],
  "Myeloid"         = "#CC6677",
  "T cell"          = "#117733",
  "B cell"          = "#DDCC77",
  "NK cell"         = "#88CCEE",
  "Endothelial cell"= "#882255",
  "Adipocyte"       = "#FFC107",
  "Granulocyte"     = "#44AA99",
  "Proliferating"   = "#AA4499"
)

# VIC subtype order for violin y-axis
vic_order <- c("FAP+ osteogenic", "Transitional", "Quiescent", "Spongiosa",
               "Contractile", "Neural crest-like", "IFN-stimulated")

# ---- Zoom regions (FAPI_donor keyed) ----
zoom_regions <- list(
  "CureAS OJC dim"                  = list(fd = "CureAS_OJC_dim",    xmin = 2323,    xmax = 3323,    ymin = 503,     ymax = 1503),
  "CureAS OJC bright"               = list(fd = "CureAS_OJC_bright", xmin = 4892,    xmax = 5892,    ymin = 1400,    ymax = 2400),
  "CureAS UCC bright"               = list(fd = "CureAS_UCC_bright", xmin = 2498.96, xmax = 3567.29, ymin = 1228.10, ymax = 2375.60),
  "CureAS UCC dim"                  = list(fd = "CureAS_UCC_dim",    xmin = 6750,    xmax = 9250,    ymin = 6000,    ymax = 8500),
  "CureAS UVN bright"               = list(fd = "CureAS_UVN_bright", xmin = 6813,    xmax = 7813,    ymin = 1818,    ymax = 2818),
  "CureAS UVN dim"                  = list(fd = "CureAS_UVN_dim",    xmin = 3452,    xmax = 4452,    ymin = 4235,    ymax = 5235),
  "CureAS XRU bright: peri-nodular" = list(fd = "CureAS_XRU_bright", xmin = 8155.97, xmax = 8746.46, ymin = 3588.51, ymax = 4179.00),
  "CureAS XRU bright: spongiosa"    = list(fd = "CureAS_XRU_bright", xmin = 1252.13, xmax = 2154.87, ymin = 3424.25, ymax = 4326.99),
  "CureAS XRU bright: cap"          = list(fd = "CureAS_XRU_bright", xmin = 2026.63, xmax = 2565.69, ymin = 4329.66, ymax = 4868.73),
  "CureAS XRU bright: immune"       = list(fd = "CureAS_XRU_bright", xmin = 2918.54, xmax = 3280.92, ymin = 3717.85, ymax = 4080.22),
  "CureAS XRU dim"                  = list(fd = "CureAS_XRU_dim",    xmin = 4746,    xmax = 5746,    ymin = 1035,    ymax = 2035)
)

# Unique patients (strip _bright / _dim suffix from FAPI_donor names)
patients <- sort(unique(sub("_(bright|dim)$", "", spatial$fapi_donors)))

# Return zoom region labels belonging to a given patient
zooms_for_patient <- function(patient) {
  names(zoom_regions)[sapply(zoom_regions, function(z) startsWith(z$fd, patient))]
}

# ---- Shared plot helpers ----

# Sparse matrix expression lookup (returns NULL if gene not found or barcodes mismatch)
get_expr <- function(mat, gene, barcodes) {
  if (is.null(mat) || !(gene %in% rownames(mat))) return(NULL)
  valid <- barcodes %in% colnames(mat)
  if (!any(valid)) return(NULL)
  out        <- rep(NA_real_, length(barcodes))
  out[valid] <- as.numeric(mat[gene, barcodes[valid]])
  out
}

# Light or dark base theme for scatter/UMAP plots
scatter_theme <- function(dark) {
  bg  <- if (dark) "#111827" else "white"
  txt <- if (dark) "#e5e7eb" else "black"
  theme_classic(base_size = 11) +
    theme(
      panel.background  = element_rect(fill = bg, colour = NA),
      plot.background   = element_rect(fill = bg, colour = NA),
      legend.background = element_rect(fill = bg, colour = NA),
      legend.key        = element_rect(fill = bg),
      axis.text         = element_text(colour = txt, size = 9),
      axis.title        = element_text(colour = txt, size = 10),
      axis.line         = element_line(colour = if (dark) "#4b5563" else "black"),
      axis.ticks        = element_line(colour = if (dark) "#4b5563" else "black"),
      legend.text       = element_text(colour = txt, size = 8),
      legend.title      = element_text(colour = txt, size = 9),
      plot.title        = element_text(colour = txt, size = 11, face = "bold")
    )
}

# Void theme for spatial plots.
void_theme <- function(dark, legend_pos = "none") {
  bg  <- if (dark) "black" else "white"
  txt <- if (dark) "#e0e0e0" else "black"
  theme_void() +
    theme(
      panel.background  = element_rect(fill = bg, colour = NA),
      plot.background   = element_rect(fill = bg, colour = NA),
      legend.background = element_rect(fill = bg, colour = NA),
      legend.text       = element_text(colour = txt, size = 8),
      legend.title      = element_text(colour = txt, size = 9, face = "bold"),
      plot.title        = element_text(colour = txt, size = 10, face = "bold"),
      plot.margin       = margin(4, 4, 4, 4),
      legend.position   = legend_pos   # must come after theme_void()
    )
}

# Violin plot helper for single-cell tabs (always light background).
# Violins are vertical (groups on x-axis, expression on y-axis).
# group_col:   column in meta_df to group by (x-axis, left → right order).
# group_order: character vector giving the desired x-axis order.
# fill_pal:    named colour vector for groups (used when split_col is NULL).
# split_col:   optional column to split each violin by (e.g. "clinical_classification").
# split_pal:   named colour vector for the split variable.
render_violin <- function(meta_df, expr_vals, gene,
                          group_col, group_order, fill_pal,
                          split_col = NULL, split_pal = NULL) {
  req(length(expr_vals) == nrow(meta_df))
  df       <- meta_df
  df$expr  <- expr_vals
  df$group <- factor(df[[group_col]], levels = group_order)

  vln_theme <- theme_classic(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.text.x      = element_text(colour = "black", size = 8,
                                      angle = 45, hjust = 1),
      axis.text.y      = element_text(colour = "black", size = 8),
      axis.title       = element_text(colour = "black", size = 9),
      legend.text      = element_text(size = 8),
      legend.title     = element_text(size = 9),
      plot.title       = element_text(size = 10, face = "bold"),
      legend.position  = "right"
    )

  if (!is.null(split_col) && !is.null(split_pal)) {
    df$split_var <- factor(df[[split_col]], levels = names(split_pal))
    p <- ggplot(df, aes(group, expr, fill = split_var)) +
      geom_violin(scale = "width", trim = TRUE,
                  position = position_dodge(0.85), linewidth = 0.15) +
      scale_fill_manual(values = split_pal, name = NULL)
  } else {
    p <- ggplot(df, aes(group, expr, fill = group)) +
      geom_violin(scale = "width", trim = TRUE, linewidth = 0.15) +
      scale_fill_manual(values = fill_pal, na.value = "grey70") +
      guides(fill = "none")
  }

  p + labs(x = NULL, y = gene, title = gene) + vln_theme
}

# ---- CSS ----
app_css <- "
  body { background-color: #f8f5f0; color: #222222; font-family: system-ui, sans-serif; }
  .navbar-default { background-color: #2b2b2b; border-color: #2b2b2b; }
  .navbar-default .navbar-brand { color: #f0f0f0; font-weight: 600; }
  .navbar-default .navbar-nav > li > a { color: #c8c8c8; }
  .navbar-default .navbar-nav > li > a:hover { color: #ffffff; background-color: #3d3d3d; }
  .navbar-default .navbar-nav > .active > a,
  .navbar-default .navbar-nav > .active > a:focus { background-color: #3d3d3d; color: #ffffff; }
  .well { background-color: #eee8de; border: 1px solid #d6cfc4; box-shadow: none; }
  label, .control-label, p, li, h3, h4 { color: #222222; }
  hr { border-color: #d6cfc4; }
  select, input[type=number], input[type=text] {
    background-color: #ffffff; color: #222222; border: 1px solid #c4bdb5; border-radius: 3px; }
  .tab-content { background-color: #f8f5f0; }
  .nav-tabs { border-bottom: 1px solid #d6cfc4; }
  .nav-tabs > li > a { color: #666660; background-color: #eee8de; border-color: #d6cfc4; }
  .nav-tabs > li > a:hover { background-color: #e4ddd2; color: #333; }
  .nav-tabs > li.active > a,
  .nav-tabs > li.active > a:focus { color: #222222; background-color: #f8f5f0; border-bottom-color: #f8f5f0; }
  .shiny-output-error { color: #c0392b; font-size: 12px; }
  .shiny-input-container .radio label, .shiny-input-container .checkbox label { color: #222222; }
  small { color: #777770; }
"

# ---- UI ----

ui <- navbarPage(
  title   = "AorticStenOmics",
  id      = "nav",
  inverse = TRUE,
  header  = tags$head(tags$style(HTML(app_css))),
  footer  = tags$div(
    style = "padding:8px; text-align:center; color:#999990; font-size:11px; background-color:#f8f5f0;",
    "AorticStenOmics companion app | Veizades Craig et al. 2026"
  ),

  # ---- Overview ----
  tabPanel("Overview",
    fluidRow(column(8, offset = 2,
    tags$h2("Cellular and Spatial Atlas of Human Aortic Stenosis"),
    tags$p(
      "This portal accompanies our multi-modal study of human aortic stenosis, ",
      "integrating single-nucleus RNA sequencing, single-cell RNA sequencing, ",
      "Xenium spatial transcriptomics, histology, and matched in vivo ",
      "[⁶⁸Ga]FAPI-46 PET-CT molecular imaging."
    ),
    tags$p(
      "Across 25 single-nucleus RNA-seq valves, 8 single-cell RNA-seq valves ",
      "with matched molecular imaging, and 8 spatially profiled valve segments (4 patients), ",
      "this resource defines the cellular populations, transcriptional programmes, ",
      "and spatial niches associated with aortic stenosis progression. The app enables ",
      "interactive exploration of gene expression, cell type annotations, and spatial ",
      "organisation, including FAP-expressing osteogenic valvular interstitial cells ",
      "enriched at active peri-nodular disease microenvironments."
    ),
    tags$h3("What’s in this data portal?"),
    tags$ul(
      tags$li("Explore gene expression across annotated single-cell and single-nucleus datasets."),
      tags$li("Visualise cell states and disease-associated gene expression using interactive violin plots and UMAPs."),
      tags$li("Inspect Xenium spatial gene expression, niches, and cell type maps from profiled valve sections."),
      tags$li(
        "Access raw and processed datasets at ",
        tags$a("PLACEHOLDER URL", href = "#", target = "_blank"), "."
      ),
      tags$li(
        "Code for this project is available on ",
        tags$a("GitHub", href = "https://github.com/sveizades/Aortic_StenOmics/"), "."
      ),
      tags$li(
        "Questions? Contact ",
        tags$a(
          "Stefan Veizades",
          href = "mailto:s.veizades@sms.ed.ac.uk"
        ),
        "."
      ),
      tags$hr())
    ))
  ),

  # ---- snRNA-seq atlas ----
  tabPanel("snRNA-seq atlas",
    sidebarLayout(
      sidebarPanel(width = 3,
        tags$h4("UBC snRNA-seq atlas"),
        radioButtons("snrna_umap_by", "Colour UMAP by",
          choices = c("Cell type (level 1)" = "lvl1",
                      "Cell type (level 2)" = "lvl2",
                      "Gene expression"     = "gene"),
          selected = "lvl2"),
        selectizeInput("snrna_gene", "Gene (UMAP + violin)", choices = NULL,
                       options = list(placeholder = "e.g. FAP, CNN1...")),
        tags$small(style = "color:#6b7280;",
          "UMAP uses MAGIC if available, otherwise RNA. Violin always uses RNA."),
        hr(),
        selectInput("snrna_vln_level", "Violin grouped by",
          choices = c("Level 1" = "lvl1", "Level 2" = "lvl2"), selected = "lvl2"),
        checkboxInput("snrna_split_disease", "Split by disease severity", value = FALSE)
      ),
      mainPanel(width = 9,
        fluidRow(
          column(6, plotOutput("snrna_umap",   height = "520px")),
          column(6, plotOutput("snrna_violin", height = "520px"))
        )
      )
    )
  ),

  # ---- snRNA-seq VICs ----
  tabPanel("snRNA-seq VICs",
    sidebarLayout(
      sidebarPanel(width = 3,
        tags$h4("UBC snRNA-seq VICs"),
        radioButtons("vics_umap_by", "Colour UMAP by",
          choices = c("VIC sub-state"   = "subtype",
                      "Gene expression" = "gene"),
          selected = "subtype"),
        selectizeInput("vics_gene", "Gene (UMAP + violin)", choices = NULL,
                       options = list(placeholder = "e.g. FAP, RUNX2...")),
        tags$small(style = "color:#6b7280;",
          "Same gene drives both the UMAP expression overlay and the violin."),
        hr(),
        checkboxInput("vics_split_disease", "Split by disease severity", value = FALSE)
      ),
      mainPanel(width = 9,
        fluidRow(
          column(6, plotOutput("vics_umap",   height = "460px")),
          column(6, plotOutput("vics_violin", height = "460px"))
        )
      )
    )
  ),

  # ---- scRNA-seq atlas ----
  tabPanel("scRNA-seq atlas",
    sidebarLayout(
      sidebarPanel(width = 3,
        tags$h4("SALTIRE3 scRNA-seq"),
        selectInput("scrna_view", "View",
          choices = c("Full atlas",
                      "VIC subclusters",
                      "Endothelial subclusters",
                      "T/NK subclusters",
                      "Myeloid subclusters"),
          selected = "Full atlas"),
        radioButtons("scrna_umap_by", "Colour UMAP by",
          choices = c("Cell type (level 1)" = "lvl1",
                      "Cell type (level 2)" = "lvl2",
                      "Cell type (level 3)" = "lvl3",
                      "Gene expression"     = "gene"),
          selected = "lvl2"),
        selectizeInput("scrna_gene", "Gene (UMAP + violin)", choices = NULL,
                       options = list(placeholder = "e.g. FAP, RUNX2...")),
        tags$small(style = "color:#6b7280;",
          "UMAP/PAGA expression uses view-specific MAGIC if available, otherwise RNA. Violin always uses RNA."),
        hr(),
        selectInput("scrna_vln_level", "Violin grouped by",
          choices = c("Level 2" = "l2", "Level 3" = "l3"), selected = "l2")
      ),
      mainPanel(width = 9,
        fluidRow(
          column(6, plotOutput("scrna_paga",   height = "520px")),
          column(6, plotOutput("scrna_violin", height = "520px"))
        )
      )
    )
  ),

  # ---- Spatial ----
  tabPanel("Spatial",
    # Controls bar above plots
    tags$div(
      class = "well",
      style = "margin:8px 15px 6px 15px; padding:8px 15px;",
      fluidRow(
        column(3,
          selectInput("sp_patient", "Patient",
                      choices = patients, selected = "CureAS_UCC")
        ),
        column(3,
          selectInput("sp_zoom", "Zoom region", choices = character(0))
        ),
        column(3,
          selectInput("sp_color_by", "Colour by",
            choices = c("Cell type (level 2)" = "lvl2",
                        "Cell type (level 1)" = "lvl1",
                        "Spatial niche"       = "niche",
                        "Gene"                = "gene"),
            selected = "niche")
        ),
        column(3,
          conditionalPanel("input.sp_color_by == 'gene'",
            selectizeInput("sp_gene", "Gene", choices = NULL,
                           options = list(placeholder = "Type gene name..."))
          )
        )
      )
    ),
    # Side-by-side plot panels
    fluidRow(
      style = "margin:0 6px;",
      column(6,
        tags$div(style = "background-color:black;",
          plotOutput("sp_overview", height = "580px"))
      ),
      column(6,
        tags$div(style = "background-color:black;",
          plotOutput("sp_zoom", height = "580px"))
      )
    ),
    # Single shared legend below both plots
    fluidRow(
      style = "margin:0 6px;",
      column(12,
        tags$div(style = "background-color:black;",
          plotOutput("sp_legend", height = "110px"))
      )
    )
  )
)

# ---- Server ----

server <- function(input, output, session) {

  # Gene selectize: server-side for all tabs
  updateSelectizeInput(session, "snrna_gene", choices = snrna$genes_rna,
                       server = TRUE, selected = "FAP")
  updateSelectizeInput(session, "vics_gene",
                       choices = sort(union(snrna$genes_rna, union(vics$genes_magic, vics$genes_rna))),
                       server = TRUE, selected = "FAP")
  updateSelectizeInput(session, "sp_gene", choices = spatial$xenium_genes, server = TRUE)

  observeEvent(input$scrna_view, {
    view <- input$scrna_view
    default_gene <- switch(view,
      "Endothelial subclusters" = "CDH11",
      "T/NK subclusters"       = "GZMK",
      "Myeloid subclusters"    = "CCL4",
      "FAP"
    )
    magic_genes <- if (view == "Full atlas") {
      scrna$genes_magic
    } else if (!is.null(scrna$subcluster_genes_magic[[view]])) {
      scrna$subcluster_genes_magic[[view]]
    } else {
      character(0)
    }
    choices <- sort(union(scrna$genes_rna, magic_genes))
    updateSelectizeInput(session, "scrna_gene", choices = choices,
                         server = TRUE,
                         selected = if (default_gene %in% choices) default_gene else choices[1])
  }, ignoreInit = FALSE)

  # Spatial: populate zoom choices from the selected patient (covers both bright + dim sections)
  observeEvent(input$sp_patient, {
    zooms <- zooms_for_patient(input$sp_patient)
    updateSelectInput(session, "sp_zoom", choices = zooms,
                      selected = if (length(zooms)) zooms[1] else character(0))
  }, ignoreInit = FALSE)

  # Derive the FAPI_donor (tissue section) from the currently selected zoom region
  sp_section <- reactive({
    req(input$sp_zoom)
    zoom_regions[[input$sp_zoom]]$fd
  })



  # ---- snRNA-seq atlas ----

  snrna_umap_data <- reactive({
    df  <- snrna$meta
    by <- input$snrna_umap_by

    if (by == "gene") {
      gene <- req(input$snrna_gene)
      vals <- get_expr(snrna$expr_magic, gene, df$cell_barcode)
      if (is.null(vals)) vals <- get_expr(snrna$expr_rna, gene, df$cell_barcode)
      req(!is.null(vals))
      df$color_val  <- vals
      df$gene_label <- gene
    } else {
      col <- if (by == "lvl1") "annotations_level1" else "annotations_level2"
      df$color_val <- df[[col]]
    }
    list(df = df, by = by)
  })

  output$snrna_umap <- renderPlot({
    d  <- snrna_umap_data()
    df <- d$df
    if (d$by == "gene") {
      q99 <- quantile(df$color_val, 0.99, na.rm = TRUE)
      ggplot(df, aes(UMAP_1, UMAP_2, colour = color_val)) +
        geom_point(size = 1, stroke = 0, shape = 16) +
        scale_colour_scico(palette = "batlow", direction = -1, name = df$gene_label[1], na.value = "grey30",
                           limits = c(0, q99), oob = scales::squish) +
        labs(x = "UMAP 1", y = "UMAP 2", title = df$gene_label[1]) +
        coord_equal() + scatter_theme(FALSE)
    } else {
      pal <- if (d$by == "lvl1") A1_celltype_cols else snrna_level2_cols
      ggplot(df, aes(UMAP_1, UMAP_2, colour = color_val)) +
        geom_point(size = 1, stroke = 0, shape = 16) +
        scale_colour_manual(values = pal, na.value = "grey70",
                            name = if (d$by == "lvl1") "Cell type" else "Cell type (level 2)") +
        guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1)) +
        labs(x = "UMAP 1", y = "UMAP 2") +
        coord_equal() + scatter_theme(FALSE)
    }
  })

  # Violin for the selected gene, grouped by annotation level, optionally split by disease
  output$snrna_violin <- renderPlot({
    gene <- req(input$snrna_gene)
    vals <- get_expr(snrna$expr_rna, gene, snrna$meta$cell_barcode)
    req(!is.null(vals))

    by_lvl    <- input$snrna_vln_level
    grp_col   <- if (by_lvl == "lvl1") "annotations_level1" else "annotations_level2"
    grp_order <- if (by_lvl == "lvl1") {
      c("Valvular interstitial cell", "Endothelial cell", "Myeloid", "T cell",
        "B cell", "NK cell", "Granulocyte", "Adipocyte", "Proliferating")
    } else {
      c(names(A2_celltype_cols)[startsWith(names(A2_celltype_cols), "VIC")],
        "Myeloid", "T cell", "B cell", "NK cell",
        "Endothelial cell", "Granulocyte", "Adipocyte", "Proliferating")
    }
    fill_pal <- if (by_lvl == "lvl1") A1_celltype_cols else snrna_level2_cols

    if (input$snrna_split_disease) {
      render_violin(snrna$meta, vals, gene, grp_col, grp_order, fill_pal,
                    split_col = "clinical_classification", split_pal = disease_cols)
    } else {
      render_violin(snrna$meta, vals, gene, grp_col, grp_order, fill_pal)
    }
  })

  # ---- snRNA-seq VICs ----

  vics_view_data <- reactive({
    if (!is.null(vics$meta)) {
      meta <- vics$meta
      meta$annotations_level2 <- ifelse(
        as.character(meta$annotations_level2) %in% names(VIC_celltype_names),
        unname(VIC_celltype_names[as.character(meta$annotations_level2)]),
        as.character(meta$annotations_level2)
      )
      return(meta)
    }

    validate(need(!is.null(vics$umap),
                  "Rerun 07_shiny_app/prepare_data.R to create the lean snRNA-seq VIC data."))

    idx  <- match(vics$umap$cell_barcode, snrna$meta$cell_barcode)
    meta <- snrna$meta[idx, , drop = FALSE]
    meta$UMAP_1 <- vics$umap$UMAP_1
    meta$UMAP_2 <- vics$umap$UMAP_2
    meta$annotations_level2 <- ifelse(
      as.character(meta$annotations_level2) %in% names(VIC_celltype_names),
      unname(VIC_celltype_names[as.character(meta$annotations_level2)]),
      as.character(meta$annotations_level2)
    )
    meta
  })

  vics_umap_data <- reactive({
    df <- vics_view_data()
    by <- input$vics_umap_by
    if (by == "gene") {
      gene <- req(input$vics_gene)
      vals <- get_expr(vics$expr_magic, gene, df$cell_barcode)
      if (is.null(vals)) vals <- get_expr(snrna$expr_rna, gene, df$cell_barcode)
      req(!is.null(vals))
      df$color_val  <- vals
      df$gene_label <- gene
    } else {
      df$color_val <- df$annotations_level2
    }
    list(df = df, by = by)
  })

  output$vics_umap <- renderPlot({
    d  <- vics_umap_data()
    df <- d$df
    if (d$by == "gene") {
      q99 <- quantile(df$color_val, 0.99, na.rm = TRUE)
      ggplot(df, aes(UMAP_1, UMAP_2, colour = color_val)) +
        geom_point(size = 1, stroke = 0, shape = 16) +
        scale_colour_scico(palette = "batlow", direction = -1, name = df$gene_label[1], na.value = "grey30",
                           limits = c(0, q99), oob = scales::squish) +
        labs(x = "UMAP 1", y = "UMAP 2", title = df$gene_label[1]) +
        coord_equal() + scatter_theme(FALSE)
    } else {
      ggplot(df, aes(UMAP_1, UMAP_2, colour = color_val)) +
        geom_point(size = 1, stroke = 0, shape = 16) +
        scale_colour_manual(values = UBC_atlas_colours_vics, name = "VIC sub-state") +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        labs(x = "UMAP 1", y = "UMAP 2", title = "VIC sub-states") +
        coord_equal() + scatter_theme(FALSE)
    }
  })

  output$vics_violin <- renderPlot({
    df <- vics_view_data()
    gene <- req(input$vics_gene)
    vals <- get_expr(snrna$expr_rna, gene, df$cell_barcode)
    req(!is.null(vals))

    if (input$vics_split_disease) {
      render_violin(df, vals, gene,
                    "annotations_level2", vic_order, UBC_atlas_colours_vics,
                    split_col = "clinical_classification", split_pal = disease_cols)
    } else {
      render_violin(df, vals, gene,
                    "annotations_level2", vic_order, UBC_atlas_colours_vics)
    }
  })

  # ---- scRNA-seq atlas ----

  scrna_view_data <- reactive({
    view <- req(input$scrna_view)

    if (view == "Full atlas") {
      return(list(
        view       = "Full atlas",
        df         = scrna$meta,
        expr_magic = scrna$expr_magic,
        x_col      = if ("PAGA_1" %in% colnames(scrna$meta)) "PAGA_1" else "UMAP_1",
        y_col      = if ("PAGA_2" %in% colnames(scrna$meta)) "PAGA_2" else "UMAP_2",
        x_lab      = if ("PAGA_1" %in% colnames(scrna$meta)) "PAGA 1" else "UMAP 1",
        y_lab      = if ("PAGA_2" %in% colnames(scrna$meta)) "PAGA 2" else "UMAP 2"
      ))
    }
    coords <- scrna$subcluster_umaps[[view]]
    idx    <- match(coords$cell_barcode, scrna$meta$cell_barcode)
    meta   <- scrna$meta[idx, , drop = FALSE]
    meta$subcluster_UMAP_1 <- coords$UMAP_1
    meta$subcluster_UMAP_2 <- coords$UMAP_2

    list(
      view       = view,
      df         = meta,
      expr_magic = scrna$subcluster_magic[[view]],
      x_col      = "subcluster_UMAP_1",
      y_col      = "subcluster_UMAP_2",
      x_lab      = "UMAP 1",
      y_lab      = "UMAP 2"
    )
  })

  scrna_umap_data <- reactive({
    vd <- scrna_view_data()
    df <- vd$df
    by <- input$scrna_umap_by
    if (by == "gene") {
      gene <- req(input$scrna_gene)
      vals <- get_expr(vd$expr_magic, gene, df$cell_barcode)
      if (is.null(vals)) vals <- get_expr(scrna$expr_rna, gene, df$cell_barcode)
      req(!is.null(vals))
      df$color_val  <- vals
      df$gene_label <- gene
    } else {
      col <- switch(by,
        lvl1 = "annotations_level1", lvl2 = "annotations_level2", lvl3 = "annotations_level3")
      df$color_val <- df[[col]]
    }
    list(df = df, by = by, view = vd$view,
         x_col = vd$x_col, y_col = vd$y_col, x_lab = vd$x_lab, y_lab = vd$y_lab)
  })

  output$scrna_paga <- renderPlot({
    d  <- scrna_umap_data()
    df <- d$df
    x_col <- d$x_col
    y_col <- d$y_col
    x_lab <- d$x_lab
    y_lab <- d$y_lab

    if (d$by == "gene") {
      q99 <- quantile(df$color_val, 0.99, na.rm = TRUE)
      ggplot(df, aes(.data[[x_col]], .data[[y_col]], colour = color_val)) +
        geom_point(size = 1, stroke = 0, shape = 16) +
        scale_colour_scico(palette = "batlow", direction = -1, name = df$gene_label[1], na.value = "grey30",
                           limits = c(0, q99), oob = scales::squish) +
        labs(x = x_lab, y = y_lab, title = paste(d$view, df$gene_label[1], sep = ": ")) +
        coord_equal() + scatter_theme(FALSE)
    } else {
      # scRNA-seq level 2 uses specific subtype names (Macrophage, T-CD4, VEC-endocardial,
      # etc.) matching A2_celltype_cols — not the broad snRNA-seq names (Myeloid, T cell).
      groups <- sort(unique(df$color_val))
      pal <- switch(d$by,
        lvl1 = A1_celltype_cols,
        lvl2 = A2_celltype_cols,
        lvl3 = c(A3_celltype_cols,
                 setNames(scales::hue_pal()(max(1, length(setdiff(groups, names(A3_celltype_cols))))),
                          setdiff(unique(df$annotations_level3), names(A3_celltype_cols))))
      )
      missing_cols <- setdiff(groups, names(pal))
      if (length(missing_cols) > 0) {
        pal <- c(pal, setNames(scales::hue_pal()(length(missing_cols)), missing_cols))
      }
      lbl <- switch(d$by, lvl1 = "Cell type",
                    lvl2 = "Cell type (level 2)", lvl3 = "Cell type (level 3)")
      ggplot(df, aes(.data[[x_col]], .data[[y_col]], colour = color_val)) +
        geom_point(size = 1, stroke = 0, shape = 16) +
        scale_colour_manual(values = pal, na.value = "grey70", name = lbl) +
        guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1)) +
        labs(x = x_lab, y = y_lab, title = d$view) +
        coord_equal() + scatter_theme(FALSE)
    }
  })

  output$scrna_violin <- renderPlot({
    vd <- scrna_view_data()
    df <- vd$df
    gene <- req(input$scrna_gene)
    vals <- get_expr(scrna$expr_rna, gene, df$cell_barcode)
    req(!is.null(vals))

    grp_col <- if (input$scrna_vln_level == "l2") "annotations_level2" else "annotations_level3"
    groups  <- sort(unique(df[[grp_col]]))
    pal     <- if (input$scrna_vln_level == "l2") {
      A2_celltype_cols
    } else {
      c(A3_celltype_cols,
        setNames(rep("grey70", length(groups)), groups)[!(groups %in% names(A3_celltype_cols))])
    }
    missing_cols <- setdiff(groups, names(pal))
    if (length(missing_cols) > 0) {
      pal <- c(pal, setNames(scales::hue_pal()(length(missing_cols)), missing_cols))
    }
    render_violin(df, vals, gene, grp_col, groups, pal)
  })

  # ---- Spatial ----

  # Reactive: colour values for all cells in the section derived from the zoom choice.
  sp_colour_ready <- reactive({
    df <- spatial$slide_data[[sp_section()]]
    by <- input$sp_color_by

    if (by == "gene") {
      gene <- req(input$sp_gene)
      vals <- get_expr(spatial$expr_xenium, gene, df$cell_barcode)
      req(!is.null(vals))
      df$color_val <- vals
      list(df = df, type = "continuous", feature = gene,
           vscale = "batlow", is_scico = TRUE,
           q99 = quantile(vals, 0.99, na.rm = TRUE))

    } else {
      col_field <- switch(by,
        lvl1  = "annotations_level1",
        lvl2  = "annotations_level2",
        niche = "niche"
      )
      palette <- switch(by,
        lvl1  = A1_celltype_cols,
        lvl2  = A2_celltype_cols_dark_background,
        niche = spatial_niche_cols
      )
      df$color_val <- df[[col_field]]
      list(df = df, type = "discrete", feature = by, palette = palette)
    }
  })

  # Helper: add colour scale without legend (legend rendered separately).
  # Gene expression uses scico batlow; continuous metadata keeps viridis plasma.
  attach_sp_scale <- function(p, cr) {
    if (cr$type == "continuous") {
      if (isTRUE(cr$is_scico)) {
        lims <- c(0, cr$q99)
        p + scale_colour_scico(palette = cr$vscale, name = cr$feature, na.value = "grey30",
                               limits = lims, oob = scales::squish)
      } else {
        p + scale_colour_viridis_c(option = cr$vscale, name = cr$feature, na.value = "grey30")
      }
    } else {
      p + scale_colour_manual(values = cr$palette, na.value = "grey50",
                              name = cr$feature, drop = FALSE)
    }
  }

  # Build a zoom-processed colour-ready object: rare types collapsed to "Other cells"
  sp_zoom_colour_ready <- reactive({
    cr  <- sp_colour_ready()
    zr  <- zoom_regions[[req(input$sp_zoom)]]

    df_zoom <- cr$df |>
      filter(x >= zr$xmin, x <= zr$xmax, y >= zr$ymin, y <= zr$ymax)
    req(nrow(df_zoom) > 0)

    cr_zoom <- cr
    # Collapse rare cell types to "Other cells" for annotation layers only;
    # spatial niche categories are always well-represented so skip greying there.
    if (cr$type == "discrete" && cr$feature != "niche") {
      counts     <- table(df_zoom$color_val)
      rare_types <- names(counts[counts < 30])
      if (length(rare_types) > 0) {
        df_zoom$color_val[df_zoom$color_val %in% rare_types] <- "Other cells"
        pal               <- cr$palette[!(names(cr$palette) %in% rare_types)]
        pal["Other cells"] <- "#2D3132"
        cr_zoom$palette    <- pal
      }
    }
    cr_zoom$df <- df_zoom
    cr_zoom
  })

  # Overview: full tissue with red zoom box; no legend (shown below)
  output$sp_overview <- renderPlot({
    cr <- sp_colour_ready()
    zr <- zoom_regions[[req(input$sp_zoom)]]

    p <- ggplot(cr$df, aes(x, y, colour = color_val)) +
      geom_point(size = 1.0, stroke = 0, shape = 16)
    p <- attach_sp_scale(p, cr)
    p + annotate("rect",
          xmin = zr$xmin, xmax = zr$xmax, ymin = zr$ymin, ymax = zr$ymax,
          fill = NA, colour = "red", linewidth = 0.6) +
      scale_y_reverse() +
      coord_fixed(expand = FALSE) +
      labs(title = sp_section()) +
      void_theme(TRUE, "none")
  }, bg = "transparent")

  # Zoom: cropped view; scale_y_reverse() matches overview orientation;
  # ylim uses normal order (ymin, ymax) — reversing the limits in coord_fixed()
  # does NOT flip the axis in ggplot2; the scale must do it.
  output$sp_zoom <- renderPlot({
    cz <- sp_zoom_colour_ready()
    zr <- zoom_regions[[req(input$sp_zoom)]]

    p <- ggplot(cz$df, aes(x, y, colour = color_val)) +
      geom_point(size = 2.5, stroke = 0, shape = 16)
    p <- attach_sp_scale(p, cz)
    p + scale_y_reverse() +
      coord_fixed(
        xlim = c(zr$xmin, zr$xmax),
        ylim = c(zr$ymin, zr$ymax),  # normal order; scale_y_reverse handles orientation
        expand = FALSE
      ) +
      labs(title = input$sp_zoom) +
      void_theme(TRUE, "none")
  }, bg = "transparent")

  # Shared legend below both plots.
  # For discrete colours: uses the FULL section palette (matching the overview),
  # then appends "Other cells" in grey to indicate that rare types in the zoom are collapsed.
  output$sp_legend <- renderPlot({
    cr <- sp_colour_ready()

    leg_theme <- theme_void() + theme(
      legend.position   = "bottom",
      legend.background = element_rect(fill = "black", colour = NA),
      legend.text  = element_text(colour = "#e0e0e0", size = 9),
      legend.title = element_text(colour = "#e0e0e0", size = 10),
      plot.background = element_rect(fill = "black", colour = NA)
    )

    if (cr$type == "continuous") {
      df_leg <- data.frame(v = seq(0, 1, length.out = 200), y = 1)
      bar_guide <- guide_colourbar(direction = "horizontal", title.position = "top",
                                   barwidth = unit(10, "cm"), barheight = unit(0.5, "cm"))
      colour_scale <- if (isTRUE(cr$is_scico)) {
        scale_colour_scico(palette = cr$vscale, name = cr$feature, guide = bar_guide,
                           limits = c(0, cr$q99), oob = scales::squish)
      } else {
        scale_colour_viridis_c(option = cr$vscale, name = cr$feature, guide = bar_guide)
      }
      p_leg  <- ggplot(df_leg, aes(v, y, colour = v)) +
        geom_point(size = 0) +
        colour_scale + leg_theme

    } else {
      # Full section palette: only types present in the data, in palette order
      present_types <- unique(cr$df$color_val)
      pal_full      <- cr$palette[names(cr$palette) %in% present_types]
      # Append "Other cells" indicator for annotation layers; not for spatial niche
      # (niche categories are always fully represented in every zoom region)
      if (cr$feature != "niche") {
        pal_full["Other cells"] <- "#2D3132"
      }

      types  <- names(pal_full)
      df_leg <- data.frame(type = factor(types, levels = types), x = seq_along(types), y = 1)
      p_leg  <- ggplot(df_leg, aes(x, y, colour = type)) +
        geom_point(size = 0) +
        scale_colour_manual(values = pal_full, name = NULL, drop = FALSE) +
        guides(colour = guide_legend(
          direction = "horizontal", nrow = 2,
          override.aes = list(size = 4)
        )) + leg_theme
    }

    cowplot::ggdraw(cowplot::get_legend(p_leg))
  }, bg = "transparent")
}

shinyApp(ui, server)
