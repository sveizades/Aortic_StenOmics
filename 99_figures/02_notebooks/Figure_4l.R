# Figure 4l: CellPhoneDB interactions to VIC-FAP-osteogenic cells

library(Seurat)
library(ggplot2)
library(ktplots)
library(ggpubr)
source("../../utils.R")

# --- Load data ---
cpdb_path <- "../../03_SALTIRE3_scRNA_seq/03_output/11_S3_scRNAseq_CellPhoneDB/results/statistical_analysis_"
cpdb_suffix <- "_05_09_2026_092951.txt"
read_cpdb <- function(type) read.delim(paste0(cpdb_path, type, cpdb_suffix), check.names = FALSE)

pvals    <- read_cpdb("pvalues")
means    <- read_cpdb("means")
sig_means <- read_cpdb("significant_means")
decon    <- read_cpdb("deconvoluted")

atlas_sce <- as.SingleCellExperiment(
  readRDS("../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/merged_objects.rds")
)

# --- CellPhoneDB plot & extract MP-macro → VIC-osteogenic interactions ---
cpdb_data <- plot_cpdb(
  scdata       = atlas_sce,
  cell_type1   = "VIC-FAP-osteogenic",
  cell_type2   = ".",
  celltype_key = "annotations_level2_readable",
  means        = means,
  pvals        = pvals,
  title        = "interacting interactions!"
)$data |>
  transform(pvals = replace(pvals, is.na(pvals), 1)) |>
  transform(neglog10_p = -log10(pmax(pvals, 1e-300)))

mp_data <- subset(cpdb_data, Var2 == "Macrophage-VIC-FAP-osteogenic")

# --- Selected interactions ---
top_MPmacro_to_VICosteo <- c(
  "IL1B-IL1-receptor", "IL1A-IL1-receptor",        # Inflammatory priming
  "TNF-TNFRSF1A", "TNFSF12-TNFRSF12A",
  
  "TGFB1-TGFbeta-receptor1", "TGFB1-TGFbeta-receptor2", # Pro-fibrotic / osteogenic priming
  
  "PDGFB-PDGFRA", "PDGFC-PDGFRA", "OSM-OSMR",      # Stromal expansion / inflammatory remodeling
  
  "WNT2B-FZD1-LRP6", "WNT2B-FZD7-LRP6",            # Canonical WNT / osteogenic differentiation
  
  "SPP1-integrin-a9b1-complex",                     # Osteopontin / calcific matrix signaling
  "PLAU-PLAUR",                                     # Matrix remodeling
  "VEGFA-NRP1",                                     # Angiogenic remodeling
  
  "SEMA4A-PLXND1", "SEMA4D-PLXNB2",                # Spatial organisation / migratory cues
  
  "Cholesterol-byLIPA-RORA"                         # Metabolic differentiation
)

# Readable labels for plotting
interaction_labels <- c(
  "IL1B-IL1-receptor"              = "IL-1β–IL-1R",
  "IL1A-IL1-receptor"              = "IL-1α–IL-1R",
  "TNF-TNFRSF1A"                   = "TNF–TNFRSF1A",
  "TNFSF12-TNFRSF12A"              = "TWEAK–Fn14",
  
  "TGFB1-TGFbeta-receptor1"        = "TGF-β1–TGFβR1",
  "TGFB1-TGFbeta-receptor2"        = "TGF-β1–TGFβR2",
  
  "PDGFB-PDGFRA"                   = "PDGF-B–PDGFRα",
  "PDGFC-PDGFRA"                   = "PDGF-C–PDGFRα",
  "OSM-OSMR"                       = "OSM–OSMR",
  
  "WNT2B-FZD1-LRP6"                = "WNT2B–FZD1/LRP6",
  "WNT2B-FZD7-LRP6"                = "WNT2B–FZD7/LRP6",
  
  "SPP1-integrin-a9b1-complex"     = "SPP1–integrin α9β1",
  "PLAU-PLAUR"                     = "uPA–uPAR",
  "VEGFA-NRP1"                     = "VEGF-A–NRP1",
  
  "SEMA4A-PLXND1"                  = "SEMA4A–PLXND1",
  "SEMA4D-PLXNB2"                  = "SEMA4D–PLXNB2",
  
  "Cholesterol-byLIPA-RORA"        = "Cholesterol–RORα"
)


# --- Plot ---
plt_data <- subset(mp_data, Var1 %in% top_MPmacro_to_VICosteo) |>
  transform(
    interaction_readable = interaction_labels[Var1]
  ) |>
  transform(
    interaction_readable = factor(
      interaction_readable,
      levels = interaction_readable[order(scaled_means)]
    )
  )

plt <- ggplot(plt_data, aes(x = interaction_readable, y = scaled_means)) +
  geom_bar(
    stat = "identity",
    fill = "#66C2A5",
    colour = "black",
    linewidth = 0.3
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = c(0, 0.5, 1)
  ) +
  coord_flip() +
  theme_pubr() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    axis.title.y       = element_blank(),
    axis.title.x       = element_blank(),
    axis.text.x        = element_text(colour = "black", size = 5),
    axis.text.y        = element_text(
      colour = "black",
      size = 5,
      margin = margin(r = 1)
    ),
    axis.line          = element_line(colour = "black"),
    strip.background   = element_blank(),
    strip.text         = element_blank(),
    plot.title         = element_blank(),
    plot.margin        = unit(c(1, 1, 1, 1), "mm"),
    legend.position    = "none",
    axis.ticks         = element_blank()
  )
ggsave(paste0(main_figure_output_path, "Figure_4l.svg"), plt, width = 3.8, height = 3.3, units = "cm")
write.csv(plt_data, file = paste0(source_data_output_path, "Figure_4l.csv"), row.names = TRUE)
