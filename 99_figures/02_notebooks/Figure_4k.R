# Figure 4k: CellPhoneDB interactions from VIC-FAP-osteogenic cells

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

# --- CellPhoneDB plot & extract Macrophage → VIC-osteogenic interactions ---
cpdb_data <- plot_cpdb(
  scdata       = atlas_sce,
  cell_type1   = "Macrophage",
  cell_type2   = ".",
  celltype_key = "annotations_level2_readable",
  means        = means,
  pvals        = pvals,
  title        = "interacting interactions!"
)$data |>
  transform(pvals = replace(pvals, is.na(pvals), 1)) |>
  transform(neglog10_p = -log10(pmax(pvals, 1e-300)))

vic_data <- subset(cpdb_data, Var2 == "VIC-FAP-osteogenic-Macrophage")


top18_VICosteo_to_MPmacro <- c(
  "CLU-TREM2-receptor",
  "APOE-TREM2-receptor",
  "CD44-TYROBP",                         # Macrophage activation / lipid-associated state
  
  "IL34-CSF1R",
  "CSF1-CSF1R",                          # Macrophage survival / maintenance
  
  "CXCL12-CXCR4",
  "CX3CL1-CX3CR1",                       # Recruitment / retention
  
  "TGFB3-TGFbeta-receptor1",             # Profibrotic signalling
  
  "CD47-SIRPA",                          # Phagocytic checkpoint
  
  "ANXA1-FPR1",
  "GAS6-MERTK",
  "GAS6-AXL",
  "LGALS3-MERTK",                        # Resolution / efferocytosis
  
  "C3-C3AR1",
  "C5-C5AR1",                            # Complement signalling
  
  "TNFSF10-TNFRSF10B",                   # TNF-family signalling
  "TNFSF11-TNFRSF11A",
  
  "VCAM1-integrin-a4b1-complex"          # Adhesion / immune retention
)

interaction_labels <- c(
  "CLU-TREM2-receptor"               = "CLU–TREM2",
  "APOE-TREM2-receptor"              = "APOE–TREM2",
  "CD44-TYROBP"                      = "CD44–TYROBP",
  
  "IL34-CSF1R"                       = "IL-34–CSF1R",
  "CSF1-CSF1R"                       = "CSF-1–CSF1R",
  
  "CXCL12-CXCR4"                     = "CXCL12–CXCR4",
  "CX3CL1-CX3CR1"                    = "CX3CL1–CX3CR1",
  
  "TGFB3-TGFbeta-receptor1"          = "TGF-β3–TGFβR1",
  
  "CD47-SIRPA"                       = "CD47–SIRPα",
  
  "ANXA1-FPR1"                       = "Annexin A1–FPR1",
  "GAS6-MERTK"                       = "GAS6–MERTK",
  "GAS6-AXL"                         = "GAS6–AXL",
  "LGALS3-MERTK"                     = "Galectin-3–MERTK",
  
  "C3-C3AR1"                         = "C3–C3aR",
  "C5-C5AR1"                         = "C5–C5aR1",
  
  "TNFSF10-TNFRSF10B"                = "TRAIL–DR5",
  "TNFSF11-TNFRSF11A" = "RANK-RANKL",
  
  "VCAM1-integrin-a4b1-complex"      = "VCAM-1–integrin α4β1"
)

# --- Plot data ---
plt_data <- subset(vic_data, Var1 %in% top18_VICosteo_to_MPmacro) |>
  transform(
    interaction_readable = unname(interaction_labels[Var1])
  ) |>
  transform(
    interaction_readable = factor(
      interaction_readable,
      levels = interaction_readable[order(scaled_means)]
    )
  )

# --- Plot ---
plt <- ggplot(plt_data, aes(x = interaction_readable, y = scaled_means)) +
  geom_bar(
    stat = "identity",
    fill = "#e7298a",
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
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    
    axis.text.x = element_text(
      colour = "black",
      size = 5
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 5,
      margin = margin(r = 1)
    ),
    
    axis.line  = element_line(colour = "black"),
    axis.ticks = element_blank(),
    
    strip.background = element_blank(),
    strip.text       = element_blank(),
    plot.title       = element_blank(),
    
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    legend.position = "none"
  )
ggsave(paste0(main_figure_output_path, "Figure_4k.svg"), plt, width = 3.8, height = 3.3, units = "cm")
write.csv(plt_data, file = paste0(source_data_output_path, "Figure_4k.csv"), row.names = TRUE)
