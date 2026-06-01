# prepare_data.R
# Run once before launching the Shiny app to extract and save minimal data files.
# Input:  06_data_deposit/03_output/01_data_deposition_object_preparation/*.rds
# Output: 07_shiny_app/01_data/*.rds

library(Seurat)
library(Matrix)
library(dplyr)

data_path <- "../06_data_deposit/03_output/01_data_deposition_object_preparation/"
out_path  <- "01_data/"
dir.create(out_path, showWarnings = FALSE)

# ---- Curated gene sets ----

snrna_atlas_default_genes <- c(
  # FAP+ osteogenic VIC (key finding)
  "FAP", "RUNX2", "POSTN", "COL1A1", "SPP1", "IBSP",
  # Contractile VIC
  "ACTA2", "TAGLN", "CNN1",
  # Quiescent VIC
  "THY1", "CD34", "DPT",
  # Transitional VIC
  "FN1", "VIM", "VCAN",
  # Spongiosa VIC
  "ACAN", "COMP",
  # Neural crest-like VIC
  "NGFR", "S100B",
  # IFN-stimulated VIC
  "ISG15", "IFIT1",
  # Broad cell type markers
  "PECAM1", "CDH5",          # endothelial
  "CD3D", "CD8A",            # T cell
  "CD68", "MRC1",            # myeloid
  "MS4A1", "CD79A",          # B cell
  "KLRD1",                   # NK cell
  "TPSAB1",                  # mast cell
  "ADIPOQ"                   # adipocyte
)

snrna_vics_default_genes <- c(
  # FAP+ osteogenic
  "FAP", "RUNX2", "POSTN", "COL1A1", "SPP1", "IBSP", "MMP13", "BGLAP",
  # Contractile
  "ACTA2", "TAGLN", "CNN1", "MYH11", "MYLK",
  # Quiescent
  "THY1", "CD34", "DPT", "PDGFRA",
  # Transitional
  "FN1", "VIM", "VCAN", "TGFBI",
  # Spongiosa
  "ACAN", "COMP", "SOX9",
  # Neural crest-like
  "NGFR", "S100B", "SOX10",
  # IFN-stimulated
  "ISG15", "IFIT1", "IFIT2", "MX1"
)

scrna_atlas_default_genes <- c(
  # FAP+ osteogenic VIC
  "FAP", "RUNX2", "POSTN", "COL1A1", "SPP1", "IBSP", "MMP13",
  # Other VIC
  "ACTA2", "TAGLN", "CNN1", "THY1", "VIM", "ACAN", "COMP", "CDH11",
  # Endothelial
  "PECAM1", "CDH5", "VWF", "ACKR1", "PROX1", "LYVE1", "SELE", "SOX17", "RGCC",
  # T / NK
  "CD3D", "CD4", "CD8A", "KLRD1", "NKG7", "GZMB", "GZMK", "FOXP3", "IL7R",
  # Myeloid
  "CD68", "MRC1", "CCL4",
  # B cell
  "MS4A1", "CD79A",
  # Mast
  "TPSAB1"
)

# ---- Helper: select genes to store ----
# Keeps all genes expressed in at least min_frac of cells, plus any curated
# genes that don't meet the threshold. This ensures the app can look up any
# biologically meaningful gene, not just variable features.
expressed_genes <- function(expr_mat, default_genes, min_frac = 0.01) {
  frac    <- Matrix::rowMeans(expr_mat > 0)
  passing <- names(frac[frac >= min_frac])
  sort(union(passing, intersect(default_genes, rownames(expr_mat))))
}

# ---- 1. snRNA-seq atlas (UBC) ----
message("Processing snRNA-seq atlas...")
obj <- readRDS(paste0(data_path, "snRNAseq_atlas.rds"))

meta <- data.frame(
  UMAP_1                  = Embeddings(obj, "umap")[, 1],
  UMAP_2                  = Embeddings(obj, "umap")[, 2],
  annotations_level1      = obj$annotations_level1,
  annotations_level2      = obj$annotations_level2,
  clinical_classification = obj$clinical_classification,
  donor_id                = obj$donor_id,
  cell_barcode            = colnames(obj),
  stringsAsFactors = FALSE
)

expr_magic  <- GetAssayData(obj, assay = "MAGIC", layer = "data")
expr_rna    <- GetAssayData(obj, assay = "RNA",   layer = "data")
store_genes <- expressed_genes(expr_rna, snrna_atlas_default_genes)
expr_rna_sub <- expr_rna[store_genes, ]
snrna_atlas_cells <- colnames(obj)
message("  Storing ", length(store_genes), " genes (expressed in ≥1% of cells + curated)")

saveRDS(list(
  meta          = meta,
  expr_magic    = expr_magic,
  expr_rna      = expr_rna_sub,
  genes_magic   = sort(rownames(expr_magic)),
  genes_rna     = store_genes,
  default_genes = snrna_atlas_default_genes
), paste0(out_path, "snrna_atlas.rds"))

rm(obj, meta, expr_magic, expr_rna, expr_rna_sub); gc()
message("  Saved snrna_atlas.rds")

# ---- 2. snRNA-seq VICs (UBC) ----
message("Processing snRNA-seq VICs...")
obj <- readRDS(paste0(data_path, "snRNAseq_vics.rds"))

common_cells <- intersect(colnames(obj), snrna_atlas_cells)
obj <- subset(obj, cells = common_cells)
coords <- Embeddings(obj, "umap")[common_cells, , drop = FALSE]

umap <- data.frame(
  UMAP_1       = coords[, 1],
  UMAP_2       = coords[, 2],
  cell_barcode = common_cells,
  stringsAsFactors = FALSE
)

expr_magic <- GetAssayData(obj, assay = "MAGIC", layer = "data")
expr_magic <- expr_magic[, common_cells, drop = FALSE]
message("  Storing VIC-specific UMAP and MAGIC assay; RNA is shared from snRNA-seq atlas")

saveRDS(list(
  umap          = umap,
  expr_magic    = expr_magic,
  genes_magic   = sort(rownames(expr_magic)),
  default_genes = snrna_vics_default_genes
), paste0(out_path, "snrna_vics.rds"))

rm(obj, common_cells, coords, umap, expr_magic); gc()
message("  Saved snrna_vics.rds")

# ---- 3. scRNA-seq atlas (SALTIRE3) ----
message("Processing scRNA-seq atlas...")
obj <- readRDS(paste0(data_path, "scRNAseq_atlas.rds"))

meta <- data.frame(
  PAGA_1             = Embeddings(obj, "paga")[, 1],
  PAGA_2             = Embeddings(obj, "paga")[, 2],
  annotations_level1 = obj$annotations_level1,
  annotations_level2 = obj$annotations_level2,
  annotations_level3 = obj$annotations_level3,
  donor_id           = obj$donor_id,
  cell_barcode       = colnames(obj),
  stringsAsFactors   = FALSE
)
if ("umap" %in% names(obj@reductions)) {
  meta$UMAP_1 <- Embeddings(obj, "umap")[, 1]
  meta$UMAP_2 <- Embeddings(obj, "umap")[, 2]
}

expr_rna     <- GetAssayData(obj, assay = "RNA", layer = "data")
expr_magic   <- if ("MAGIC" %in% names(obj@assays)) {
  GetAssayData(obj, assay = "MAGIC", layer = "data")
} else {
  NULL
}
store_genes  <- expressed_genes(expr_rna, scrna_atlas_default_genes)
expr_rna_sub <- expr_rna[store_genes, ]
message("  Storing ", length(store_genes), " genes")

subcluster_paths <- c(
  "VIC subclusters"         = paste0(data_path, "vics.rds"),
  "Endothelial subclusters" = paste0(data_path, "endothelial.rds"),
  "T/NK subclusters"       = paste0(data_path, "tnk_cells.rds"),
  "Myeloid subclusters"    = paste0(data_path, "myeloid_cells.rds")
)

subcluster_umaps <- lapply(subcluster_paths, function(path) {
  if (!file.exists(path)) {
    warning("Subcluster object not found: ", path)
    return(NULL)
  }

  sub_obj <- readRDS(path)
  common_cells <- intersect(colnames(sub_obj), colnames(obj))
  if (!("umap" %in% names(sub_obj@reductions)) || length(common_cells) == 0) {
    warning("No usable UMAP/barcodes in: ", path)
    return(NULL)
  }

  coords <- Embeddings(sub_obj, "umap")[common_cells, , drop = FALSE]
  data.frame(
    UMAP_1       = coords[, 1],
    UMAP_2       = coords[, 2],
    cell_barcode = common_cells,
    stringsAsFactors = FALSE
  )
})
subcluster_umaps <- subcluster_umaps[!vapply(subcluster_umaps, is.null, logical(1))]

valid_subclusters <- names(subcluster_paths)[names(subcluster_paths) %in% names(subcluster_umaps)]
subcluster_magic <- Map(function(view_name, path) {
  sub_obj <- readRDS(path)
  common_cells <- subcluster_umaps[[view_name]]$cell_barcode

  if (!("MAGIC" %in% names(sub_obj@assays))) {
    warning("No MAGIC assay found in: ", path)
    return(NULL)
  }

  magic_mat <- GetAssayData(sub_obj, assay = "MAGIC", layer = "data")
  magic_mat[, common_cells, drop = FALSE]
}, valid_subclusters, subcluster_paths[valid_subclusters])
names(subcluster_magic) <- valid_subclusters
subcluster_magic <- subcluster_magic[!vapply(subcluster_magic, is.null, logical(1))]

saveRDS(list(
  meta             = meta,
  subcluster_umaps = subcluster_umaps,
  subcluster_magic = subcluster_magic,
  expr_magic       = expr_magic,
  expr_rna         = expr_rna_sub,
  genes_magic      = if (!is.null(expr_magic)) sort(rownames(expr_magic)) else character(0),
  subcluster_genes_magic = lapply(subcluster_magic, function(x) sort(rownames(x))),
  genes_rna        = store_genes,
  default_genes    = scrna_atlas_default_genes
), paste0(out_path, "scrna_atlas.rds"))

rm(obj, meta, expr_magic, expr_rna, expr_rna_sub, subcluster_umaps, subcluster_magic); gc()
message("  Saved scrna_atlas.rds")

# ---- 4. Xenium spatial atlas ----
message("Processing Xenium spatial atlas...")
obj <- readRDS(paste0(data_path, "spatial_atlas.rds"))

fapi_slide_map <- tapply(obj$slide_id, obj$FAPI_donor, function(x) x[1])
fapi_donors    <- sort(names(fapi_slide_map))

slide_data <- lapply(setNames(fapi_donors, fapi_donors), function(fd) {
  slide <- fapi_slide_map[[fd]]

  # GetTissueCoordinates() returns coordinates in the same order as cells appear
  # in the Seurat object for that slide. The rownames of its output may be
  # unprefixed barcodes (before RenameCells applied the slide_id prefix), so we
  # must NOT use rownames(coords) to index the expression matrix or metadata.
  # Instead, filter the object by slide_id — the same approach used in figure scripts.
  coords     <- GetTissueCoordinates(obj, image = slide)
  slide_mask <- obj$slide_id == slide
  barcodes   <- colnames(obj)[slide_mask]   # prefixed barcodes matching the expression matrix
  meta       <- obj@meta.data[slide_mask, ]

  data.frame(
    x                   = coords$x,
    y                   = coords$y,
    cell_barcode        = barcodes,          # expression-matrix-compatible barcodes
    annotations_level1  = meta$annotations_level1,
    annotations_level2  = meta$annotations_level2,
    annotations_level3  = meta$annotations_level3,
    niche               = meta$niche,
    dist_to_nodule_edge = meta$dist_to_nodule_edge,
    nCount_Xenium       = meta$nCount_Xenium,
    nFeature_Xenium     = meta$nFeature_Xenium,
    donor_id            = meta$donor_id,
    stringsAsFactors    = FALSE
  )
})

expr_xenium  <- GetAssayData(obj, assay = "Xenium",  layer = "data")
expr_protein <- GetAssayData(obj, assay = "Protein", layer = "data")

saveRDS(list(
  slide_data       = slide_data,
  expr_xenium      = expr_xenium,
  expr_protein     = expr_protein,
  xenium_genes     = sort(rownames(expr_xenium)),
  protein_features = sort(rownames(expr_protein)),
  fapi_donors      = fapi_donors
), paste0(out_path, "spatial.rds"))

rm(obj, slide_data, expr_xenium, expr_protein); gc()
message("  Saved spatial.rds")

message("\nData preparation complete. Run shiny::runApp() to launch the app.")
