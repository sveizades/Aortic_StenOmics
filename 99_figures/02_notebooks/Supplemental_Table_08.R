# Supplemental Table 8: scRNA-seq differentially expressed genes

library(openxlsx)
source("../../utils.R")

# Supplemental Table 8: SALTIRE3 scRNA-seq differentially expressed genes
# FindAllMarkers output across all cell type annotation levels from the SALTIRE3 scRNA-seq atlas.
# Level 1 DEGs are from the combined atlas (notebook 10); sub-cluster DEGs are from each
# cell-type-specific sub-clustering notebook (03-09) and the re-embedded sub-atlases (notebook 10).

# Broad cell types (level 1) — combined atlas
scrna_level1 <- read.csv(
  "../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/merged_objects_annotations_level1_readable_degs.csv",
  row.names = 1
)
scrna_level1$sig <- NULL

# Sub clusters (level 2) — combined atlas
scrna_level2 <- read.csv(
  "../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/merged_objects_annotations_level2_readable_degs.csv",
  row.names = 1
)
scrna_level2$sig <- NULL

# VIC subtypes — level 3 annotation from VIC sub-clustering (notebook 03)
scrna_vics <- read.csv(
  "../../03_SALTIRE3_scRNA_seq/03_output/03_S3_scRNAseq_vic_annotation/vics_annotations_level3_readable_degs.csv",
  row.names = 1
)
scrna_vics$sig <- NULL

# Endothelial subtypes — level 3 annotation from endothelial sub-clustering (notebook 04)
scrna_endo <- read.csv(
  "../../03_SALTIRE3_scRNA_seq/03_output/04_S3_scRNAseq_endothelial_annotation/endothelial_annotations_level3_readable_degs.csv",
  row.names = 1
)
scrna_endo$sig <- NULL

# T and NK subtypes — level 3 annotation from re-embedded T/NK sub-atlas (notebook 10)
scrna_tnk <- read.csv(
  "../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/tnk_cells_annotations_level3_readable_degs.csv",
  row.names = 1
)
scrna_tnk$sig <- NULL

# Myeloid subtypes — level 3 annotation from re-embedded myeloid sub-atlas (notebook 10)
scrna_myeloid <- read.csv(
  "../../03_SALTIRE3_scRNA_seq/03_output/10_S3_scRNAseq_annotations_combining/myeloid_cells_annotations_level3_readable_degs.csv",
  row.names = 1
)
scrna_myeloid$sig <- NULL

wb <- createWorkbook()

addWorksheet(wb, "Cell type (level 1)")
writeData(wb, "Cell type (level 1)", scrna_level1)

addWorksheet(wb, "Cell type (level 2)")
writeData(wb, "Cell type (level 2)", scrna_level2)

addWorksheet(wb, "VIC subtypes")
writeData(wb, "VIC subtypes", scrna_vics)

addWorksheet(wb, "Endothelial subtypes")
writeData(wb, "Endothelial subtypes", scrna_endo)

addWorksheet(wb, "T and NK subtypes")
writeData(wb, "T and NK subtypes", scrna_tnk)

addWorksheet(wb, "Myeloid subtypes")
writeData(wb, "Myeloid subtypes", scrna_myeloid)

saveWorkbook(wb, file = paste0(supplemental_table_output_path, "Supplemental_Table_08.xlsx"), overwrite = TRUE)
