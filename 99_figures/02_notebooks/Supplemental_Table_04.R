# Supplemental Table 4: snRNA-seq broad cluster differentially expressed genes

library(openxlsx)
source("../../utils.R")

# Supplemental Table 4: snRNA-seq broad cluster differentially expressed genes
# FindAllMarkers output from the combined UBC snRNA-seq atlas (notebook 11).
# Level 1: broad cell type labels; Level 2: cell type sub-labels.

snrna_level1 <- read.csv(
  "../../02_UBC_snRNA_seq/03_output/11_UBC_snRNAseq_annotation_combining/atlas_annotations_level1_degs.csv",
  row.names = 1
)

snrna_level1$cluster <- plyr::mapvalues(snrna_level1$cluster, from = c("proliferating", "mast", "adipocytes", "endothelial", "NK_cell", "T_cell", "B_cell", "myeloid", "vics"), to = c("Proliferating", "Mast cells", "Adipocytes", "Endothelial cells", "NK cells", "T cells", "B cells", "Myeloid", "VICs"))
snrna_level1$sig <- NULL

wb <- createWorkbook()

addWorksheet(wb, "Broad cell types (level 1)")
writeData(wb, "Broad cell types (level 1)", snrna_level1)

saveWorkbook(wb, file = paste0(supplemental_table_output_path, "Supplemental_Table_04.xlsx"), overwrite = TRUE)
