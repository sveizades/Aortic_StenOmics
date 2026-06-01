# Supplemental Table 6: VIC marker gene GO biological process pathways

library(openxlsx)
source("../../utils.R")

# Supplemental Table 6: Gene ontology biological processes pathway list from VIC sub-cluster marker genes
# GO-BP pathways enriched in the marker genes for each snRNA-seq VIC sub-cluster (clusterProfiler).
# Source workbook contains one sheet per VIC subtype; copied here sheet-by-sheet.

vic_pathway_file <- "../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/pathways_vics_annotations_level2_readable.xlsx"
sheet_names <- getSheetNames(vic_pathway_file)

wb <- createWorkbook()

for (sheet in sheet_names) {
  dat <- read.xlsx(vic_pathway_file, sheet = sheet)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, dat)
}

saveWorkbook(wb, file = paste0(supplemental_table_output_path, "Supplemental_Table_06.xlsx"), overwrite = TRUE)
