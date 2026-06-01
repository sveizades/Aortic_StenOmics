# Supplemental Table 5: snRNA-seq VIC subcluster differentially expressed genes

library(openxlsx)
source("../../utils.R")

# Supplemental Table 5: snRNA-seq valvular interstitial cell sub-cluster differentially expressed genes
# FindAllMarkers output from the UBC snRNA-seq VIC sub-clustering analysis (notebook 04).
# Markers used to define the level-2 (readable) VIC sub-cluster annotations.

snrna_vic_degs <- read.csv(
  "../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/vics_annotations_level2_readable_degs.csv",
  row.names = 1
)

snrna_vic_degs$sig <- NULL

wb <- createWorkbook()

addWorksheet(wb, "VIC subtypes")
writeData(wb, "VIC subtypes", snrna_vic_degs)

saveWorkbook(wb, file = paste0(supplemental_table_output_path, "Supplemental_Table_05.xlsx"), overwrite = TRUE)
