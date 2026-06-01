# Supplemental Table 3: Bulk RNA-seq GO biological process pathways

library(openxlsx)
source("../../utils.R")

# Supplemental Table 3: Gene ontology biological processes pathway list from bulk RNA sequencing
# Pathways enriched in marker genes for each disease-severity group (clusterProfiler GO-BP).
# One sheet per group: control, mild/moderate, severe.

control_pathways    <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/control_pathways.csv",    row.names = 1)
mildmoderate_pathways <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/mildmoderate_pathways.csv", row.names = 1)
severe_pathways     <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/severe_pathways.csv",     row.names = 1)

wb <- createWorkbook()

addWorksheet(wb, "Control pathways")
writeData(wb, "Control pathways", control_pathways)

addWorksheet(wb, "Mild-moderate pathways")
writeData(wb, "Mild-moderate pathways", mildmoderate_pathways)

addWorksheet(wb, "Severe pathways")
writeData(wb, "Severe pathways", severe_pathways)

saveWorkbook(wb, file = paste0(supplemental_table_output_path, "Supplemental_Table_03.xlsx"), overwrite = TRUE)
