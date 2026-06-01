# Supplemental Table 2: Bulk RNA-seq differentially expressed genes

library(openxlsx)
source("../../utils.R")

# Supplemental Table 2: Bulk RNA sequencing differentially expressed genes
# Four pairwise comparisons from DESeq2:
#   AS vs. no-AS, severe vs. control, mild/moderate vs. control, severe vs. mild/moderate

bulk_DE_AS_v_noAS  <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_AS_v_noAS.csv",  row.names = 1)
bulk_DE_mm_v_ctrl  <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_mm_v_ctrl.csv",  row.names = 1)
bulk_DE_sev_v_ctrl <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_sev_v_ctrl.csv", row.names = 1)
bulk_DE_sev_v_mm   <- read.csv("../../01_bulk_seq/03_output/01_processing_bulkRNAseq/DE_sev_v_mm.csv",   row.names = 1)

wb <- createWorkbook()

addWorksheet(wb, "AS vs no-AS")
writeData(wb, "AS vs no-AS", bulk_DE_AS_v_noAS)

addWorksheet(wb, "Severe vs control")
writeData(wb, "Severe vs control", bulk_DE_sev_v_ctrl)

addWorksheet(wb, "Mild-moderate vs control")
writeData(wb, "Mild-moderate vs control", bulk_DE_mm_v_ctrl)

addWorksheet(wb, "Severe vs mild-moderate")
writeData(wb, "Severe vs mild-moderate", bulk_DE_sev_v_mm)

saveWorkbook(wb, file = paste0(supplemental_table_output_path, "Supplemental_Table_02.xlsx"), overwrite = TRUE)
