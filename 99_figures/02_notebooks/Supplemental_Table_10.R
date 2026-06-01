# Supplemental Table 10: Spatial niche differentially expressed genes

library(openxlsx)
source("../../utils.R")

# Supplemental Table 10: Xenium spatial niche differentially expressed genes
# FindAllMarkers output for cells grouped by scNiche spatial niche assignment.
# DEGs were computed in Supplemental_Figure_10k.R and deposited in source_data_output_path.

niche_degs <- read.csv(
  "../03_output/04_source_data/Supplemental_Figure_10k_atlas_lim_nodular_tissue_niche_degs.csv",
  row.names = 1
)
niche_degs$sig <- NULL

wb <- createWorkbook()

addWorksheet(wb, "Spatial niches")
writeData(wb, "Spatial niches", niche_degs)

saveWorkbook(wb, file = paste0(supplemental_table_output_path, "Supplemental_Table_10.xlsx"), overwrite = TRUE)
