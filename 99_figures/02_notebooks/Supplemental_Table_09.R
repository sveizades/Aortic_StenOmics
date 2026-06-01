# Supplemental Table 9: Xenium differentially expressed genes

library(Seurat)
library(openxlsx)
library(future)
source("../../utils.R")

# Supplemental Table 9: Xenium spatial transcriptomics differentially expressed genes
# FindAllMarkers run on the Xenium atlas (atlas_lim) grouped by annotation level.

atlas_lim <- readRDS("../../04_SALTIRE3_Xenium/03_output/02_SALTIRE3_xenium_histology/atlas_lim.rds")

# Run FindAllMarkers in parallel for each annotation level
plan(multisession)

Idents(atlas_lim) <- "annotations_level1"
xenium_level1 <- FindAllMarkers(atlas_lim, only.pos = TRUE, densify = TRUE, max.cells.per.ident = 2000)

Idents(atlas_lim) <- "annotations_level2"
xenium_level2 <- FindAllMarkers(atlas_lim, only.pos = TRUE, densify = TRUE, max.cells.per.ident = 2000)

Idents(atlas_lim) <- "annotations_level3"
xenium_level3 <- FindAllMarkers(atlas_lim, only.pos = TRUE, densify = TRUE, max.cells.per.ident = 2000)

plan(sequential)

wb <- createWorkbook()

addWorksheet(wb, "Broad cell types (level 1)")
writeData(wb, "Broad cell types (level 1)", xenium_level1)

addWorksheet(wb, "Subclusters (level 2)")
writeData(wb, "Subclusters (level 2)", xenium_level2)

addWorksheet(wb, "Fine-grained subtypes (level 3)")
writeData(wb, "Fine-grained subtypes (level 3)", xenium_level3)

saveWorkbook(wb, file = paste0(supplemental_table_output_path, "Supplemental_Table_09.xlsx"), overwrite = TRUE)
