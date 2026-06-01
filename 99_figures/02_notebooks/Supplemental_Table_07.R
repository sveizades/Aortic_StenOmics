# Supplemental Table 7: Multivariate analysis of FAPI PET TBRmax

library(openxlsx)
source("../../utils.R")

# Supplemental Table 7: Multivariate analysis of FAPI PET TBRmax
# Cox proportional-hazards and linear regression results testing TBRmax as a predictor
# of aortic stenosis severity (generated in the statistics pipeline).
# Source file is copied directly as the formatting is produced by the statistics notebook.

fapi_mva_wb <- loadWorkbook(
  "../../05_statistics/03_output/02_FAPI_statistics/Supplemental_Table_7_TBRmax_MVA.xlsx"
)

saveWorkbook(fapi_mva_wb, file = paste0(supplemental_table_output_path, "Supplemental_Table_07.xlsx"), overwrite = TRUE)
