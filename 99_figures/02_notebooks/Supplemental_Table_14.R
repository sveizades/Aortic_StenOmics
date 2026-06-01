# Supplemental Table 11: Statistical tests reported in the manuscript

library(openxlsx)

source("../../utils.R")

safe_read_csv <- function(path, ...) {
  if (!file.exists(path)) {
    return(data.frame(note = paste("Missing source file:", path)))
  }
  read.csv(path, check.names = FALSE, ...)
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

write_sheet <- function(wb, sheet, data) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data)
  freezePane(wb, sheet, firstActiveRow = 2)
  setColWidths(wb, sheet, cols = seq_along(data), widths = "auto")
}

catalogue <- data.frame(
  figure_or_table = c(
    "Figure 1f; Supplemental Table 2",
    "Supplemental Table 3",
    "Figure 2d; Supplemental Figure 4a-b",
    "Figure 2e; Supplemental Table 6",
    "Figure 2l",
    "Figure 3c",
    "Figure 3d-g",
    "Supplemental Table 7",
    "Figure 4g",
    "Figure 4h",
    "Figure 4k-l",
    "Figure 4m",
    "Figure 5b-d",
    "Figure 5f",
    "Supplemental Figure 3j",
    "Supplemental Figure 10f",
    "Supplemental Figure 10g; Supplemental Table 9",
    "Supplemental Figure 10h",
    "Supplemental Figure 10l",
    "Supplemental Tables 4, 5, 8 and 10"
  ),
  analysis = c(
    "Bulk RNA-seq differential expression",
    "Bulk RNA-seq GO biological process enrichment",
    "snRNA-seq VIC subtype proportions by disease severity",
    "snRNA-seq VIC marker GO biological process enrichment",
    "Multiplex immunofluorescence regional quantification",
    "[68Ga]FAPI-46 TBRmax trend across disease severity",
    "[68Ga]FAPI-46 TBRmax correlations with imaging measures",
    "Multivariable predictors of valvular [68Ga]FAPI-46 TBRmax",
    "Spatial enrichment or depletion relative to calcific nodules",
    "Spatial niche cell-type enrichment",
    "CellPhoneDB ligand-receptor interactions",
    "drug2cell candidate perturbagen enrichment",
    "Follow-up imaging and haemodynamic correlations",
    "Mediation of TBRmax association with haemodynamic progression by CT calcium progression",
    "Aortic stenosis module score across snRNA-seq broad cell types",
    "Cross-platform VIC subtype expression correlations",
    "Xenium spatial transcriptomics marker genes",
    "Distance to calcific nodule edge by Xenium VIC subtype",
    "Osteogenic module score by spatial niche in Xenium VICs",
    "Single-cell, single-nucleus and niche marker genes"
  ),
  statistical_test = c(
    "DESeq2 Wald test",
    "clusterProfiler enrichGO over-representation test",
    "limma linear model with robust empirical Bayes moderation on logit-transformed proportions",
    "clusterProfiler enrichGO over-representation test",
    "Linear mixed-effects model with emmeans pairwise contrasts",
    "Linear model trend test",
    "Pearson product-moment correlation",
    "Multivariable linear regression with HC3 robust standard errors",
    "Binomial generalised linear mixed-effects model",
    "scNiche enrichment analysis",
    "CellPhoneDB permutation-based statistical test",
    "drug2cell rank_genes_groups enrichment test",
    "Pearson product-moment correlation",
    "Causal mediation analysis with nonparametric bootstrap",
    "Pairwise Wilcoxon rank-sum test",
    "Pearson product-moment correlation",
    "Seurat FindAllMarkers Wilcoxon rank-sum test",
    "Kruskal-Wallis test and pairwise Wilcoxon rank-sum test",
    "Pairwise Wilcoxon rank-sum test",
    "Seurat FindAllMarkers Wilcoxon rank-sum test"
  ),
  model_or_comparison = c(
    "Clinical classification and AS versus no-AS contrasts",
    "GO BP enrichment among bulk RNA-seq DE genes",
    "Progressive linear disease-severity trend and severe-specific contrast",
    "GO BP enrichment among VIC subtype marker genes",
    "logit(percent stained area) ~ region + (1 | slide_id)",
    "TBRmax ~ ordinal disease severity score",
    "TBRmax versus CT calcium score, peak aortic jet velocity, non-calcific leaflet volume and calcific leaflet volume",
    "log(TBRmax) ~ age + sex + hypertension + diabetes + calcific volume + non-calcific volume",
    "Cell-type membership ~ z-scaled distance to nodule edge + (1 | donor_id)",
    "Cell-type over-representation within inferred spatial niches",
    "Macrophage and VIC-FAP-osteogenic ligand-receptor pairs",
    "Drug target gene enrichment across scRNA-seq VIC subtypes",
    "Baseline TBRmax, annualised CT calcium score change and annualised peak aortic jet velocity change",
    "Baseline TBRmax -> annualised CT calcium score change -> annualised peak aortic jet velocity change",
    "Aortic stenosis module score, comparisons against VICs shown in manuscript",
    "Average subtype expression correlation between Xenium, snRNA-seq and scRNA-seq",
    "Cluster marker genes for Xenium cell states and niches",
    "Distance to nodule edge across VIC subtypes; osteogenic versus contractile comparison shown",
    "Osteogenic module score in Peri-nodular niche versus other niches",
    "Marker genes for reported cell state annotations"
  ),
  multiple_testing = c(
    "Benjamini-Hochberg FDR",
    "Benjamini-Hochberg FDR",
    "Benjamini-Hochberg FDR",
    "Benjamini-Hochberg FDR",
    "Holm adjustment for pairwise contrasts",
    "None",
    "None",
    "None",
    "Benjamini-Hochberg FDR",
    "scNiche-reported enrichment p-values",
    "CellPhoneDB p < 0.05 threshold for displayed interactions",
    "Adjusted p-values from drug2cell output",
    "None",
    "Bootstrap confidence intervals and p-values",
    "Benjamini-Hochberg FDR",
    "None",
    "Seurat-adjusted p-values",
    "Benjamini-Hochberg FDR for pairwise comparison",
    "Benjamini-Hochberg FDR",
    "Seurat-adjusted p-values"
  ),
  source_code = c(
    "01_bulk_seq/02_notebooks/01_bulk_seq_data_analysis.Rmd",
    "01_bulk_seq/02_notebooks/01_bulk_seq_data_analysis.Rmd",
    "02_UBC_snRNA_seq/02_notebooks/04_UBC_snRNAseq_vic_annotation.Rmd",
    "02_UBC_snRNA_seq/02_notebooks/04_UBC_snRNAseq_vic_annotation.Rmd",
    "05_statistics/02_notebooks/01_vic_mIF_statistics.Rmd",
    "05_statistics/02_notebooks/02_FAPI_statistics.Rmd",
    "05_statistics/02_notebooks/02_FAPI_statistics.Rmd",
    "05_statistics/02_notebooks/02_FAPI_statistics.Rmd",
    "99_figures/02_notebooks/Figure_4g.R",
    "04_SALTIRE3_Xenium/02_notebooks/01_S3_Xenium_preprocessing.Rmd",
    "99_figures/02_notebooks/Figure_4k.R; 99_figures/02_notebooks/Figure_4l.R",
    "99_figures/02_notebooks/Figure_4m.R",
    "05_statistics/02_notebooks/03_FAPI_follow_up.Rmd",
    "05_statistics/02_notebooks/03_FAPI_follow_up.Rmd",
    "02_UBC_snRNA_seq/02_notebooks/11_UBC_snRNAseq_annotation_combining.Rmd",
    "99_figures/02_notebooks/Supplemental_Figure_10f-g.R",
    "99_figures/02_notebooks/Supplemental_Table_09.R",
    "99_figures/02_notebooks/Supplemental_Figure_10h.R",
    "99_figures/02_notebooks/Supplemental_Figure_10l.R",
    "Supplemental_Table_04.R; Supplemental_Table_05.R; Supplemental_Table_08.R; Supplemental_Table_10.R"
  ),
  result_location = c(
    "99_figures/03_output/03_supplemental_tables/Supplemental_Table_02.xlsx",
    "99_figures/03_output/03_supplemental_tables/Supplemental_Table_03.xlsx",
    "Sheets in this workbook: F2d_SF4_VIC_prop",
    "99_figures/03_output/03_supplemental_tables/Supplemental_Table_06.xlsx",
    "Sheet in this workbook: F2l_mIF",
    "Sheet in this workbook: F3_FAPI_stats",
    "Sheet in this workbook: F3_FAPI_stats",
    "Sheet in this workbook: ST7_MVA",
    "Sheet in this workbook: F4g_nodule_GLMM",
    "99_figures/03_output/04_source_data/Figure_4h.csv",
    "Sheet in this workbook: F4k_l_CPDB",
    "Sheet in this workbook: F4m_drug2cell",
    "Sheet in this workbook: F5_follow_up",
    "Sheet in this workbook: F5_follow_up",
    "Sheet in this workbook: SF3j_AS_score",
    "99_figures/03_output/04_source_data/Supplemental_Figure_10f.xlsx",
    "99_figures/03_output/03_supplemental_tables/Supplemental_Table_09.xlsx",
    "Sheet in this workbook: SF10h_distance",
    "Sheet in this workbook: SF10l_osteo_score",
    "99_figures/03_output/03_supplemental_tables/"
  )
)

vic_progressive <- safe_read_csv("../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/VIC_proportions_statistics.csv")
vic_progressive$contrast <- "Progressive disease-severity trend"
vic_severe <- safe_read_csv("../../02_UBC_snRNA_seq/03_output/04_UBC_snRNAseq_vic_annotation/VIC_proportions_severe_statistics.csv")
vic_severe$contrast <- "Severe-specific contrast"
vic_prop_stats <- rbind(vic_progressive, vic_severe)
vic_prop_stats <- vic_prop_stats[, c("contrast", setdiff(names(vic_prop_stats), "contrast"))]

mif_stats <- safe_read_csv("../../05_statistics/03_output/01_vic_mIF_statistics/mIF_region_LMM_pairwise_with_means.csv")

fapi_trend <- safe_read_csv("../../05_statistics/03_output/02_FAPI_statistics/Figure_3c_ASdiagnosis_FAPI_trend_stats.csv")
fapi_trend <- subset(fapi_trend, term == "severity_score")
fapi_trend_out <- data.frame(
  figure_panel = "Figure 3c",
  analysis = "TBRmax trend across disease severity",
  test = "Linear model trend test",
  effect = fapi_trend$estimate,
  statistic = fapi_trend$statistic,
  ci_lower = fapi_trend$conf.low,
  ci_upper = fapi_trend$conf.high,
  p_value = fapi_trend$p.value,
  p_label = format_p(fapi_trend$p.value),
  n = NA
)

correlation_files <- data.frame(
  figure_panel = c("Figure 3d", "Figure 3e", "Figure 3f", "Figure 3g"),
  analysis = c(
    "TBRmax versus log1p CT calcium score",
    "Peak aortic jet velocity versus TBRmax",
    "Non-calcific leaflet volume versus TBRmax",
    "log1p calcific leaflet volume versus TBRmax"
  ),
  path = c(
    "../../05_statistics/03_output/02_FAPI_statistics/Figure_3d_CTcalcium_TBRmax_correlation.csv",
    "../../05_statistics/03_output/02_FAPI_statistics/Figure_3e_Vmax_TBRmax_correlation.csv",
    "../../05_statistics/03_output/02_FAPI_statistics/Figure_3f_NonCalcificThickening_TBRmax_correlation.csv",
    "../../05_statistics/03_output/02_FAPI_statistics/Figure_3g_CalcificThickening_TBRmax_correlation.csv"
  )
)
fapi_correlations <- do.call(rbind, lapply(seq_len(nrow(correlation_files)), function(i) {
  dat <- safe_read_csv(correlation_files$path[i])
  data.frame(
    figure_panel = correlation_files$figure_panel[i],
    analysis = correlation_files$analysis[i],
    test = dat$Method,
    effect = dat$Correlation,
    statistic = NA,
    ci_lower = dat$CI_lower,
    ci_upper = dat$CI_upper,
    p_value = dat$p_value,
    p_label = format_p(dat$p_value),
    n = dat$n
  )
}))
fapi_stats <- rbind(fapi_trend_out, fapi_correlations)

mva <- read.xlsx("../../05_statistics/03_output/02_FAPI_statistics/Supplemental_Table_7_TBRmax_MVA.xlsx", sheet = 1)

follow_up_files <- data.frame(
  figure_panel = c("Figure 5b", "Figure 5c", "Figure 5d"),
  analysis = c(
    "Baseline TBRmax versus annualised CT calcium score change",
    "Annualised peak aortic jet velocity change versus annualised CT calcium score change",
    "Baseline TBRmax versus annualised peak aortic jet velocity change"
  ),
  path = c(
    "../../05_statistics/03_output/03_FAPI_follow_up/TBRmax_CTcalc_correlation.csv",
    "../../05_statistics/03_output/03_FAPI_follow_up/jet_CTcalc_correlation.csv",
    "../../05_statistics/03_output/03_FAPI_follow_up/TBRmax_jet_correlation.csv"
  )
)
follow_up_cor <- do.call(rbind, lapply(seq_len(nrow(follow_up_files)), function(i) {
  dat <- safe_read_csv(follow_up_files$path[i])
  data.frame(
    figure_panel = follow_up_files$figure_panel[i],
    analysis = follow_up_files$analysis[i],
    term = "Pearson correlation",
    estimate = dat$estimate,
    statistic = dat$statistic,
    p_value = dat$p.value,
    p_label = format_p(dat$p.value),
    conf_low = dat$conf.low,
    conf_high = dat$conf.high,
    stringsAsFactors = FALSE
  )
}))
mediation <- safe_read_csv("../../05_statistics/03_output/03_FAPI_follow_up/mediation.csv")
mediation$figure_panel <- "Figure 5f"
mediation$analysis <- "Causal mediation analysis"
follow_up_stats <- rbind(
  follow_up_cor,
  data.frame(
    figure_panel = mediation$figure_panel,
    analysis = mediation$analysis,
    term = mediation$term,
    estimate = mediation$estimate,
    statistic = NA,
    p_value = mediation$p.value,
    p_label = format_p(mediation$p.value),
    conf_low = NA,
    conf_high = NA,
    stringsAsFactors = FALSE
  )
)

fig4g <- safe_read_csv("../../99_figures/03_output/04_source_data/Figure_4g.csv")
sf3j <- safe_read_csv("../../02_UBC_snRNA_seq/03_output/11_UBC_snRNAseq_annotation_combining/Supplementary_Table_CellType_Signature_Pairwise.csv")

sf10h_kw <- safe_read_csv("../../99_figures/03_output/04_source_data/Supplemental_Figure_10h_distance_to_nodule_edge_VIC_kw_stats.csv")
sf10h_cells <- safe_read_csv("../../99_figures/03_output/04_source_data/Supplemental_Figure_10h.csv")
sf10h_pairwise <- data.frame(
  comparison = "VIC-FAP-osteogenic versus VIC-contractile",
  p_value = pairwise.wilcox.test(
    sf10h_cells$distance_to_nodule_edge,
    sf10h_cells$annotations_level2,
    p.adjust.method = "BH"
  )$p.value["VIC-FAP-osteogenic", "VIC-contractile"],
  stringsAsFactors = FALSE
)
sf10h_pairwise$p_label <- format_p(sf10h_pairwise$p_value)
sf10h <- rbind(
  data.frame(
    test = "Kruskal-Wallis",
    comparison = "All Xenium VIC subtypes",
    statistic = sf10h_kw$statistic,
    df = sf10h_kw$df,
    p_value = sf10h_kw$p,
    p_label = format_p(sf10h_kw$p),
    stringsAsFactors = FALSE
  ),
  data.frame(
    test = "Pairwise Wilcoxon rank-sum test with BH adjustment",
    comparison = sf10h_pairwise$comparison,
    statistic = NA,
    df = NA,
    p_value = sf10h_pairwise$p_value,
    p_label = sf10h_pairwise$p_label,
    stringsAsFactors = FALSE
  )
)

sf10l_cells <- safe_read_csv("../../99_figures/03_output/04_source_data/Supplemental_Figure_10l.csv")
sf10l_cells <- subset(sf10l_cells, !is.na(niche) & !is.na(osteogenic_score1))
sf10l_comparisons <- setdiff(unique(sf10l_cells$niche), "Peri-nodular")
sf10l <- do.call(rbind, lapply(sf10l_comparisons, function(niche_i) {
  p <- wilcox.test(
    sf10l_cells$osteogenic_score1[sf10l_cells$niche == "Peri-nodular"],
    sf10l_cells$osteogenic_score1[sf10l_cells$niche == niche_i]
  )$p.value
  data.frame(
    comparison = paste("Peri-nodular versus", niche_i),
    p_value = p,
    stringsAsFactors = FALSE
  )
}))
sf10l$p_adj <- p.adjust(sf10l$p_value, method = "BH")
sf10l$p_label <- format_p(sf10l$p_adj)

fig4k <- safe_read_csv("../../99_figures/03_output/04_source_data/Figure_4k.csv")
fig4k$figure_panel <- "Figure 4k"
fig4l <- safe_read_csv("../../99_figures/03_output/04_source_data/Figure_4l.csv")
fig4l$figure_panel <- "Figure 4l"
cpdb <- rbind(fig4k, fig4l)
cpdb <- cpdb[, c("figure_panel", "interaction_readable", "Var1", "Var2", "scaled_means", "pvals", "neglog10_p", "classification")]

drug2cell <- safe_read_csv("../../99_figures/03_output/04_source_data/Figure_4m.csv")
drug2cell <- drug2cell[, c("group", "drug_name", "names", "scores", "logfoldchanges", "pvals", "pvals_adj")]

wb <- createWorkbook()

write_sheet(wb, "Test_catalogue", catalogue)
write_sheet(wb, "F2d_SF4_VIC_prop", vic_prop_stats)
write_sheet(wb, "F2l_mIF", mif_stats)
write_sheet(wb, "F3_FAPI_stats", fapi_stats)
write_sheet(wb, "ST7_MVA", mva)
write_sheet(wb, "F4g_nodule_GLMM", fig4g)
write_sheet(wb, "F4k_l_CPDB", cpdb)
write_sheet(wb, "F4m_drug2cell", drug2cell)
write_sheet(wb, "F5_follow_up", follow_up_stats)
write_sheet(wb, "SF3j_AS_score", sf3j)
write_sheet(wb, "SF10h_distance", sf10h)
write_sheet(wb, "SF10l_osteo_score", sf10l)

saveWorkbook(
  wb,
  file = paste0(supplemental_table_output_path, "Supplemental_Table_14.xlsx"),
  overwrite = TRUE
)
