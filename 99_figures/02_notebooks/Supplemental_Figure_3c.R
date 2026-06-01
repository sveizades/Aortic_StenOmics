# Supplemental Figure 3c: DAPI intensity violin plots per donor for UBC snRNA-seq nuclei quality control.
library(flowCore)
library(ggcyto)
library(scales)
library(ggrastr)
library(sp)
library(scales)

source("../../utils.R")

ff_C4J <- read.FCS("../01_data/Specimen_001_C4J_Dapi.fcs", transformation = FALSE)
expr_C4J <- as.data.frame(exprs(ff_C4J))
ff_MK0 <- read.FCS("../01_data/Specimen_001_MK0 Sort.fcs", transformation = FALSE)
expr_MK0 <- as.data.frame(exprs(ff_MK0))
ff_HFU <- read.FCS("../01_data/Specimen_001_HFU Sort.fcs", transformation = FALSE)
expr_HFU <- as.data.frame(exprs(ff_HFU))

expr_C4J <- expr_C4J[expr_C4J$`V 450/50-A` > 20000 & expr_C4J$`V 450/50-A` < 80000, ]
expr_MK0 <- expr_MK0[expr_MK0$`V 450/50-A` > 20000 & expr_MK0$`V 450/50-A` < 80000, ]
expr_HFU <- expr_HFU[expr_HFU$`V 450/50-A` > 20000 & expr_HFU$`V 450/50-A` < 80000, ]

plot_df <- data.frame(sample = c(rep("C4J\n(excluded)", nrow(expr_C4J)), rep("MK0", nrow(expr_MK0)), rep("HFU", nrow(expr_HFU))), DAPI = c(expr_C4J$`V 450/50-A`, expr_MK0$`V 450/50-A`, expr_HFU$`V 450/50-A`))

ubc_donor_colours <- c(AKY = "#91cd2b", `C4J\n(excluded)` = "grey1", C5H = "#02dd83", CBN = "#b30090", DA5 = "#9ea800", GW1 = "#8984ff", HFU = "#e6a600", I4Q = "#0058b4", INV = "#f17f0b", JI1 = "#0087cc", KS4 = "#b53600", LTO = "#27dec3", M6Y = "#e71b95", MK0 = "#7d7600", O2N = "#ff6ed7", PGS = "#8a4a00", PJN = "#d8a4ff", PUM = "#734a0f", Q8H = "#ff529d", QZJ = "#8a7745", R9L = "#d50035", RW7 = "#75406b", SFA = "#fcb787", UKG = "#a31a1a", V5R = "#ff9ccd", XZD = "#913040", Y84 = "#e89da0", YDA = "#ff6f82")

plt <- ggplot(plot_df, aes(x = sample, y = DAPI, fill = sample)) +
  geom_violin(trim = FALSE, alpha = 0.7, colour = NA) +
  scale_fill_manual(values = ubc_donor_colours) +
  scale_y_log10(labels = label_number(scale = 1/1000, suffix = "k")) +
  labs(x = NULL, y = expression(DAPI~(log[10]))) +
  theme_classic(base_size = 5) +
  theme(
    axis.title.y = element_text(size = 5),
    axis.text.x  = element_text(size = 5),
    axis.text.y  = element_text(size = 5),
    axis.line  = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(0.5, "mm"),
    plot.margin = margin(1, 1, 1, 1, "mm"), legend.position = "none"
  )

write.csv(plot_df, file = paste0(source_data_output_path, "Supplemental_Figure_3c.csv"), row.names = TRUE)

ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_3c.svg"),
       plot = plt,
       width = 28, height = 28, units = "mm")
