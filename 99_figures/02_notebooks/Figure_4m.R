# Figure 4m: Drug2cell candidate perturbagens across VIC sub-states

library(dplyr)
library(ggplot2)
library(dplyr)
library(forcats)
library(stringr)
library(cowplot)

source("../../utils.R")
drug_de <- read.csv("../../03_SALTIRE3_scRNA_seq/01_data/drug2cell_DE_results.csv")

drugs <- c(
  "CHEMBL314854|FINGOLIMOD",
  "CHEMBL2336071|SIPONIMOD",
  "CHEMBL3707247|OZANIMOD",
  "CHEMBL1096146|PONESIMOD",
  "CHEMBL1201575|EFALIZUMAB",
  "CHEMBL4297734|CRIZANLIZUMAB",
  "CHEMBL1200699|DOXYCYCLINE",
  "CHEMBL1237023|DENOSUMAB",
  "CHEMBL3137301|SACUBITRIL"
)
ct_labels <- c(
  "VIC-quiescent"     = "Quiescent",
  "VIC-spongiosa"  = "Spongiosa",
  "VIC-contractile"           = "Contractile",
  "VIC-neural-crest"  = "Neural-crest-like",
  "VIC-transitional"  = "Transitional",
  "VIC-IFN"           = "IFN-stimulated",
  "VIC-FAP-osteogenic"    = "FAP+ osteogenic"
)
drug_labels <- c(
  "CHEMBL314854|FINGOLIMOD"              = "Fingolimod",
  "CHEMBL2336071|SIPONIMOD"              = "Siponimod",
  "CHEMBL3707247|OZANIMOD"              = "Ozanimod",
  "CHEMBL1096146|PONESIMOD"              = "Ponesimod",
  "CHEMBL1201575|EFALIZUMAB"              = "Efalizumab",
  "CHEMBL4297734|CRIZANLIZUMAB"              = "Crizanlizumab",
  "CHEMBL1200699|DOXYCYCLINE" = "Doxycycline",
  "CHEMBL1237023|DENOSUMAB"           = "Denosumab",
  "CHEMBL3137301|SACUBITRIL"          = "Sacubitril"
)


plot_df <- drug_de %>%
  filter(names %in% drugs) %>%
  mutate(
    drug_name = recode(names, !!!drug_labels)   # <- use same labels everywhere
  ) %>%
  group_by(names) %>%
  mutate(
    score_scaled = as.numeric(scale(scores))    # <- avoid matrix column from scale()
  ) %>%
  ungroup()


drug_order <- plot_df %>%
  filter(group == "VIC-FAP-osteogenic") %>%
  arrange(desc(logfoldchanges)) %>%
  pull(drug_name)


plot_df <- plot_df %>%
  mutate(
    drug_name = factor(drug_name, levels = unique(drug_order)),
    group = factor(
      group,
      levels = c(
        "VIC-IFN","VIC-neural-crest","VIC-spongiosa",
        "VIC-contractile","VIC-FAP-osteogenic","VIC-transitional","VIC-quiescent"
      )
    ),
    neglog10_padj = -log10(ifelse(pvals_adj == 0, 1e-300, pvals_adj)),
    neglog10_padj_plot = pmin(neglog10_padj, 5),
    score_scaled = pmax(score_scaled, 0),
    score_scaled = pmin(score_scaled, 1.5)
  )

plot_df$score_scaled[plot_df$score_scaled < 0] <- 0
plot_df$score_scaled[plot_df$score_scaled > 1.5] <- 1.5

p <- ggplot(plot_df, aes(x = drug_name, y = group)) +
  geom_point(aes(size = neglog10_padj_plot, fill = score_scaled), shape = 21, stroke = 0) +
  scale_size(range = c(0, 2.5)) +   # <- control dot sizes here
  scale_fill_distiller(
    palette = "Greens",
    direction = 1
  ) +
  scale_y_discrete(labels = ct_labels) +
  coord_equal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=5, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.05),
        axis.text.y = element_text(size=5, colour = "black", margin = margin(r = 1)),
        axis.title=element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),plot.background =element_blank(),panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"),legend.position = "none")

ggsave(plot = p, filename = paste0(main_figure_output_path, "Figure_4m.svg"), width = 5, height = 3, units = "cm")

# Save scale bar
df_legend <- data.frame(x = 1, y = 1, fill = 1:5)
# Create a plot with just the color scale
legend_plot <- ggplot(df_legend, aes(x = x, y = y, fill = fill)) +
  geom_tile() +  # just to trigger the fill scale
  scale_fill_distiller(palette = "Greens",direction = 1)+
  theme_void() +  # remove background and axes
  theme(legend.position = "right")
# Save only the legend as a plot
ggsave(paste0(main_figure_output_path, "Figure_4m-legend.svg"), plot = legend_plot, width = 3, height = 5)


size_legend_plot <- ggplot(plot_df, aes(x = drug_name, y = group)) +
  geom_point(aes(size = neglog10_padj_plot), shape = 21, fill = "grey50", stroke = 0.3) +
  scale_size(
    range = c(0, 2.5),
    name  = expression(-log[10](p[adj])),
    breaks = c(1, 2, 3, 4, 5),
    labels = c("1", "2", "3", "4", "≥5")
  ) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 5),
    legend.text      = element_text(size = 5),
    legend.key.size  = unit(3, "mm")
  )

size_legend <- get_legend(size_legend_plot)

ggsave(
  paste0(main_figure_output_path, "Figure_4m-size-legend.svg"),
  plot   = size_legend,
  width  = 2,
  height = 2,
  units  = "cm"
)

write.csv(plot_df, file = paste0(source_data_output_path, "Figure_4m.csv"), row.names = TRUE)
