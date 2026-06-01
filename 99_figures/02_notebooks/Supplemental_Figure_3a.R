# Supplemental Figure 3a: Flow cytometry gating strategy density plots for UBC snRNA-seq nuclei isolation.
library(flowCore)
library(ggcyto)
library(scales)
library(ggrastr)
library(sp)
library(openxlsx)

source("../../utils.R")

make_bin2d_source_data <- function(df, xch, ych, bins = 350, xlim = NULL, ylim = NULL) {

  p <- ggplot(df, aes(x = .data[[xch]], y = .data[[ych]])) +
    stat_bin2d(bins = bins)

  gb <- ggplot_build(p)$data[[1]]

  out <- gb[, c("xmin", "xmax", "ymin", "ymax", "count", "density", "ncount", "ndensity")]
  out$x_channel <- xch
  out$y_channel <- ych
  out$log1p_count <- log1p(out$count)

  if (!is.null(xlim)) {
    out <- out[out$xmin >= xlim[1] & out$xmax <= xlim[2], ]
  }
  if (!is.null(ylim)) {
    out <- out[out$ymin >= ylim[1] & out$ymax <= ylim[2], ]
  }

  out
}

# Read an FCS
ff <- read.FCS("../01_data/Specimen_001_UKG sort_007.fcs", transformation = FALSE)
expr <- as.data.frame(exprs(ff))
xch <- "V450/50-A"
ych <- "SSC-A"

# choose *few* labeled decades (adjust to your data range)
x_breaks <- 10^(1:5)

tmp <- ggplot_build(ggplot(expr, aes(x = .data[[xch]], y = .data[[ych]])) + stat_bin2d(bins = 400))$data[[1]]$count
lim_hi <- quantile(log1p(tmp), 0.995, na.rm = TRUE) # tune 0.99-0.999

plt <- ggplot(expr, aes(x = .data[[xch]], y = .data[[ych]])) +
  ggrastr::rasterise(
    stat_bin2d(aes(fill = after_stat(log1p(count))), bins = 350),
    dpi = 600
  ) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, lim_hi), oob = squish) +
  scale_x_log10(
    name = "DAPI",
    breaks = x_breaks,
    labels = trans_format("log10", math_format(10^.x)),
    minor_breaks = NULL,        # kill minor breaks (cleaner)
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "SSC-A",
    breaks = c(0, 5e4, 1e5, 1.5e5, 2e5, 2.5e5),
    labels = label_number(scale = 1e-3, suffix = "K"),
    expand = c(0, 0), limits = c(0, 270000)
  ) +
  theme_classic(base_size = 5) +
  theme(
    axis.title.x = element_text(size = 5),
    axis.title.y = element_text(size = 5),
    axis.text.x  = element_text(size = 3),
    axis.text.y  = element_text(size = 3),
    axis.line  = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(0.5, "mm"),
    plot.margin = margin(1, 1, 1, 1, "mm"), legend.position = "none"
  )

gate <- data.frame(
  x = c(45000, 45000, 170000, 170000, 45000),
  y = c(0, 250000, 250000, 0, 0)
)

plt <- plt + geom_path(
  data = gate, aes(x = x, y = y),
  inherit.aes = FALSE,
  colour = "black", linewidth = 0.2
)

ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_3a-1.svg"),
       plot = plt,
       width = 28, height = 28, units = "mm")

source_panel1 <- make_bin2d_source_data(
  df = expr,
  xch = "V450/50-A",
  ych = "SSC-A",
  bins = 350,
  ylim = c(0, 270000)
)

xmin <- 45000
xmax <- 170000
ymin <- 0
ymax <- 250000
expr$gate1_dapi <-
  expr[[xch]] >= xmin & expr[[xch]] <= xmax &
  expr[[ych]] >= ymin & expr[[ych]] <= ymax

expr_g1 <- expr[expr$gate1_dapi, , drop = FALSE]
ych <- "SSC-H"
xch <- "SSC-A"

tmp <- ggplot_build(ggplot(expr_g1, aes(x = .data[[xch]], y = .data[[ych]])) + stat_bin2d(bins = 400))$data[[1]]$count
lim_hi <- quantile(log1p(tmp), 0.995, na.rm = TRUE) # tune 0.99-0.999

plt <- ggplot(expr_g1, aes(x = .data[[xch]], y = .data[[ych]])) +
  ggrastr::rasterise(
    stat_bin2d(aes(fill = after_stat(log1p(count))), bins = 350),
    dpi = 600
  ) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, lim_hi), oob = squish) +
  scale_x_continuous(
    name = "SSC-A",
    breaks = c(0, 5e4, 1e5, 1.5e5, 2e5, 2.5e5),
    labels = label_number(scale = 1e-3, suffix = "K"),
    expand = c(0, 0), limits = c(0, 270000)
  ) +
  scale_y_continuous(
    name = "SSC-H",
    breaks = c(0, 5e4, 1e5, 1.5e5, 2e5, 2.5e5),
    labels = label_number(scale = 1e-3, suffix = "K"),
    expand = c(0, 0), limits = c(0, 270000)
  ) +
  theme_classic(base_size = 5) +
  theme(
    axis.title.x = element_text(size = 5),
    axis.title.y = element_text(size = 5),
    axis.text.x  = element_text(size = 3),
    axis.text.y  = element_text(size = 3),
    axis.line  = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(0.5, "mm"),
    plot.margin = margin(1, 1, 1, 1, "mm"), legend.position = "none"
  )

gate <- data.frame(
  x = c(220000, 150000, 0.1, 0.1, 10000, 220000),
  y = c(80000, 150000, 30000, 0.1, 0.1, 80000)
)

plt <- plt + geom_path(
  data = gate, aes(x = x, y = y),
  inherit.aes = FALSE,
  colour = "black", linewidth = 0.2
)

ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_3a-2.svg"),
  plot = plt,
  width = 28, height = 28, units = "mm")

source_panel2 <- make_bin2d_source_data(
  df = expr_g1,
  xch = "SSC-A",
  ych = "SSC-H",
  bins = 350,
  xlim = c(0, 270000),
  ylim = c(0, 270000)
)

# Your polygon (must be closed OR sp will treat it as closed automatically)
gate <- data.frame(
  x = c(220000, 150000, 0, 0, 10000, 220000),
  y = c(80000, 150000, 30000, 0, 0, 80000)
)

# Choose the axes used for this gate
x2 <- "SSC-A"
y2 <- "SSC-H"

# Run point-in-polygon
pip <- sp::point.in.polygon(
  point.x = expr_g1[[x2]],
  point.y = expr_g1[[y2]],
  pol.x   = gate$x,
  pol.y   = gate$y
)

# pip returns:
# 0 = outside
# 1 = inside
# 2 = boundary
# 3 = vertex

expr_g1$gate2_poly <- pip > 0   # include boundary + vertex

# Subset
expr_g2 <- expr_g1[expr_g1$gate2_poly, , drop = FALSE]

ych <- "SSC-A"
xch <- "FSC-A"
tmp <- ggplot_build(ggplot(expr_g2, aes(x = .data[[xch]], y = .data[[ych]])) + stat_bin2d(bins = 400))$data[[1]]$count
lim_hi <- quantile(log1p(tmp), 0.995, na.rm = TRUE) # tune 0.99-0.999

plt <- ggplot(expr_g2, aes(x = .data[[xch]], y = .data[[ych]])) +
  ggrastr::rasterise(
    stat_bin2d(aes(fill = after_stat(log1p(count))), bins = 350),
    dpi = 600
  ) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, lim_hi), oob = squish) +
  scale_x_continuous(
    name = "FSC-A",
    breaks = c(0, 5e4, 1e5, 1.5e5, 2e5, 2.5e5),
    labels = label_number(scale = 1e-3, suffix = "K"),
    expand = c(0, 0), limits = c(0, 270000)
  ) +
  scale_y_continuous(
    name = "SSC-A",
    breaks = c(0, 5e4, 1e5, 1.5e5, 2e5, 2.5e5),
    labels = label_number(scale = 1e-3, suffix = "K"),
    expand = c(0, 0), limits = c(0, 270000)
  ) +
  theme_classic(base_size = 5) +
  theme(
    axis.title.x = element_text(size = 5),
    axis.title.y = element_text(size = 5),
    axis.text.x  = element_text(size = 3),
    axis.text.y  = element_text(size = 3),
    axis.line  = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(0.5, "mm"),
    plot.margin = margin(1, 1, 1, 1, "mm"), legend.position = "none"
  )

gate <- data.frame(
  x = c(
    10000,   # bottom-left
    10000,   # left mid
    45000,   # upper-left slope
    150000,  # top-right
    150000,  # right mid
    135000,  # lower-right slope
    50000,   # lower mid
    10000    # close
  ),
  y = c(
    5000,   # bottom-left
    60000,   # left mid
    220000,  # top-left
    220000,  # top-right
    70000,   # right mid
    2000,   # lower-right
    5000,   # lower mid
    5000    # close
  )
)

plt <- plt + geom_path(
  data = gate, aes(x = x, y = y),
  inherit.aes = FALSE,
  colour = "black", linewidth = 0.2
)

ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_3a-3.svg"),
  plot = plt,
  width = 28, height = 28, units = "mm")

source_panel3 <- make_bin2d_source_data(
  df = expr_g2,
  xch = "FSC-A",
  ych = "SSC-A",
  bins = 350,
  xlim = c(0, 270000),
  ylim = c(0, 270000)
)

source_data_tables <- list(
  "Panel 1 DAPI vs SSC-A" = source_panel1,
  "Panel 2 SSC-A vs SSC-H" = source_panel2,
  "Panel 3 FSC-A vs SSC-A" = source_panel3
)

wb <- createWorkbook()

for (sheet_name in names(source_data_tables)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, source_data_tables[[sheet_name]])
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(source_data_tables[[sheet_name]]), widths = "auto")
}

saveWorkbook(
  wb,
  file.path(source_data_output_path, "Supplemental_Figure_3a.xlsx"),
  overwrite = TRUE
)
