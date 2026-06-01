# Supplemental Figure 8a: Flow cytometry gating strategy panels (FSC/SSC, singlet, Calcein Violet vs DRAQ7) from Specimen_001_LCC_04.fcs
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

# Load data
ff <- read.FCS("../01_data/Specimen_001_LCC_04.fcs", transformation = FALSE)
expr <- as.data.frame(exprs(ff))
xch <- "FSC-A"
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
  scale_x_continuous(
    name = xch,
    breaks = c(0, 5e4, 1e5, 1.5e5, 2e5, 2.5e5),
    labels = label_number(scale = 1e-3, suffix = "K"),
    expand = c(0, 0), limits = c(0, 270000)
  ) +
  scale_y_log10(
    name = ych,
    breaks = x_breaks,
    labels = trans_format("log10", math_format(10^.x)),
    minor_breaks = NULL,
    expand = c(0, 0)
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
    20000,   # left-bottom
    25000,   # left-top (shorter than before)
    65000,   # top-left corner
    255000,  # top-right
    255000,  # right-mid
    150000,  # lower-right slope
    65000,   # lower-mid
    20000    # close
  ),
  y = c(
    20,    # ~6e3
    3000,   # ~9e4  (top is now 1e5-ish, not 1e6+)
    130000,  # ~1.3e5
    170000,  # ~1.7e5
    3000,   # ~3.5e4
    200,   # ~1.5e4
    100,    # ~8.5e3
    20
  )
)

plt <- plt + geom_path(
  data = gate, aes(x = x, y = y),
  inherit.aes = FALSE,
  colour = "black", linewidth = 0.2
)
ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_8a-1.svg"),
  plot = plt,
  width = 28, height = 28, units = "mm")

# Run point-in-polygon
pip <- sp::point.in.polygon(
  point.x = expr[[xch]],
  point.y = expr[[ych]],
  pol.x   = gate$x,
  pol.y   = gate$y
)

# pip returns:
# 0 = outside
# 1 = inside
# 2 = boundary
# 3 = vertex

expr$gate1_poly <- pip > 0   # include boundary + vertex

# Subset
expr_g1 <- expr[expr$gate1_poly, , drop = FALSE]

xch <- "FSC-A"
ych <- "FSC-H"

tmp <- ggplot_build(ggplot(expr_g1, aes(x = .data[[xch]], y = .data[[ych]])) + stat_bin2d(bins = 400))$data[[1]]$count
lim_hi <- quantile(log1p(tmp), 0.995, na.rm = TRUE) # tune 0.99-0.999

plt <- ggplot(expr_g1, aes(x = .data[[xch]], y = .data[[ych]])) +
  ggrastr::rasterise(
    stat_bin2d(aes(fill = after_stat(log1p(count))), bins = 350),
    dpi = 600
  ) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, lim_hi), oob = squish) +
  scale_x_continuous(
    name = xch,
    breaks = c(0, 5e4, 1e5, 1.5e5, 2e5, 2.5e5),
    labels = label_number(scale = 1e-3, suffix = "K"),
    expand = c(0, 0), limits = c(0, 270000)
  ) +
  scale_y_continuous(
    name = ych,
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
    10000,   # left vertical bend
    215000,  # top-right
    225000,  # right drop
    120000,  # lower-right slope
    10000    # close polygon
  ),
  y = c(
    5000,    # bottom-left
    95000,   # upper-left
    195000,  # top-right
    120000,  # mid-right
    60000,   # lower slope
    5000
  )
)

plt <- plt + geom_path(
  data = gate, aes(x = x, y = y),
  inherit.aes = FALSE,
  colour = "black", linewidth = 0.2
)
ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_8a-2.svg"),
  plot = plt,
  width = 28, height = 28, units = "mm")

# Run point-in-polygon
pip <- sp::point.in.polygon(
  point.x = expr_g1[[xch]],
  point.y = expr_g1[[ych]],
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

xch <- "SSC-A"
ych <- "SSC-H"

tmp <- ggplot_build(ggplot(expr_g2, aes(x = .data[[xch]], y = .data[[ych]])) + stat_bin2d(bins = 400))$data[[1]]$count
lim_hi <- quantile(log1p(tmp), 0.995, na.rm = TRUE) # tune 0.99-0.999

plt <- ggplot(expr_g2, aes(x = .data[[xch]], y = .data[[ych]])) +
  ggrastr::rasterise(
    stat_bin2d(aes(fill = after_stat(log1p(count))), bins = 350),
    dpi = 600
  ) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, lim_hi), oob = squish) +
  scale_x_continuous(
    name = xch,
    breaks = c(0, 5e4, 1e5, 1.5e5, 2e5, 2.5e5),
    labels = label_number(scale = 1e-3, suffix = "K"),
    expand = c(0, 0), limits = c(0, 270000)
  ) +
  scale_y_continuous(
    name = ych,
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
    10000,       # bottom-left
    0,
    0,       # left-top
    215000,  # top-right
    225000,  # right-mid (drop)
    130000,  # lower-right (slope back)
    10000        # close
  ),
  y = c(
    0,   # bottom-left
    0,
    50000,   # left-top  (set > bottom-left so it's a real edge)
    195000,  # top-right
    140000,  # right-mid
    90000,   # lower-right
    0    # close
  )
)

plt <- plt + geom_path(
  data = gate, aes(x = x, y = y),
  inherit.aes = FALSE,
  colour = "black", linewidth = 0.2
)
ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_8a-3.svg"),
  plot = plt,
  width = 28, height = 28, units = "mm")

# Run point-in-polygon
pip <- sp::point.in.polygon(
  point.x = expr_g2[[xch]],
  point.y = expr_g2[[ych]],
  pol.x   = gate$x,
  pol.y   = gate$y
)

# pip returns:
# 0 = outside
# 1 = inside
# 2 = boundary
# 3 = vertex

expr_g2$gate3_poly <- pip > 0   # include boundary + vertex

# Subset
expr_g3 <- expr_g2[expr_g2$gate3_poly, , drop = FALSE]

xch <- "V 450/50-A"
ych <- "R 780/60-A"

# choose *few* labeled decades (adjust to your data range)
x_breaks <- 10^(1:5)

tmp <- ggplot_build(ggplot(expr_g3, aes(x = .data[[xch]], y = .data[[ych]])) + stat_bin2d(bins = 400))$data[[1]]$count
lim_hi <- quantile(log1p(tmp), 0.995, na.rm = TRUE) # tune 0.99-0.999
logicle_labels <- function(x) {
  sapply(x, function(val) {
    if (val == 0) {
      return(expression(0))
    } else if (val < 0) {
      bquote(-10^.(log10(abs(val))))
    } else {
      bquote(10^.(log10(val)))
    }
  })
}
plt <- ggplot(expr_g3, aes(x = .data[[xch]], y = .data[[ych]])) +
  ggrastr::rasterise(
    stat_bin2d(aes(fill = after_stat(log1p(count))), bins = 350),
    dpi = 600
  ) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, lim_hi), oob = squish) +
  scale_x_logicle(
    name = "Calcein Violet",
    breaks = x_breaks,
    labels = logicle_labels,
    minor_breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_y_logicle(
    name = "DRAQ7",
    breaks = x_breaks,
    labels = logicle_labels,
    minor_breaks = NULL,
    expand = c(0, 0)
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
    200,   # left-bottom
    200,   # left-top (shorter than before)
    200000,   # top-left corner
    200000,  # top-right
    200    # close
  ),
  y = c(
    -300,    # ~6e3
    6000,   # ~9e4  (top is now 1e5-ish, not 1e6+)
    6000,  # ~1.3e5
    -300,  # ~1.7e5
    -300
  )
)

plt <- plt + geom_path(
  data = gate, aes(x = x, y = y),
  inherit.aes = FALSE,
  colour = "black", linewidth = 0.2
)
ggsave(paste0(supplemental_figure_output_path, "Supplemental_Figure_8a-4.svg"),
  plot = plt,
  width = 28, height = 28, units = "mm")

source_panel1 <- make_bin2d_source_data(
  df = expr,
  xch = "FSC-A",
  ych = "SSC-A",
  bins = 350,
  xlim = c(0, 270000)
)

source_panel2 <- make_bin2d_source_data(
  df = expr_g1,
  xch = "FSC-A",
  ych = "FSC-H",
  bins = 350,
  xlim = c(0, 270000),
  ylim = c(0, 270000)
)

source_panel3 <- make_bin2d_source_data(
  df = expr_g2,
  xch = "SSC-A",
  ych = "SSC-H",
  bins = 350,
  xlim = c(0, 270000),
  ylim = c(0, 270000)
)
source_panel4 <- make_bin2d_source_data(
  df = expr_g3,
  xch = "V 450/50-A",
  ych = "R 780/60-A",
  bins = 350
)

source_data_tables <- list(
  "Panel 1 FSC-A vs SSC-A" = source_panel1,
  "Panel 2 FSC-A vs FSC-H" = source_panel2,
  "Panel 3 SSC-A vs SSC-H" = source_panel3,
  "Panel 4 Calcein Violet vs DRAQ7" = source_panel4
)

wb <- createWorkbook()

for (sheet_name in names(source_data_tables)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, source_data_tables[[sheet_name]])
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(
    wb,
    sheet_name,
    cols = 1:ncol(source_data_tables[[sheet_name]]),
    widths = "auto"
  )
}

saveWorkbook(
  wb,
  file.path(source_data_output_path, "Supplemental_Figure_8a.xlsx"),
  overwrite = TRUE
)
