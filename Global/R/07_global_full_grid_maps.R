#!/usr/bin/env Rscript

# =============================================================================
# 07 - GLOBAL FEISTY FOUR-PANEL MAPS
# =============================================================================
# Purpose: use the completed one-degree summary from script 06 to create the
# four comparable map figures requested for the global analysis:
#   1. SSB at the beginning of the projection
#   2. Yield at the beginning of the projection
#   3. Absolute change in SSB
#   4. Absolute change in yield
#
# All panels within a figure use the same colour scale. Values beyond the
# pooled 98th percentile are capped for display only; the underlying results
# are unchanged.
#
# Input:  outputs/tables/global_full_grid_period_summary.csv
# Output: four PNG figures and a CSV recording each display limit
# =============================================================================

# ---- 1. Packages, paths and input checks ------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

input_file <- "outputs/tables/global_full_grid_period_summary.csv"
figures_dir <- "outputs/figures"

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Read and validate the full-grid period summary ----------------------

period_results <- fread(input_file)

required_columns <- c(
  "cell_id", "scenario", "lon", "lat", "metric", "group",
  "2015-2019", "2083-2087", "absolute_change"
)
missing_columns <- setdiff(required_columns, names(period_results))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

map_results <- period_results[scenario == "F040"]

expected_groups <- c("total", "smallPel", "demersals", "largePel")
expected_metrics <- c("SSB", "Yield")

if (!setequal(unique(map_results$group), expected_groups)) {
  stop("Unexpected functional groups in the F040 results")
}
if (!setequal(unique(map_results$metric), expected_metrics)) {
  stop("Unexpected metrics in the F040 results")
}

cell_count <- uniqueN(map_results$cell_id)
expected_rows <- cell_count * length(expected_groups) * length(expected_metrics)
if (nrow(map_results) != expected_rows) {
  stop(
    "Incomplete F040 map results: found ", nrow(map_results),
    " rows, expected ", expected_rows
  )
}

map_results[, panel := factor(
  group,
  levels = expected_groups,
  labels = c("Total", "Small pelagics", "Demersals", "Large pelagics")
)]

world_map <- map_data("world")

# ---- 3. Reusable function for one four-panel map ----------------------------

make_global_map <- function(
    data,
    value_column,
    title,
    subtitle,
    legend_title,
    output_file,
    change = FALSE
) {
  values <- data[[value_column]]
  values <- values[is.finite(values)]

  if (change) {
    scale_limit <- as.numeric(quantile(abs(values), 0.98, na.rm = TRUE))
  } else {
    scale_limit <- as.numeric(quantile(values, 0.98, na.rm = TRUE))
  }

  if (!is.finite(scale_limit) || scale_limit <= 0) {
    scale_limit <- 1
  }

  if (change) {
    fill_scale <- scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-scale_limit, scale_limit),
      oob = squish,
      labels = label_number(accuracy = 0.1),
      name = legend_title
    )
  } else {
    fill_scale <- scale_fill_viridis_c(
      option = "C",
      begin = 0.03,
      end = 0.95,
      limits = c(0, scale_limit),
      oob = squish,
      labels = label_number(accuracy = 0.1),
      name = legend_title
    )
  }

  plot_object <- ggplot() +
    geom_raster(
      data = data,
      aes(x = lon, y = lat, fill = .data[[value_column]]),
      interpolate = FALSE
    ) +
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = "grey88",
      colour = "grey55",
      linewidth = 0.08
    ) +
    facet_wrap(~panel, ncol = 2) +
    fill_scale +
    coord_quickmap(
      xlim = c(-180, 180),
      ylim = c(-80, 90),
      expand = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_bw(base_family = "Arial", base_size = 11) +
    theme(
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey75"),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 16, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 10.5, colour = "grey30", margin = margin(b = 10)),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8, colour = "grey25"),
      legend.position = "right",
      legend.title = element_text(size = 9.5),
      legend.text = element_text(size = 8.5),
      legend.key.height = grid::unit(1.4, "cm"),
      plot.margin = margin(12, 16, 10, 12)
    )

  ggsave(
    filename = file.path(figures_dir, output_file),
    plot = plot_object,
    width = 13.2,
    height = 8.6,
    dpi = 300,
    bg = "white"
  )

  invisible(list(plot = plot_object, scale_limit = scale_limit))
}

# ---- 4. Create the four requested figures -----------------------------------

ssb_data <- map_results[metric == "SSB"]
yield_data <- map_results[metric == "Yield"]

figure_info <- list(
  ssb_start = make_global_map(
    ssb_data,
    "2015-2019",
    "Global FEISTY: spawning-stock biomass at the beginning of the projection",
    "Mean 2015-2019 under F = 0.40 yr-1; common colour scale capped at the pooled 98th percentile",
    "SSB\n(g wet weight m-2)",
    "global_F040_reference_SSB_four_panel.png"
  ),
  yield_start = make_global_map(
    yield_data,
    "2015-2019",
    "Global FEISTY: yield at the beginning of the projection",
    "Mean 2015-2019 under F = 0.40 yr-1; common colour scale capped at the pooled 98th percentile",
    "Yield\n(g wet weight m-2 yr-1)",
    "global_F040_reference_yield_four_panel.png"
  ),
  ssb_change = make_global_map(
    ssb_data,
    "absolute_change",
    "Global FEISTY: absolute change in spawning-stock biomass",
    "Mean 2083-2087 minus mean 2015-2019 under F = 0.40 yr-1; common symmetric scale capped at the pooled 98th percentile",
    "Change in SSB\n(g wet weight m-2)",
    "global_F040_absolute_SSB_change_four_panel.png",
    change = TRUE
  ),
  yield_change = make_global_map(
    yield_data,
    "absolute_change",
    "Global FEISTY: absolute change in yield",
    "Mean 2083-2087 minus mean 2015-2019 under F = 0.40 yr-1; common symmetric scale capped at the pooled 98th percentile",
    "Change in yield\n(g wet weight m-2 yr-1)",
    "global_F040_absolute_yield_change_four_panel.png",
    change = TRUE
  )
)

# ---- 5. Save the colour limits used in every figure -------------------------

scale_summary <- data.table(
  figure = names(figure_info),
  display_limit = vapply(figure_info, function(x) x$scale_limit, numeric(1)),
  cap_percentile = 98,
  scenario = "F040",
  start_period = "2015-2019",
  end_period = "2083-2087",
  cells = cell_count
)

fwrite(scale_summary, file.path(figures_dir, "global_F040_map_display_scales.csv"))

cat("Global map figures completed\n")
cat("Cells represented:", format(cell_count, big.mark = ","), "\n")
cat("Scenario: F = 0.40 yr-1\n")
print(scale_summary)
