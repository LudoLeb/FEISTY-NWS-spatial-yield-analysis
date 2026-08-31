#!/usr/bin/env Rscript

# =============================================================================
# FEISTY GLOBAL: OBSERVED VERSUS SIMULATED YIELD BY LME
# =============================================================================
#
# Objective
# ---------
# Compare observed catch and simulated FEISTY yield in each Large Marine
# Ecosystem (LME). One point in the final scatter plot represents one LME.
#
# Comparison setup
# ----------------
# - Observed catch: mean 1995-2004
# - Simulated yield: mean 2015-2019 under the constant F = 0.40 scenario
# - X axis: observed catch
# - Y axis: simulated yield
# - Four panels: small pelagics, large pelagics, demersals and their total
# - Logarithmic axes and a 1:1 reference line
#
# Important limitation
# --------------------
# The observed and simulated periods do not overlap. This figure is therefore
# a broad spatial comparison among LMEs, not a strict same-year validation.
# The global simulation also uses the same F in every cell, whereas real
# fishing pressure varies among LMEs.
#
# Required inputs
# ---------------
# 1. outputs/tables/global_full_grid_yield_2015_2019_for_LME.csv
# 2. grid_LME_1deg.RData, containing an object named `degdata`
# 3. catch_Ftype_LME_year.RData, containing an object named `tot`
#
# The two RData paths can be supplied without editing this script:
#   FEISTY_LME_GRID=/path/grid_LME_1deg.RData
#   FEISTY_LME_CATCH=/path/catch_Ftype_LME_year.RData
#
# Outputs
# -------
# - outputs/tables/global_LME_observed_vs_simulated_yield.csv
# - outputs/tables/global_LME_observed_vs_simulated_yield_statistics.csv
# - outputs/figures/global_LME_observed_vs_simulated_yield.png
# =============================================================================


# =============================================================================
# 1. PACKAGES AND USER SETTINGS
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# These settings are kept together so the comparison can be adjusted easily.
observed_years <- 1995:2004
simulated_period <- "2015-2019"
simulated_scenario <- "F040"
minimum_display_yield <- 1 # tonnes yr-1; keeps logarithmic axes readable

simulation_file <- Sys.getenv(
  "FEISTY_LME_SIMULATED_YIELD",
  unset = "outputs/tables/global_full_grid_yield_2015_2019_for_LME.csv"
)
lme_grid_file <- Sys.getenv(
  "FEISTY_LME_GRID",
  unset = "data/raw/grid_LME_1deg.RData"
)
observed_catch_file <- Sys.getenv(
  "FEISTY_LME_CATCH",
  unset = "data/raw/catch_Ftype_LME_year.RData"
)

output_table <- "outputs/tables/global_LME_observed_vs_simulated_yield.csv"
output_statistics <- paste0(
  "outputs/tables/global_LME_observed_vs_simulated_yield_statistics.csv"
)
output_figure <- "outputs/figures/global_LME_observed_vs_simulated_yield.png"

required_files <- c(simulation_file, lme_grid_file, observed_catch_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing input file(s): ", paste(missing_files, collapse = ", "),
    ". See data/README.md for the expected local files."
  )
}

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 2. LOAD THE COMPLETED FEISTY FULL-GRID RESULTS
# =============================================================================

message("Reading simulated FEISTY yield: ", simulation_file)
simulation <- fread(simulation_file)

required_simulation_columns <- c("lon", "lat", "group", simulated_period)
missing_simulation_columns <- setdiff(
  required_simulation_columns,
  names(simulation)
)
if (length(missing_simulation_columns) > 0L) {
  stop(
    "The simulation table is missing: ",
    paste(missing_simulation_columns, collapse = ", ")
  )
}

groups <- c("smallPel", "largePel", "demersals", "total")
simulation <- simulation[group %in% groups, .(
  lon,
  lat,
  group,
  simulated_yield_g_m2_yr = get(simulated_period)
)]


# =============================================================================
# 3. LOAD THE 1-DEGREE LME GRID AND MATCH IT TO FEISTY CELLS
# =============================================================================

lme_environment <- new.env(parent = emptyenv())
load(lme_grid_file, envir = lme_environment)
if (!exists("degdata", envir = lme_environment, inherits = FALSE)) {
  stop("The LME grid file must contain an object named `degdata`.")
}

lme_grid <- as.data.table(lme_environment$degdata)
required_lme_columns <- c("latlon", "area_msq", "ECO_NAME", "LME_nb")
missing_lme_columns <- setdiff(required_lme_columns, names(lme_grid))
if (length(missing_lme_columns) > 0L) {
  stop("The LME grid is missing: ", paste(missing_lme_columns, collapse = ", "))
}

# `latlon` stores latitude and longitude in one text field. Split it into the
# same coordinates used by the FEISTY one-degree grid.
coordinates <- tstrsplit(trimws(as.character(lme_grid$latlon)), "[[:space:]]+")
lme_grid[, `:=`(
  lat = as.numeric(coordinates[[1]]),
  lon = as.numeric(coordinates[[2]])
)]

if (anyDuplicated(lme_grid[, .(lon, lat)])) {
  stop("The LME grid contains duplicated longitude-latitude coordinates.")
}

simulation_lme <- merge(
  simulation,
  lme_grid[, .(lon, lat, LME_nb, ECO_NAME, area_msq)],
  by = c("lon", "lat"),
  all = FALSE
)

if (nrow(simulation_lme) == 0L) {
  stop("No FEISTY cells could be matched to the LME grid.")
}

message(
  "Matched ", uniqueN(simulation_lme[, .(lon, lat)]),
  " FEISTY cells across ", uniqueN(simulation_lme$LME_nb), " LMEs."
)


# =============================================================================
# 4. AGGREGATE SIMULATED YIELD WITHIN EACH LME
# =============================================================================

# FEISTY yield is in g wet weight m-2 yr-1. Multiplying by cell area gives
# grams per year; dividing by 1e6 converts grams to tonnes.
simulated_lme <- simulation_lme[, .(
  simulated_tonnes_yr = sum(
    simulated_yield_g_m2_yr * area_msq / 1e6,
    na.rm = TRUE
  ),
  matched_cells = uniqueN(paste(lon, lat))
), by = .(LME_nb, ECO_NAME, group)]


# =============================================================================
# 5. PREPARE THE OBSERVED CATCH FOR THE SAME FEISTY GROUPS
# =============================================================================

catch_environment <- new.env(parent = emptyenv())
load(observed_catch_file, envir = catch_environment)
if (!exists("tot", envir = catch_environment, inherits = FALSE)) {
  stop("The observed catch file must contain an object named `tot`.")
}

observed <- as.data.table(catch_environment$tot)
required_observed_columns <- c("LME_nb", "Ftype", "Year", "Catch_tonnes")
missing_observed_columns <- setdiff(required_observed_columns, names(observed))
if (length(missing_observed_columns) > 0L) {
  stop(
    "The observed catch table is missing: ",
    paste(missing_observed_columns, collapse = ", ")
  )
}

# The observed functional types that directly correspond to FEISTY groups.
group_key <- data.table(
  Ftype = c("spel", "lpel", "dem"),
  group = c("smallPel", "largePel", "demersals")
)

valid_lmes <- sort(unique(lme_grid$LME_nb))
observed <- observed[
  LME_nb %in% valid_lmes &
    Year %in% observed_years &
    Ftype %in% group_key$Ftype
]

# Rows absent from the aggregated catch file represent zero recorded catch for
# that LME, group and year. Completing the grid ensures that each LME mean uses
# the same ten years.
observed_complete <- CJ(
  LME_nb = valid_lmes,
  Year = observed_years,
  Ftype = group_key$Ftype,
  unique = TRUE
)
observed_complete <- merge(
  observed_complete,
  observed[, .(LME_nb, Year, Ftype, Catch_tonnes)],
  by = c("LME_nb", "Year", "Ftype"),
  all.x = TRUE
)
observed_complete[is.na(Catch_tonnes), Catch_tonnes := 0]
observed_complete <- merge(observed_complete, group_key, by = "Ftype")

observed_groups <- observed_complete[, .(
  observed_tonnes_yr = mean(Catch_tonnes)
), by = .(LME_nb, group)]

# The observed total contains only the three groups represented by FEISTY.
# Other reported categories are deliberately excluded to compare like with like.
observed_total <- observed_complete[, .(
  annual_total = sum(Catch_tonnes)
), by = .(LME_nb, Year)][, .(
  observed_tonnes_yr = mean(annual_total)
), by = LME_nb]
observed_total[, group := "total"]

observed_lme <- rbindlist(list(observed_groups, observed_total), use.names = TRUE)


# =============================================================================
# 6. MERGE OBSERVED AND SIMULATED VALUES AND CALCULATE DIAGNOSTICS
# =============================================================================

comparison <- merge(
  simulated_lme,
  observed_lme,
  by = c("LME_nb", "group"),
  all = FALSE
)

comparison[, observed_period := paste0(min(observed_years), "-", max(observed_years))]
comparison[, simulated_period := simulated_period]
comparison[, simulated_scenario := simulated_scenario]
comparison[, ratio_simulated_observed := fifelse(
  observed_tonnes_yr > 0,
  simulated_tonnes_yr / observed_tonnes_yr,
  NA_real_
)]
comparison[, plot_eligible :=
  observed_tonnes_yr >= minimum_display_yield &
    simulated_tonnes_yr >= minimum_display_yield]

setcolorder(comparison, c(
  "LME_nb", "ECO_NAME", "group", "observed_period", "simulated_period",
  "simulated_scenario", "observed_tonnes_yr", "simulated_tonnes_yr",
  "ratio_simulated_observed", "matched_cells", "plot_eligible"
))
setorder(comparison, group, LME_nb)
fwrite(comparison, output_table)

plot_data <- comparison[plot_eligible == TRUE]
statistics <- plot_data[, .(
  n_lmes = .N,
  pearson_log10 = cor(
    log10(observed_tonnes_yr),
    log10(simulated_tonnes_yr)
  ),
  spearman = cor(observed_tonnes_yr, simulated_tonnes_yr, method = "spearman")
), by = group]
fwrite(statistics, output_statistics)


# =============================================================================
# 7. CREATE THE FOUR-PANEL SCATTER PLOT
# =============================================================================

group_labels <- c(
  smallPel = "Small pelagics",
  largePel = "Large pelagics",
  demersals = "Demersals",
  total = "Total of the three groups"
)
group_colours <- c(
  smallPel = "#1B9E77",
  largePel = "#D95F02",
  demersals = "#7570B3",
  total = "#1F4E79"
)

plot_data[, group := factor(group, levels = names(group_labels))]

# Use identical limits on both axes so the diagonal truly represents equality.
maximum_yield <- max(
  c(plot_data$observed_tonnes_yr, plot_data$simulated_tonnes_yr),
  na.rm = TRUE
)
axis_maximum <- 10^ceiling(log10(maximum_yield))
axis_limits <- c(minimum_display_yield, axis_maximum)
axis_breaks <- 10^seq(
  floor(log10(minimum_display_yield)),
  ceiling(log10(axis_maximum)),
  by = 2
)

plot_object <- ggplot(
  plot_data,
  aes(
    x = observed_tonnes_yr,
    y = simulated_tonnes_yr,
    colour = group
  )
) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.7, linetype = "dashed") +
  geom_point(size = 2.5, alpha = 0.78) +
  facet_wrap(
    ~group,
    ncol = 2,
    labeller = as_labeller(group_labels)
  ) +
  scale_colour_manual(values = group_colours, guide = "none") +
  scale_x_log10(
    limits = axis_limits,
    breaks = axis_breaks,
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  scale_y_log10(
    limits = axis_limits,
    breaks = axis_breaks,
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  coord_equal() +
  labs(
    title = "Observed and simulated yield across Large Marine Ecosystems",
    subtitle = paste0(
      "Observed mean ", min(observed_years), "-", max(observed_years),
      "; FEISTY mean ", simulated_period, " under F = 0.40 yr-1",
      "\nEach point is one LME; dashed line indicates equal observed and simulated yield"
    ),
    x = "Observed catch (tonnes yr-1, log scale)",
    y = "Simulated FEISTY yield (tonnes yr-1, log scale)",
    caption = paste0(
      "Observed total includes only small pelagics, large pelagics and demersals. ",
      "Values below ", minimum_display_yield, " tonne yr-1 are not displayed."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", colour = "#1F4E79", size = 16),
    plot.subtitle = element_text(size = 10.5, lineheight = 1.15),
    strip.background = element_rect(fill = "grey94", colour = "grey70"),
    strip.text = element_text(face = "bold", colour = "#1F4E79"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0, colour = "grey35", size = 8.5)
  )

ggsave(
  filename = output_figure,
  plot = plot_object,
  width = 11,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("LME comparison complete.")
message("Table: ", output_table)
message("Statistics: ", output_statistics)
message("Figure: ", output_figure)
