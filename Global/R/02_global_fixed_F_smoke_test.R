# =============================================================================
# 02 - FIXED-F FISHING SENSITIVITY: TWELVE-CELL SMOKE TEST
# =============================================================================
# Purpose: add one simple constant fishing scenario to the same 12 cells used
# in script 01 and verify that SSB and yield outputs behave as expected.
#
# Input:  derived forcing and F = 0 results produced by script 01
# Output: F = 0.30 annual results, period means, SSB ratios and a figure
# Note:   F = 0.30 is only an initial sensitivity; script 03 selects F = 0.40
# =============================================================================

# ---- 1. Packages and required inputs ----------------------------------------

suppressPackageStartupMessages({
  library(FEISTY)
  library(data.table)
  library(ggplot2)
})

input_file <- "data/derived/global_smoke_test_forcing_12_cells.rds"
baseline_file <- "outputs/tables/global_F0_smoke_test_annual_timeseries.csv"

if (!file.exists(input_file) || !file.exists(baseline_file)) {
  stop("Run R/01_global_climate_baseline_smoke_test.R first.")
}

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ---- 2. Read the twelve selected cells --------------------------------------

forcing <- as.data.table(readRDS(input_file))
forcing[, time := as.Date(time)]
setorder(forcing, site, time)

# A simple controlled fishing sensitivity. F is an annual instantaneous rate
# and is held constant across time, cells and the three commercial groups.
f_value <- 0.3

# Global correction factors retained from the established calibration workflow.
factor_production <- 1.4 / 1.238
factor_bottom_flux <- 5 / 69
factor_small_large <- 0.5
factor_biomass <- 1 / 1.295
target_groups <- c("smallPel", "largePel", "demersals")

# ---- 3. Run FEISTY with the same fixed F in each commercial group -----------

run_one_site <- function(position) {
  setorder(position, time)
  n_time <- nrow(position)

  p <- setupTimeseries(
    p = setupVertical2(depth = position$depth[1]),
    tStep_ts = 1 / 12,
    tSpin = 40,
    szprod_ts = position$zprod * factor_production * factor_small_large,
    szbio_ts = position$zbio * factor_small_large * factor_biomass,
    lzbio_ts = position$zbio * (1 - factor_small_large) * factor_biomass,
    lzprod_ts = position$zprod * (1 - factor_small_large) * factor_production,
    dfbot = position$bprodin * factor_bottom_flux,
    Tp = position$Tp,
    Tm = position$Tm,
    Tb = position$Tb,
    Fsmp_ts = rep(f_value, n_time),
    Fdem_ts = rep(f_value, n_time),
    Flgp_ts = rep(f_value, n_time)
  )

  sim <- simulateFEISTY(p)
  output_rows <- 2:(n_time + 1)
  group_names <- p$groupnames[5:(4 + p$nGroups)]

  extract_metric <- function(values, metric_name) {
    values <- as.data.table(values[output_rows, , drop = FALSE])
    setnames(values, group_names)
    values[, `:=`(
      site = position$site[1],
      habitat = position$habitat[1],
      lon = position$lon[1],
      lat = position$lat[1],
      time = position$time,
      metric = metric_name
    )]
    melt(
      values,
      id.vars = c("site", "habitat", "lon", "lat", "time", "metric"),
      measure.vars = target_groups,
      variable.name = "group",
      value.name = "value"
    )
  }

  rbindlist(list(
    extract_metric(sim$SSB, "SSB"),
    extract_metric(sim$yield, "Yield")
  ))
}

message("Running fixed-F smoke test with F = ", f_value, " yr-1")
results_monthly <- rbindlist(lapply(split(forcing, by = "site"), run_one_site))
results_monthly[, year := as.integer(format(time, "%Y"))]

# ---- 4. Calculate annual results and comparison periods ---------------------

results_annual <- results_monthly[, .(
  value = mean(value, na.rm = TRUE)
), by = .(site, habitat, lon, lat, year, metric, group)]

commercial_total <- results_annual[, .(
  value = sum(value)
), by = .(site, habitat, lon, lat, year, metric)]
commercial_total[, group := "total"]
results_annual <- rbindlist(list(results_annual, commercial_total), use.names = TRUE)
results_annual[, scenario := "Constant F = 0.3"]

baseline <- fread(baseline_file)[metric == "SSB"]
baseline[, scenario := "F = 0"]
baseline <- baseline[, .(site, lon, lat, year, metric, group, value, scenario)]

ssb_comparison <- rbindlist(list(
  baseline,
  results_annual[metric == "SSB", .(site, lon, lat, year, metric, group, value, scenario)]
))

fished_period <- results_annual[
  year %between% c(2015L, 2019L) | year %between% c(2083L, 2087L),
  .(
    value = mean(value),
    period = ifelse(year <= 2019L, "2015-2019", "2083-2087")
  ),
  by = .(site, habitat, lon, lat, metric, group, year)
][, .(value = mean(value)), by = .(site, habitat, lon, lat, metric, group, period)]

fished_period <- dcast(
  fished_period,
  site + habitat + lon + lat + metric + group ~ period,
  value.var = "value"
)
fished_period[, absolute_change := `2083-2087` - `2015-2019`]
fished_period[, percent_change := fifelse(
  abs(`2015-2019`) > .Machine$double.eps,
  100 * absolute_change / `2015-2019`,
  NA_real_
)]

# ---- 5. Check the 30% SSB0 threshold ----------------------------------------
# Evaluate the threshold against the climate-only reference in each
# year and cell. The total is the sum of the three commercial groups.
ssb_fished <- results_annual[metric == "SSB"]
ssb_zero <- fread(baseline_file)[metric == "SSB"]
ssb_ratio <- merge(
  ssb_fished,
  ssb_zero[, .(site, year, group, ssb_F0 = value)],
  by = c("site", "year", "group"),
  all.x = TRUE
)
ssb_ratio[, ssb_fraction_of_F0 := fifelse(
  ssb_F0 > .Machine$double.eps,
  value / ssb_F0,
  NA_real_
)]

# ---- 6. Save tables and diagnostic plots ------------------------------------

fwrite(results_annual, "outputs/tables/global_F030_smoke_test_annual_timeseries.csv")
fwrite(fished_period, "outputs/tables/global_F030_smoke_test_period_summary.csv")
fwrite(ssb_ratio, "outputs/tables/global_F030_smoke_test_SSB_fraction_of_F0.csv")

mean_ssb <- ssb_comparison[group %in% target_groups, .(
  value = mean(value)
), by = .(year, group, scenario)]

mean_yield <- results_annual[
  metric == "Yield" & group %in% target_groups,
  .(value = mean(value)),
  by = .(year, group)
]

colours <- c(
  smallPel = "#0072B2",
  largePel = "#D55E00",
  demersals = "#009E73"
)

p_ssb <- ggplot(mean_ssb, aes(year, value, colour = group, linetype = scenario)) +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = colours) +
  scale_linetype_manual(values = c("F = 0" = "solid", "Constant F = 0.3" = "dashed")) +
  labs(x = NULL, y = expression("SSB (g wet weight " * m^-2 * ")"), linetype = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_yield <- ggplot(mean_yield, aes(year, value, colour = group)) +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = colours) +
  labs(x = "Year", y = expression("Yield (g wet weight " * m^-2 * yr^-1 * ")")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

png(
  "outputs/figures/global_F030_smoke_test_SSB_and_yield.png",
  width = 2400,
  height = 2100,
  res = 300
)
grid::grid.newpage()
layout <- grid::grid.layout(nrow = 2, ncol = 1)
grid::pushViewport(grid::viewport(layout = layout))
print(p_ssb, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_yield, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

message("Fixed-F smoke test completed successfully.")
