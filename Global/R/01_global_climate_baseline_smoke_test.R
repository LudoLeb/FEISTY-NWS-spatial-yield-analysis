# =============================================================================
# 01 - GLOBAL CLIMATE-ONLY BASELINE: TWELVE-CELL SMOKE TEST
# =============================================================================
# Purpose: check the global forcing and FEISTY setup on 12 representative
# cells before any fishing scenario or large calculation is attempted.
#
# Input:  global_feisty_projection.RDS (set FEISTY_GLOBAL_FORCING if needed)
# Output: selected cells, annual F = 0 results, period means and a SSB figure
# Next:   run 02_global_fixed_F_smoke_test.R and then script 03
# =============================================================================

# ---- 1. Packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(FEISTY)
  library(data.table)
  library(ggplot2)
})

# ---- 2. Input file and output folders ---------------------------------------

input_file <- Sys.getenv("FEISTY_GLOBAL_FORCING", unset = "data/raw/global_feisty_projection.RDS")

if (!file.exists(input_file)) {
  stop(
    "Global forcing file not found. Set FEISTY_GLOBAL_FORCING to the path of ",
    "global_feisty_projection.RDS."
  )
}

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ---- 3. Read and validate the forcing ---------------------------------------

message("Reading global forcing: ", input_file)
forcing <- as.data.table(readRDS(input_file))

required_columns <- c(
  "lon", "lat", "time", "bprodin", "zbio", "zprod",
  "Tm", "Tp", "Tb", "depth"
)
missing_columns <- setdiff(required_columns, names(forcing))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

forcing[, time := as.Date(time)]
setorder(forcing, lon, lat, time)

expected_start <- as.Date("2015-01-01")
expected_end <- as.Date("2087-12-01")
if (min(forcing$time) != expected_start || max(forcing$time) != expected_end) {
  warning(
    "Unexpected forcing period: ", min(forcing$time), " to ", max(forcing$time)
  )
}

# ---- 4. Select representative shelf and open-ocean cells --------------------
# Representative targets span shelves, open oceans and climate zones.
targets <- data.table(
  site = c(
    "North Sea shelf", "Bering Sea shelf", "Patagonian shelf",
    "Benguela shelf", "East China Sea shelf", "NW Atlantic shelf",
    "Tropical Atlantic open ocean", "Tropical Pacific open ocean",
    "Tropical Indian open ocean", "North Pacific open ocean",
    "Southern Ocean open ocean", "Arctic open ocean"
  ),
  habitat = rep(c("Shelf", "Open ocean"), each = 6),
  target_lon = c(2, -170, -62, 14, 124, -65, -30, -140, 80, 160, 0, 0),
  target_lat = c(56, 60, -45, -25, 30, 42, 0, 0, -10, 40, -60, 75),
  max_depth = c(rep(500, 6), rep(Inf, 6))
)

grid <- unique(forcing[, .(lon, lat, depth)])

nearest_cell <- function(target_lon, target_lat, max_depth) {
  candidate_grid <- grid[depth <= max_depth]
  if (nrow(candidate_grid) == 0) {
    stop("No candidate cell is available for the requested depth range.")
  }
  # Longitude distance is scaled by latitude for a simple local approximation.
  lon_scale <- max(cos(target_lat * pi / 180), 0.05)
  distance2 <- ((candidate_grid$lon - target_lon) * lon_scale)^2 +
    (candidate_grid$lat - target_lat)^2
  candidate_grid[which.min(distance2)]
}

selected <- rbindlist(lapply(seq_len(nrow(targets)), function(i) {
  cell <- nearest_cell(
    targets$target_lon[i],
    targets$target_lat[i],
    targets$max_depth[i]
  )
  data.table(
    site = targets$site[i],
    habitat = targets$habitat[i],
    target_lon = targets$target_lon[i],
    target_lat = targets$target_lat[i],
    lon = cell$lon,
    lat = cell$lat,
    depth = cell$depth
  )
}))

test_forcing <- merge(
  forcing,
  selected[, .(site, habitat, lon, lat)],
  by = c("lon", "lat"),
  all = FALSE
)
setorder(test_forcing, site, time)

if (any(test_forcing[, .N, by = site]$N != 876L)) {
  stop("At least one selected cell does not have the expected 876 monthly records.")
}

saveRDS(
  test_forcing,
  "data/derived/global_smoke_test_forcing_12_cells.rds",
  compress = "xz"
)
fwrite(selected, "outputs/tables/global_smoke_test_selected_cells.csv")

# ---- 5. FEISTY settings and conversion factors ------------------------------
# Global correction factors retained from the established calibration workflow.
factor_production <- 1.4 / 1.238
factor_bottom_flux <- 5 / 69
factor_small_large <- 0.5
factor_biomass <- 1 / 1.295

target_groups <- c("smallPel", "largePel", "demersals")

# ---- 6. Run FEISTY without fishing in one selected cell ---------------------

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
    Fsmp_ts = rep(0, n_time),
    Fdem_ts = rep(0, n_time),
    Flgp_ts = rep(0, n_time)
  )

  sim <- simulateFEISTY(p)
  output_rows <- 2:(n_time + 1)
  group_names <- p$groupnames[5:(4 + p$nGroups)]

  extract_metric <- function(values, metric_name) {
    values <- as.data.table(values[output_rows, , drop = FALSE])
    setnames(values, group_names)
    values[, `:=`(
      site = position$site[1],
      lon = position$lon[1],
      lat = position$lat[1],
      time = position$time,
      metric = metric_name
    )]
    melt(
      values,
      id.vars = c("site", "lon", "lat", "time", "metric"),
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

message(
  "Running FEISTY at ", uniqueN(test_forcing$site),
  " representative shelf and open-ocean cells"
)
results_monthly <- rbindlist(lapply(split(test_forcing, by = "site"), run_one_site))

results_monthly[, year := as.integer(format(time, "%Y"))]
# ---- 7. Convert monthly results to annual and five-year summaries -----------

results_annual <- results_monthly[, .(
  value = mean(value, na.rm = TRUE)
), by = .(site, lon, lat, year, metric, group)]

commercial_total <- results_annual[, .(
  value = sum(value)
), by = .(site, lon, lat, year, metric)]
commercial_total[, group := "total"]
results_annual <- rbindlist(list(results_annual, commercial_total), use.names = TRUE)

period_summary <- results_annual[
  year %between% c(2015L, 2019L) | year %between% c(2083L, 2087L),
  .(
    value = mean(value, na.rm = TRUE),
    period = ifelse(year <= 2019L, "2015-2019", "2083-2087")
  ),
  by = .(site, lon, lat, metric, group, year)
][, .(value = mean(value)), by = .(site, lon, lat, metric, group, period)]

period_wide <- dcast(
  period_summary,
  site + lon + lat + metric + group ~ period,
  value.var = "value"
)
period_wide[, absolute_change := `2083-2087` - `2015-2019`]
period_wide[, percent_change := fifelse(
  abs(`2015-2019`) > .Machine$double.eps,
  100 * absolute_change / `2015-2019`,
  NA_real_
)]

# ---- 8. Save tables and the diagnostic figure -------------------------------

fwrite(results_annual, "outputs/tables/global_F0_smoke_test_annual_timeseries.csv")
fwrite(period_wide, "outputs/tables/global_F0_smoke_test_period_summary.csv")

plot_data <- results_annual[metric == "SSB" & group %in% target_groups]
mean_line <- plot_data[, .(value = mean(value)), by = .(year, group)]

p <- ggplot(mean_line, aes(year, value, colour = group)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c(
    smallPel = "#0072B2",
    largePel = "#D55E00",
    demersals = "#009E73"
  )) +
  scale_x_continuous(breaks = seq(2020, 2080, 20)) +
  labs(
    title = "Global FEISTY climate-only smoke test",
    subtitle = "Mean SSB across 6 shelf and 6 open-ocean cells; F = 0",
    x = "Year",
    y = expression("SSB (g wet weight " * m^-2 * ")"),
    colour = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  "outputs/figures/global_F0_smoke_test_timeseries.png",
  p,
  width = 8,
  height = 4.8,
  dpi = 300
)

message("Smoke test completed successfully.")
