# =============================================================================
# 04 - STRATIFIED GLOBAL PILOT: F = 0 AND F = 0.40
# =============================================================================
# Purpose: validate the selected fishing scenario on a larger and spatially
# balanced sample before running a production grid.
#
# Method: select shelf/open-ocean cells across latitude and longitude strata
# Input:  global_feisty_projection.RDS
# Output: selected cells, resumable checkpoints, annual/period results,
#         the 30% SSB0 check and a weighted global time series
# Result: 134 cells and 268 successful scenario simulations
# =============================================================================

# ---- 1. Packages, input and output folders ----------------------------------

suppressPackageStartupMessages({
  library(FEISTY)
  library(data.table)
  library(parallel)
  library(ggplot2)
})

input_file <- Sys.getenv("FEISTY_GLOBAL_FORCING", unset = "data/raw/global_feisty_projection.RDS")
if (!file.exists(input_file)) {
  stop("Set FEISTY_GLOBAL_FORCING to global_feisty_projection.RDS.")
}

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/checkpoints/pilot", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ---- 2. Scenarios and FEISTY settings ---------------------------------------

scenarios <- data.table(
  scenario = c("F0", "F040"),
  f_value = c(0, 0.40)
)
cells_per_stratum <- 2L
random_seed <- 20260825L
target_groups <- c("smallPel", "largePel", "demersals")

factor_production <- 1.4 / 1.238
factor_bottom_flux <- 5 / 69
factor_small_large <- 0.5
factor_biomass <- 1 / 1.295

# ---- 3. Read the forcing and select a stratified sample ---------------------

message("Reading global forcing")
forcing <- as.data.table(readRDS(input_file))
forcing[, time := as.Date(time)]
setorder(forcing, lon, lat, time)

grid <- unique(forcing[, .(lon, lat, depth)])
grid[, habitat := fifelse(depth <= 500, "Shelf", "Open ocean")]
grid[, lat_band := cut(
  lat,
  breaks = c(-90, -60, -30, 0, 30, 60, 90),
  include.lowest = TRUE,
  right = FALSE
)]
grid[, lon_sector := cut(
  lon,
  breaks = c(-180, -120, -60, 0, 60, 120, 180),
  include.lowest = TRUE,
  right = FALSE
)]
grid[, stratum := interaction(habitat, lat_band, lon_sector, drop = TRUE)]
grid[, stratum_cells := .N, by = stratum]

set.seed(random_seed)
selected <- grid[, .SD[sample.int(.N, min(.N, cells_per_stratum))], by = stratum]
selected[, pilot_id := seq_len(.N)]
selected[, selected_in_stratum := .N, by = stratum]
selected[, area_factor := pmax(cos(lat * pi / 180), 0.01)]
selected[, sample_weight := area_factor * stratum_cells / selected_in_stratum]
setorder(selected, pilot_id)

pilot_forcing <- merge(
  forcing,
  selected[, .(
    pilot_id, lon, lat, depth, habitat, lat_band, lon_sector,
    stratum, stratum_cells, selected_in_stratum, sample_weight
  )],
  by = c("lon", "lat", "depth"),
  all = FALSE
)
setorder(pilot_forcing, pilot_id, time)

counts <- pilot_forcing[, .N, by = pilot_id]
if (any(counts$N != 876L)) {
  stop("At least one pilot cell does not contain 876 monthly records.")
}

saveRDS(
  pilot_forcing,
  "data/derived/global_pilot_forcing_stratified.rds",
  compress = "xz"
)
fwrite(selected, "outputs/tables/global_pilot_selected_cells.csv")

positions <- split(pilot_forcing, by = "pilot_id")
rm(forcing, pilot_forcing, grid)
invisible(gc())

jobs <- rbindlist(lapply(seq_len(nrow(scenarios)), function(i) {
  data.table(
    pilot_id = as.integer(names(positions)),
    scenario = scenarios$scenario[i],
    f_value = scenarios$f_value[i]
  )
}))
setorder(jobs, scenario, pilot_id)
jobs[, checkpoint := file.path(
  "outputs/checkpoints/pilot",
  sprintf("%s_cell_%04d.rds", scenario, pilot_id)
)]
jobs[, complete := file.exists(checkpoint)]

# ---- 4. Run one cell-scenario job and save its checkpoint -------------------

run_one_job <- function(job_index) {
  job <- jobs[job_index]
  if (job$complete) {
    return(job$checkpoint)
  }

  position <- copy(positions[[as.character(job$pilot_id)]])
  setorder(position, time)
  n_time <- nrow(position)

  invisible(capture.output({
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
      Fsmp_ts = rep(job$f_value, n_time),
      Fdem_ts = rep(job$f_value, n_time),
      Flgp_ts = rep(job$f_value, n_time)
    )
    sim <- simulateFEISTY(p)
  }))

  output_rows <- 2:(n_time + 1)
  group_names <- p$groupnames[5:(4 + p$nGroups)]
  years <- as.integer(format(position$time, "%Y"))

  extract_annual <- function(values, metric_name) {
    values <- as.data.table(values[output_rows, , drop = FALSE])
    setnames(values, group_names)
    values[, year := years]
    annual <- values[, lapply(.SD, mean), by = year, .SDcols = target_groups]
    long <- melt(
      annual,
      id.vars = "year",
      measure.vars = target_groups,
      variable.name = "group",
      value.name = "value"
    )
    long[, metric := metric_name]
    long
  }

  out <- rbindlist(list(
    extract_annual(sim$SSB, "SSB"),
    extract_annual(sim$yield, "Yield")
  ))

  totals <- out[, .(value = sum(value)), by = .(year, metric)]
  totals[, group := "total"]
  out <- rbindlist(list(out, totals), use.names = TRUE)
  out[, `:=`(
    pilot_id = job$pilot_id,
    scenario = job$scenario,
    f_value = job$f_value,
    lon = position$lon[1],
    lat = position$lat[1],
    depth = position$depth[1],
    habitat = position$habitat[1],
    stratum = as.character(position$stratum[1]),
    sample_weight = position$sample_weight[1]
  )]

  saveRDS(out, job$checkpoint, compress = "xz")
  job$checkpoint
}

pending <- which(!jobs$complete)
if (length(pending) > 0) {
  workers <- min(4L, detectCores(logical = TRUE), length(pending))
  message(
    "Running ", length(pending), " pilot jobs across ", nrow(selected),
    " cells with ", workers, " workers"
  )
  if (.Platform$OS.type == "windows" || workers == 1L) {
    invisible(lapply(pending, run_one_job))
  } else {
    invisible(mclapply(
      pending,
      run_one_job,
      mc.cores = workers,
      mc.preschedule = TRUE
    ))
  }
} else {
  message("All pilot checkpoints are already available.")
}

missing_checkpoints <- jobs[!file.exists(checkpoint)]
if (nrow(missing_checkpoints) > 0) {
  stop(nrow(missing_checkpoints), " pilot jobs did not produce checkpoints.")
}

# ---- 5. Combine checkpoints and calculate final summaries -------------------

annual <- rbindlist(lapply(jobs$checkpoint, readRDS), use.names = TRUE)
setorder(annual, scenario, pilot_id, metric, group, year)
fwrite(annual, "outputs/tables/global_pilot_annual_timeseries.csv")

period_summary <- annual[
  year %between% c(2015L, 2019L) | year %between% c(2083L, 2087L),
  .(
    value = mean(value),
    period = fifelse(year <= 2019L, "2015-2019", "2083-2087")
  ),
  by = .(
    pilot_id, scenario, f_value, lon, lat, depth, habitat,
    stratum, sample_weight, metric, group, year
  )
][, .(
  value = mean(value)
), by = .(
  pilot_id, scenario, f_value, lon, lat, depth, habitat,
  stratum, sample_weight, metric, group, period
)]

period_wide <- dcast(
  period_summary,
  pilot_id + scenario + f_value + lon + lat + depth + habitat +
    stratum + sample_weight + metric + group ~ period,
  value.var = "value"
)
period_wide[, absolute_change := `2083-2087` - `2015-2019`]
period_wide[, percent_change := fifelse(
  abs(`2015-2019`) > .Machine$double.eps,
  100 * absolute_change / `2015-2019`,
  NA_real_
)]
fwrite(period_wide, "outputs/tables/global_pilot_period_summary.csv")

global_annual <- annual[, .(
  weighted_mean = weighted.mean(value, sample_weight, na.rm = TRUE)
), by = .(scenario, f_value, year, metric, group)]
fwrite(global_annual, "outputs/tables/global_pilot_weighted_annual_summary.csv")

final_ssb <- period_summary[metric == "SSB" & period == "2083-2087"]
f0 <- final_ssb[scenario == "F0", .(
  ssb_F0 = weighted.mean(value, sample_weight)
), by = group]
f040 <- final_ssb[scenario == "F040", .(
  ssb_F040 = weighted.mean(value, sample_weight)
), by = group]
constraint <- merge(f0, f040, by = "group")
constraint[, ssb_fraction_of_F0 := ssb_F040 / ssb_F0]
constraint[, passes_30pct := ssb_fraction_of_F0 >= 0.30]
fwrite(constraint, "outputs/tables/global_pilot_F040_SSB_constraint.csv")

# ---- 6. Save the weighted global time-series figure -------------------------

plot_data <- global_annual[group %in% target_groups]
p <- ggplot(plot_data, aes(year, weighted_mean, colour = group, linetype = scenario)) +
  geom_line(linewidth = 0.65) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c(
    smallPel = "#0072B2",
    largePel = "#D55E00",
    demersals = "#009E73"
  )) +
  scale_linetype_manual(values = c(F0 = "solid", F040 = "dashed")) +
  labs(
    title = "FEISTY global stratified pilot",
    subtitle = paste0(nrow(selected), " cells; solid F = 0, dashed F = 0.40"),
    x = "Year",
    y = "Weighted mean model output",
    colour = "Group",
    linetype = "Scenario"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(
  "outputs/figures/global_pilot_weighted_timeseries.png",
  p,
  width = 8,
  height = 6.2,
  dpi = 300
)

message("Pilot completed: ", nrow(selected), " cells and ", nrow(jobs), " jobs")
