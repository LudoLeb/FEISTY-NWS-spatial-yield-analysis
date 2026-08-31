# =============================================================================
# 05 - INTERMEDIATE PRODUCTION RUN ON A REPRESENTATIVE 5-DEGREE GRID
# =============================================================================
# Purpose: validate F = 0 and F = 0.40 on a much larger grid while keeping the
# runtime below that of the complete one-degree calculation.
#
# Method: retain one shelf and one open-ocean cell per 5-degree block
# Input:  global_feisty_projection.RDS
# Output: selected cells, resumable chunks, annual/period summaries,
#         the 30% SSB0 check and a weighted time-series figure
# Next:   run script 06 only after this intermediate calculation succeeds
# =============================================================================

# ---- 1. Packages, input and runtime settings --------------------------------

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

workers_default <- max(1L, min(8L, detectCores(logical = TRUE) - 2L))
workers <- as.integer(Sys.getenv("FEISTY_GLOBAL_WORKERS", unset = workers_default))
chunk_size <- as.integer(Sys.getenv("FEISTY_GLOBAL_CHUNK_SIZE", unset = "40"))
block_size <- 5

if (is.na(workers) || workers < 1L) stop("FEISTY_GLOBAL_WORKERS must be positive.")
if (is.na(chunk_size) || chunk_size < 1L) stop("FEISTY_GLOBAL_CHUNK_SIZE must be positive.")

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/checkpoints/production_5deg", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ---- 2. Scenarios and FEISTY conversion factors -----------------------------

scenarios <- data.table(
  scenario = c("F0", "F040"),
  f_value = c(0, 0.40)
)
target_groups <- c("smallPel", "largePel", "demersals")

factor_production <- 1.4 / 1.238
factor_bottom_flux <- 5 / 69
factor_small_large <- 0.5
factor_biomass <- 1 / 1.295

# ---- 3. Read the forcing and select the 5-degree grid -----------------------

message("Reading global forcing: ", input_file)
forcing <- as.data.table(readRDS(input_file))
forcing[, time := as.Date(time)]
setorder(forcing, lon, lat, time)

if (min(forcing$time) != as.Date("2015-01-01") ||
    max(forcing$time) != as.Date("2087-12-01")) {
  warning("Unexpected forcing period: ", min(forcing$time), " to ", max(forcing$time))
}

# One shelf and one open-ocean representative are retained, when available,
# in each 5-degree block. The nearest cell to the block centre is selected.
grid <- unique(forcing[, .(lon, lat, depth)])
grid[, habitat := fifelse(depth <= 500, "Shelf", "Open ocean")]
grid[, lon_block := floor((lon + 180) / block_size)]
grid[, lat_block := floor((lat + 90) / block_size)]
grid[, block_centre_lon := -180 + (lon_block + 0.5) * block_size]
grid[, block_centre_lat := -90 + (lat_block + 0.5) * block_size]
grid[, distance_to_centre :=
  ((lon - block_centre_lon) * pmax(cos(block_centre_lat * pi / 180), 0.05))^2 +
  (lat - block_centre_lat)^2
]
grid[, stratum_cells := .N, by = .(lon_block, lat_block, habitat)]

setorder(grid, lon_block, lat_block, habitat, distance_to_centre, depth)
selected <- grid[, .SD[1], by = .(lon_block, lat_block, habitat)]
selected[, production_id := seq_len(.N)]
selected[, area_factor := pmax(cos(lat * pi / 180), 0.01)]
selected[, sample_weight := area_factor * stratum_cells]
setorder(selected, production_id)

production_forcing <- merge(
  forcing,
  selected[, .(
    production_id, lon, lat, depth, habitat, lon_block, lat_block,
    stratum_cells, sample_weight
  )],
  by = c("lon", "lat", "depth"),
  all = FALSE
)
setorder(production_forcing, production_id, time)

counts <- production_forcing[, .N, by = production_id]
if (any(counts$N != 876L)) {
  stop("At least one production cell does not contain 876 monthly records.")
}

saveRDS(
  production_forcing,
  "data/derived/global_production_5deg_forcing.rds",
  compress = "gzip"
)
fwrite(selected, "outputs/tables/global_production_5deg_selected_cells.csv")

positions <- split(production_forcing, by = "production_id")
rm(forcing, production_forcing, grid)
invisible(gc())

selected[, chunk_id := ceiling(production_id / chunk_size)]
chunk_ids <- sort(unique(selected$chunk_id))

# ---- 4. Functions used for one FEISTY cell ----------------------------------

extract_annual <- function(values, output_rows, group_names, years, metric_name) {
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

run_one_cell <- function(production_id) {
  position <- copy(positions[[as.character(production_id)]])
  setorder(position, time)
  n_time <- nrow(position)
  years <- as.integer(format(position$time, "%Y"))

  tryCatch({
    scenario_results <- lapply(seq_len(nrow(scenarios)), function(scenario_index) {
      scenario <- scenarios[scenario_index]

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
          Fsmp_ts = rep(scenario$f_value, n_time),
          Fdem_ts = rep(scenario$f_value, n_time),
          Flgp_ts = rep(scenario$f_value, n_time)
        )
        sim <- simulateFEISTY(p)
      }))

      output_rows <- 2:(n_time + 1)
      group_names <- p$groupnames[5:(4 + p$nGroups)]
      out <- rbindlist(list(
        extract_annual(sim$SSB, output_rows, group_names, years, "SSB"),
        extract_annual(sim$yield, output_rows, group_names, years, "Yield")
      ))
      totals <- out[, .(value = sum(value)), by = .(year, metric)]
      totals[, group := "total"]
      out <- rbindlist(list(out, totals), use.names = TRUE)
      out[, `:=`(
        scenario = scenario$scenario,
        f_value = scenario$f_value
      )]
      out
    })

    data <- rbindlist(scenario_results, use.names = TRUE)
    data[, `:=`(
      production_id = production_id,
      lon = position$lon[1],
      lat = position$lat[1],
      depth = position$depth[1],
      habitat = position$habitat[1],
      sample_weight = position$sample_weight[1]
    )]
    list(data = data, error = NULL)
  }, error = function(e) {
    list(
      data = NULL,
      error = data.table(
        production_id = production_id,
        lon = position$lon[1],
        lat = position$lat[1],
        message = conditionMessage(e)
      )
    )
  })
}

message(
  "Selected ", nrow(selected), " production cells in ", length(chunk_ids),
  " chunks; using ", workers, " workers"
)

# ---- 5. Run resumable parallel chunks ---------------------------------------

for (chunk_id in chunk_ids) {
  checkpoint <- file.path(
    "outputs/checkpoints/production_5deg",
    sprintf("chunk_%04d.rds", chunk_id)
  )
  if (file.exists(checkpoint)) {
    message("Chunk ", chunk_id, "/", max(chunk_ids), " already complete")
    next
  }

  ids <- selected$production_id[selected$chunk_id == chunk_id]
  message(
    "Running chunk ", chunk_id, "/", max(chunk_ids),
    " (", length(ids), " cells)"
  )

  active_workers <- min(workers, length(ids))
  if (.Platform$OS.type == "windows" || active_workers == 1L) {
    results <- lapply(ids, run_one_cell)
  } else {
    results <- mclapply(
      ids,
      run_one_cell,
      mc.cores = active_workers,
      mc.preschedule = TRUE
    )
  }

  checkpoint_object <- list(
    data = rbindlist(lapply(results, `[[`, "data"), use.names = TRUE, fill = TRUE),
    errors = rbindlist(lapply(results, `[[`, "error"), use.names = TRUE, fill = TRUE)
  )
  saveRDS(checkpoint_object, checkpoint, compress = "gzip")
  rm(results, checkpoint_object)
  invisible(gc())
}

checkpoint_files <- file.path(
  "outputs/checkpoints/production_5deg",
  sprintf("chunk_%04d.rds", chunk_ids)
)
if (any(!file.exists(checkpoint_files))) {
  stop("The production run ended with missing chunk checkpoints.")
}

checkpoint_objects <- lapply(checkpoint_files, readRDS)
errors <- rbindlist(lapply(checkpoint_objects, `[[`, "errors"), fill = TRUE)
if (nrow(errors) > 0) {
  fwrite(errors, "outputs/tables/global_production_5deg_errors.csv")
  warning(nrow(errors), " cells returned errors; see the error table.")
}

# ---- 6. Combine checkpoints and calculate summaries -------------------------

annual <- rbindlist(lapply(checkpoint_objects, `[[`, "data"), use.names = TRUE)
setorder(annual, scenario, production_id, metric, group, year)
saveRDS(annual, "outputs/tables/global_production_5deg_annual_timeseries.rds", compress = "gzip")

period_summary <- annual[
  year %between% c(2015L, 2019L) | year %between% c(2083L, 2087L),
  .(
    value = mean(value),
    period = fifelse(year <= 2019L, "2015-2019", "2083-2087")
  ),
  by = .(
    production_id, scenario, f_value, lon, lat, depth, habitat,
    sample_weight, metric, group, year
  )
][, .(
  value = mean(value)
), by = .(
  production_id, scenario, f_value, lon, lat, depth, habitat,
  sample_weight, metric, group, period
)]

period_wide <- dcast(
  period_summary,
  production_id + scenario + f_value + lon + lat + depth + habitat +
    sample_weight + metric + group ~ period,
  value.var = "value"
)
period_wide[, absolute_change := `2083-2087` - `2015-2019`]
period_wide[, percent_change := fifelse(
  abs(`2015-2019`) > .Machine$double.eps,
  100 * absolute_change / `2015-2019`,
  NA_real_
)]
fwrite(period_wide, "outputs/tables/global_production_5deg_period_summary.csv")

global_annual <- annual[, .(
  weighted_mean = weighted.mean(value, sample_weight, na.rm = TRUE)
), by = .(scenario, f_value, year, metric, group)]
fwrite(global_annual, "outputs/tables/global_production_5deg_weighted_annual_summary.csv")

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
fwrite(constraint, "outputs/tables/global_production_5deg_F040_SSB_constraint.csv")

# ---- 7. Save the weighted global time-series figure -------------------------

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
    title = "FEISTY global 5-degree production run",
    subtitle = paste0(nrow(selected), " representative cells; F = 0 and F = 0.40"),
    x = "Year",
    y = "Weighted mean model output",
    colour = "Group",
    linetype = "Scenario"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(
  "outputs/figures/global_production_5deg_weighted_timeseries.png",
  p,
  width = 8,
  height = 6.2,
  dpi = 300
)

message(
  "Production run completed: ", nrow(selected), " cells; ",
  nrow(errors), " errors"
)
