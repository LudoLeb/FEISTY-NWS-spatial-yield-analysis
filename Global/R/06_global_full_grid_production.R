# =============================================================================
# 06 - COMPLETE GLOBAL PRODUCTION RUN ON THE ONE-DEGREE OCEAN GRID
# =============================================================================
# Purpose: run the final F = 0 and F = 0.40 scenarios on all 41,008 cells.
#
# Input:  global_feisty_projection.RDS, January 2015 to December 2087
# Method: monthly FEISTY simulations, parallel chunks and resumable checkpoints
# Output: cell list, period summaries, weighted annual time series,
#         full-grid 30% SSB0 check and a global diagnostic figure
# Runtime: many hours; use monitor_full_grid_progress.R to follow progress
# Next:   script 07 creates maps and script 08 compares yields among LMEs
# =============================================================================

# ---- 1. Packages, input and runtime settings --------------------------------

suppressPackageStartupMessages({
  library(FEISTY)
  library(data.table)
  library(parallel)
  library(ggplot2)
})

input_file <- Sys.getenv(
  "FEISTY_GLOBAL_FORCING",
  unset = "data/raw/global_feisty_projection.RDS"
)
if (!file.exists(input_file)) {
  stop("Set FEISTY_GLOBAL_FORCING to global_feisty_projection.RDS.")
}

workers_default <- max(1L, min(7L, detectCores(logical = TRUE) - 1L))
workers <- as.integer(Sys.getenv("FEISTY_GLOBAL_WORKERS", unset = workers_default))
chunk_size <- as.integer(Sys.getenv("FEISTY_GLOBAL_CHUNK_SIZE", unset = "70"))
stop_after_chunks <- as.integer(Sys.getenv("FEISTY_GLOBAL_STOP_AFTER_CHUNKS", unset = "0"))

if (is.na(workers) || workers < 1L) stop("FEISTY_GLOBAL_WORKERS must be positive.")
if (is.na(chunk_size) || chunk_size < 1L) stop("FEISTY_GLOBAL_CHUNK_SIZE must be positive.")
if (is.na(stop_after_chunks) || stop_after_chunks < 0L) {
  stop("FEISTY_GLOBAL_STOP_AFTER_CHUNKS cannot be negative.")
}

# ---- 2. Output folders, scenarios and FEISTY factors ------------------------

checkpoint_dir <- "outputs/checkpoints/full_grid"
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

scenarios <- data.table(
  scenario = c("F0", "F040"),
  f_value = c(0, 0.40)
)
target_groups <- c("smallPel", "largePel", "demersals")

# Global FEISTY conversion factors used in the supplied exercise workflow.
factor_production <- 1.4 / 1.238
factor_bottom_flux <- 5 / 69
factor_small_large <- 0.5
factor_biomass <- 1 / 1.295

# ---- 3. Read and validate the complete forcing ------------------------------

message("Reading full global forcing: ", input_file)
forcing <- as.data.table(readRDS(input_file))
forcing[, time := as.Date(time)]
setorder(forcing, lon, lat, depth, time)
forcing[, cell_id := rleid(lon, lat, depth)]

if (min(forcing$time) != as.Date("2015-01-01") ||
    max(forcing$time) != as.Date("2087-12-01")) {
  warning("Unexpected forcing period: ", min(forcing$time), " to ", max(forcing$time))
}

grid <- forcing[!duplicated(cell_id), .(cell_id, lon, lat, depth)]
grid[, habitat := fifelse(depth <= 500, "Shelf", "Open ocean")]
grid[, area_weight := pmax(cos(lat * pi / 180), 0.01)]
grid[, chunk_id := ceiling(cell_id / chunk_size)]
setorder(grid, cell_id)

expected_months <- 876L
if (nrow(forcing) != nrow(grid) * expected_months) {
  stop("The forcing does not contain exactly 876 monthly records per grid cell.")
}
if (uniqueN(forcing$time) != expected_months) {
  stop("The forcing does not contain the expected 876 distinct months.")
}

fwrite(grid, "outputs/tables/global_full_grid_cells.csv")
setkey(forcing, cell_id)
chunk_ids <- sort(unique(grid$chunk_id))

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

run_one_cell <- function(cell_id) {
  position <- copy(chunk_positions[[as.character(cell_id)]])
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
      out[, `:=`(scenario = scenario$scenario, f_value = scenario$f_value)]
      out
    })

    data <- rbindlist(scenario_results, use.names = TRUE)
    data[, `:=`(
      cell_id = cell_id,
      lon = position$lon[1],
      lat = position$lat[1],
      depth = position$depth[1],
      habitat = position$habitat[1],
      area_weight = position$area_weight[1]
    )]
    list(data = data, error = NULL)
  }, error = function(e) {
    list(
      data = NULL,
      error = data.table(
        cell_id = cell_id,
        lon = position$lon[1],
        lat = position$lat[1],
        message = conditionMessage(e)
      )
    )
  })
}

message(
  "Full grid contains ", nrow(grid), " cells in ", length(chunk_ids),
  " chunks; using ", workers, " workers"
)

# ---- 5. Run or resume the parallel chunks -----------------------------------

chunks_completed_this_session <- 0L
for (chunk_id in chunk_ids) {
  checkpoint <- file.path(checkpoint_dir, sprintf("chunk_%04d.rds", chunk_id))
  ids <- grid$cell_id[grid$chunk_id == chunk_id]

  if (file.exists(checkpoint)) {
    valid_checkpoint <- tryCatch({
      existing <- readRDS(checkpoint)
      expected_period_rows <- length(ids) * nrow(scenarios) * 2L * 4L * 2L
      identical(sort(existing$cell_ids), sort(ids)) &&
        setequal(unique(existing$period$cell_id), ids) &&
        nrow(existing$period) == expected_period_rows &&
        nrow(existing$weighted_annual) == nrow(scenarios) * 73L * 2L * 4L &&
        nrow(existing$errors) == 0L &&
        all(is.finite(existing$period$value)) &&
        all(is.finite(existing$weighted_annual$weighted_sum))
    }, error = function(e) FALSE)
    if (valid_checkpoint) {
      message("Chunk ", chunk_id, "/", max(chunk_ids), " already complete")
      next
    }
    invalid_name <- paste0(checkpoint, ".invalid_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    file.rename(checkpoint, invalid_name)
  }

  message(
    "Running chunk ", chunk_id, "/", max(chunk_ids),
    " (", length(ids), " cells)"
  )

  chunk_forcing <- copy(forcing[J(ids), nomatch = 0L])
  metadata_index <- match(chunk_forcing$cell_id, grid$cell_id)
  chunk_forcing[, `:=`(
    habitat = grid$habitat[metadata_index],
    area_weight = grid$area_weight[metadata_index]
  )]
  chunk_positions <- split(chunk_forcing, by = "cell_id")

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

  chunk_annual <- rbindlist(
    lapply(results, `[[`, "data"),
    use.names = TRUE,
    fill = TRUE
  )
  errors <- rbindlist(
    lapply(results, `[[`, "error"),
    use.names = TRUE,
    fill = TRUE
  )

  period <- chunk_annual[
    year %between% c(2015L, 2019L) | year %between% c(2083L, 2087L),
    .(value = mean(value)),
    by = .(
      cell_id, scenario, f_value, lon, lat, depth, habitat, area_weight,
      metric, group,
      period = fifelse(year <= 2019L, "2015-2019", "2083-2087")
    )
  ]
  weighted_annual <- chunk_annual[, .(
    weighted_sum = sum(value * area_weight),
    weight_sum = sum(area_weight)
  ), by = .(scenario, f_value, year, metric, group)]

  checkpoint_object <- list(
    cell_ids = ids,
    period = period,
    weighted_annual = weighted_annual,
    errors = errors
  )
  saveRDS(checkpoint_object, checkpoint, compress = "gzip")

  rm(
    results, checkpoint_object, chunk_annual, errors, period,
    weighted_annual, chunk_forcing, chunk_positions
  )
  invisible(gc())

  chunks_completed_this_session <- chunks_completed_this_session + 1L
  if (stop_after_chunks > 0L &&
      chunks_completed_this_session >= stop_after_chunks) {
    message(
      "Requested validation stop after ", chunks_completed_this_session,
      " newly completed chunk(s)."
    )
    quit(save = "no", status = 0L)
  }
}

checkpoint_files <- file.path(
  checkpoint_dir,
  sprintf("chunk_%04d.rds", chunk_ids)
)
if (any(!file.exists(checkpoint_files))) {
  stop("The full-grid run ended with missing chunk checkpoints.")
}

# ---- 6. Combine and validate every checkpoint -------------------------------
# Release the multi-gigabyte forcing object before consolidating results.
rm(forcing)
invisible(gc())

checkpoint_objects <- lapply(checkpoint_files, readRDS)
errors <- rbindlist(lapply(checkpoint_objects, `[[`, "errors"), fill = TRUE)
if (nrow(errors) > 0L) {
  fwrite(errors, "outputs/tables/global_full_grid_errors.csv")
  warning(nrow(errors), " cells returned errors; see the error table.")
}

# ---- 7. Create cell-period, annual and SSB-constraint summaries -------------

period_summary <- rbindlist(
  lapply(checkpoint_objects, `[[`, "period"),
  use.names = TRUE,
  fill = TRUE
)
setorder(period_summary, scenario, cell_id, metric, group, period)
fwrite(period_summary, "outputs/tables/global_full_grid_period_summary_long.csv")

period_wide <- dcast(
  period_summary,
  cell_id + scenario + f_value + lon + lat + depth + habitat +
    area_weight + metric + group ~ period,
  value.var = "value"
)
period_wide[, absolute_change := `2083-2087` - `2015-2019`]
period_wide[, percent_change := fifelse(
  abs(`2015-2019`) > .Machine$double.eps,
  100 * absolute_change / `2015-2019`,
  NA_real_
)]
fwrite(period_wide, "outputs/tables/global_full_grid_period_summary.csv")

# Keep a much smaller table with only the values needed by the LME comparison.
# This file can be shared on GitHub, unlike the complete period summary.
lme_yield_input <- period_wide[
  scenario == "F040" & metric == "Yield",
  .(lon, lat, group, `2015-2019`)
]
fwrite(
  lme_yield_input,
  "outputs/tables/global_full_grid_yield_2015_2019_for_LME.csv"
)

weighted_chunks <- rbindlist(
  lapply(checkpoint_objects, `[[`, "weighted_annual"),
  use.names = TRUE,
  fill = TRUE
)
global_annual <- weighted_chunks[, .(
  weighted_sum = sum(weighted_sum),
  weight_sum = sum(weight_sum)
), by = .(scenario, f_value, year, metric, group)]
global_annual[, weighted_mean := weighted_sum / weight_sum]
fwrite(global_annual, "outputs/tables/global_full_grid_weighted_annual_summary.csv")

final_ssb <- period_summary[metric == "SSB" & period == "2083-2087"]
f0 <- final_ssb[scenario == "F0", .(
  ssb_F0 = weighted.mean(value, area_weight)
), by = group]
f040 <- final_ssb[scenario == "F040", .(
  ssb_F040 = weighted.mean(value, area_weight)
), by = group]
constraint <- merge(f0, f040, by = "group")
constraint[, ssb_fraction_of_F0 := ssb_F040 / ssb_F0]
constraint[, passes_30pct := ssb_fraction_of_F0 >= 0.30]
fwrite(constraint, "outputs/tables/global_full_grid_F040_SSB_constraint.csv")

# ---- 8. Save the full-grid global time-series figure ------------------------

plot_data <- global_annual[group %in% target_groups]
p <- ggplot(
  plot_data,
  aes(year, weighted_mean, colour = group, linetype = scenario)
) +
  geom_line(linewidth = 0.65) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c(
    smallPel = "#0072B2",
    largePel = "#D55E00",
    demersals = "#009E73"
  )) +
  scale_linetype_manual(values = c(F0 = "solid", F040 = "dashed")) +
  labs(
    title = "FEISTY global full-grid production run",
    subtitle = paste0(nrow(grid), " ocean cells; F = 0 and F = 0.40"),
    x = "Year",
    y = "Area-weighted mean model output",
    colour = "Group",
    linetype = "Scenario"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(
  "outputs/figures/global_full_grid_weighted_timeseries.png",
  p,
  width = 8,
  height = 6.2,
  dpi = 300
)

writeLines(
  paste0(
    "Completed ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "; cells=", nrow(grid),
    "; errors=", nrow(errors)
  ),
  "outputs/checkpoints/full_grid/COMPLETE.txt"
)
message("Full-grid run completed: ", nrow(grid), " cells; ", nrow(errors), " errors")
