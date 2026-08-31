# =============================================================================
# 03 - SELECT A CONSTANT GLOBAL F ON THE TWELVE-CELL TEST GRID
# =============================================================================
# Purpose: test F = 0.05 to 0.60 and retain the highest value for which all
# three commercial groups keep at least 30% of their no-fishing SSB.
#
# Input:  derived forcing and F = 0 results produced by script 01
# Output: results by site/group, aggregated SSB checks and an F-yield figure
# Result: F = 0.40 passes; F = 0.45 fails for large pelagics
# Next:   validate F = 0.40 with scripts 04, 05 and 06
# =============================================================================

# ---- 1. Packages and required inputs ----------------------------------------

suppressPackageStartupMessages({
  library(FEISTY)
  library(data.table)
  library(parallel)
  library(ggplot2)
})

forcing_file <- "data/derived/global_smoke_test_forcing_12_cells.rds"
baseline_file <- "outputs/tables/global_F0_smoke_test_annual_timeseries.csv"

if (!file.exists(forcing_file) || !file.exists(baseline_file)) {
  stop("Run R/01_global_climate_baseline_smoke_test.R first.")
}

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

forcing <- as.data.table(readRDS(forcing_file))
forcing[, time := as.Date(time)]
setorder(forcing, site, time)

# ---- 2. Values of F and biological constraint -------------------------------

f_values <- seq(0.05, 0.60, by = 0.05)
target_groups <- c("smallPel", "largePel", "demersals")
final_years <- 2083:2087
ssb_threshold <- 0.30

# Very small cell-level reference biomasses create unstable ratios while
# contributing negligibly to the aggregate. Cell-level checks are retained as
# diagnostics, but the 30% decision is made on aggregated SSB, consistently
# with the regional NWS analysis.
minimum_reference_ssb <- 0.01

factor_production <- 1.4 / 1.238
factor_bottom_flux <- 5 / 69
factor_small_large <- 0.5
factor_biomass <- 1 / 1.295

# ---- 3. Prepare the final-period no-fishing reference -----------------------

baseline <- fread(baseline_file)[
  metric == "SSB" & group %in% target_groups & year %in% final_years,
  .(ssb_F0 = mean(value)),
  by = .(site, group)
]

positions <- split(forcing, by = "site")
jobs <- CJ(f_value = f_values, site = names(positions), sorted = TRUE)
results_file <- "outputs/tables/global_constant_F_scan_site_group_results.csv"

existing_results <- NULL
if (file.exists(results_file)) {
  previous <- fread(results_file)[f_value %in% f_values]
  required_previous <- c(
    "f_value", "site", "habitat", "lon", "lat", "group",
    "final_ssb", "final_yield"
  )
  if (all(required_previous %in% names(previous))) {
    existing_results <- previous[, ..required_previous]
    completed_jobs <- unique(existing_results[, .(f_value, site)])
    jobs <- jobs[!completed_jobs, on = .(f_value, site)]
  }
}

# ---- 4. Run one combination of cell and fishing pressure --------------------

run_one_job <- function(job_index) {
  f_value <- jobs$f_value[job_index]
  site_name <- jobs$site[job_index]
  position <- copy(positions[[site_name]])
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

  invisible(capture.output(sim <- simulateFEISTY(p)))
  output_rows <- 2:(n_time + 1)
  final_rows <- output_rows[as.integer(format(position$time, "%Y")) %in% final_years]
  group_names <- p$groupnames[5:(4 + p$nGroups)]

  ssb <- colMeans(sim$SSB[final_rows, , drop = FALSE])
  yield <- colMeans(sim$yield[final_rows, , drop = FALSE])
  names(ssb) <- group_names
  names(yield) <- group_names

  data.table(
    f_value = f_value,
    site = site_name,
    habitat = position$habitat[1],
    lon = position$lon[1],
    lat = position$lat[1],
    group = target_groups,
    final_ssb = ssb[target_groups],
    final_yield = yield[target_groups]
  )
}

if (nrow(jobs) > 0) {
  workers <- min(4L, detectCores(logical = TRUE), nrow(jobs))
  message(
    "Running ", nrow(jobs), " missing F-by-site jobs using ",
    workers, " workers"
  )

  if (.Platform$OS.type == "windows" || workers == 1L) {
    new_results <- rbindlist(lapply(seq_len(nrow(jobs)), run_one_job))
  } else {
    new_results <- rbindlist(mclapply(
      seq_len(nrow(jobs)),
      run_one_job,
      mc.cores = workers,
      mc.preschedule = TRUE
    ))
  }
} else {
  message("All requested F-by-site jobs are already available.")
  new_results <- NULL
}

scan_results <- rbindlist(
  Filter(Negate(is.null), list(existing_results, new_results)),
  use.names = TRUE
)

scan_results <- merge(
  scan_results,
  baseline,
  by = c("site", "group"),
  all.x = TRUE
)
scan_results[, ssb_fraction_of_F0 := fifelse(
  ssb_F0 > .Machine$double.eps,
  final_ssb / ssb_F0,
  NA_real_
)]
scan_results[, threshold_evaluable := ssb_F0 >= minimum_reference_ssb]
scan_results[, threshold_pass := fifelse(
  threshold_evaluable,
  ssb_fraction_of_F0 >= ssb_threshold,
  NA
)]

# ---- 5. Aggregate SSB, apply the threshold and select F ---------------------

group_summary <- scan_results[, .(
  aggregate_ssb_fraction = sum(final_ssb) / sum(ssb_F0),
  aggregate_yield = sum(final_yield)
), by = .(f_value, group)]

threshold_summary <- group_summary[, .(
  min_group_ssb_fraction = min(aggregate_ssb_fraction),
  limiting_group = group[which.min(aggregate_ssb_fraction)],
  all_groups_pass = all(aggregate_ssb_fraction >= ssb_threshold)
), by = f_value]

yield_summary <- scan_results[, .(
  total_yield = sum(final_yield)
), by = .(f_value, site)][, .(
  mean_total_yield = mean(total_yield)
), by = f_value]

scan_summary <- merge(threshold_summary, yield_summary, by = "f_value")

eligible <- scan_summary[all_groups_pass == TRUE]
recommended_f <- if (nrow(eligible) > 0) {
  eligible$f_value[which.max(eligible$mean_total_yield)]
} else {
  NA_real_
}
scan_summary[, recommended_on_smoke_test := f_value == recommended_f]

# ---- 6. Save tables and the selection figure --------------------------------

fwrite(scan_results, results_file)
fwrite(group_summary, "outputs/tables/global_constant_F_scan_aggregate_group_results.csv")
fwrite(scan_summary, "outputs/tables/global_constant_F_scan_summary.csv")

p <- ggplot(scan_summary, aes(f_value, mean_total_yield)) +
  geom_line(colour = "#0072B2", linewidth = 0.8) +
  geom_point(aes(fill = all_groups_pass), shape = 21, size = 3, colour = "black") +
  scale_fill_manual(
    values = c(`TRUE` = "#009E73", `FALSE` = "#D55E00"),
    labels = c(`TRUE` = "All aggregate groups above 30% SSB0", `FALSE` = "At least one aggregate group below 30% SSB0")
  ) +
  scale_x_continuous(breaks = f_values) +
  labs(
    title = "Constant-F smoke-test scan",
    subtitle = "Mean final-period yield across 6 shelf and 6 open-ocean cells",
    x = expression("Constant fishing mortality F (yr"^-1 * ")"),
    y = expression("Mean total yield (g wet weight " * m^-2 * yr^-1 * ")"),
    fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(
  "outputs/figures/global_constant_F_scan_smoke_test.png",
  p,
  width = 8,
  height = 5,
  dpi = 300
)

if (is.na(recommended_f)) {
  message("No scanned F value satisfied the aggregate 30% SSB0 threshold for every group.")
} else {
  message("Recommended smoke-test constant F: ", recommended_f, " yr-1")
}
