# ============================================================
# FEISTY - FULL NWS GRID: PRESENT-FUTURE YIELD CHANGE
#
# Objective
# ---------
# Map the spatial change in yield between 2020-2029 and
# 2090-2100 under the best valid regional fishing strategy.
#
# The regional strategies come from:
# NWS_decadal_optimisation_best_approach.csv
#
# Each selected strategy has already passed the regional rule:
# all three fished groups remain at or above 30% of SSB0.
#
# Outputs
# -------
# Absolute and percentage yield change for:
# - total yield (the sum of the three groups below)
# - smallPel
# - demersals
# - largePel
#
# Important interpretation
# ------------------------
# These maps show the local contribution to yield under the
# regionally selected maximum-valid strategy. They are not a
# separate cell-by-cell optimisation of fishing mortality.
#
# Full-grid constraint check
# --------------------------
# A third, unfished simulation is run for every cell. This allows
# the 30% SSB0 rule to be recalculated on the full spatial grid.
# Maps are produced only if all regional full-grid ratios pass.
#
# Runtime safety
# --------------
# The default is a small test. Change full_run to TRUE only after
# the test completes successfully. Full results are checkpointed
# by chunk and can be resumed after an interruption.
# ============================================================


# ============================================================
# 0. RUN SETTINGS
# ============================================================

full_run <- FALSE

test_cells_per_region <- 1L
chunk_size <- 200L

detected_cores <- parallel::detectCores(logical = FALSE)

if (is.na(detected_cores)) {
  detected_cores <- 2L
}

n_cores <- max(
  1L,
  min(8L, detected_cores - 2L)
)


# ============================================================
# 1. LIBRARIES
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(FEISTY)
  library(maps)
})


# ============================================================
# 2. PATHS
# ============================================================

# Run this script from the repository root (opening FEISTY_NWS.Rproj
# in RStudio does this automatically).
base_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- file.path(base_dir, "data")
outputs_dir <- file.path(base_dir, "outputs")
rds_dir <- file.path(outputs_dir, "rds")
tables_dir <- file.path(outputs_dir, "tables")
figures_dir <- file.path(outputs_dir, "figures")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

environment_file <- file.path(
  data_dir,
  "nws_for_feisty_proj_year.RDS"
)

ffmsy_file <- file.path(
  data_dir,
  "all_FFMSY.RDS"
)

classified_cells_file <- file.path(
  rds_dir,
  "NWS_full_grid_four_regions_cells.RDS"
)

best_strategy_file <- file.path(
  tables_dir,
  "NWS_decadal_optimisation_best_approach.csv"
)

required_files <- c(
  environment_file,
  ffmsy_file,
  classified_cells_file,
  best_strategy_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n- ",
    paste(missing_files, collapse = "\n- ")
  )
}


# ============================================================
# 3. ANALYSIS PARAMETERS
# ============================================================

factorProduction <- 1.4 / 1.238
factorBottomflux <- 5 / 69
factorSmallLarge <- 0.5
factorBiomass <- 1 / 1.295

temperature_slope <- 0.063
SSB_threshold <- 0.30

target_groups <- c(
  "smallPel",
  "demersals",
  "largePel"
)

study_regions <- c(
  "Bay of Biscay",
  "Celtic Sea",
  "English Channel",
  "North Sea"
)

period_definition <- tribble(
  ~decade,      ~climate_period,          ~start_year, ~end_year,
  "2020-2029", "Reference 2020-2029",          2020L,      2029L,
  "2090-2100", "Future 2090-2100",             2090L,      2100L
)


# ============================================================
# 4. SELECT AND VERIFY THE REGIONAL STRATEGIES
# ============================================================

best_strategy <- read.csv(
  best_strategy_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  as_tibble() %>%
  mutate(
    region = as.character(region),
    decade = as.character(decade),
    best_approach = as.character(best_approach),
    optimal_F_or_Fcoef = as.numeric(optimal_F_or_Fcoef),
    minimum_SSB_ratio = as.numeric(minimum_SSB_ratio)
  ) %>%
  inner_join(period_definition, by = "decade") %>%
  filter(region %in% study_regions) %>%
  arrange(match(region, study_regions), start_year)


# ============================================================
# FINAL BAY OF BISCAY STRATEGY
# ============================================================
# Values validated on the spatial grid:
# 2020-2029: constant F = 0.33
# 2090-2100: FFMSY temperature coefficient = 0.54
best_strategy <- best_strategy %>%
  mutate(
    optimal_F_or_Fcoef = case_when(
      region == "Bay of Biscay" &
        decade == "2020-2029" ~ 0.33,

      region == "Bay of Biscay" &
        decade == "2090-2100" ~ 0.54,

      TRUE ~ optimal_F_or_Fcoef
    )
  )

required_strategy_columns <- c(
  "region",
  "decade",
  "best_approach",
  "optimal_F_or_Fcoef",
  "minimum_SSB_ratio",
  "climate_period",
  "start_year",
  "end_year"
)

missing_strategy_columns <- setdiff(
  required_strategy_columns,
  names(best_strategy)
)

if (length(missing_strategy_columns) > 0) {
  stop(
    "Missing strategy columns: ",
    paste(missing_strategy_columns, collapse = ", ")
  )
}

if (nrow(best_strategy) != length(study_regions) * 2L) {
  stop(
    "Expected exactly eight selected regional strategies but found ",
    nrow(best_strategy),
    "."
  )
}

if (anyDuplicated(best_strategy[c("region", "decade")])) {
  stop("The selected strategy table has duplicated region-decade rows.")
}

if (any(!best_strategy$best_approach %in% c("constant_F", "FFMSY_temp"))) {
  stop("An unknown fishing approach is present in the strategy table.")
}

if (any(best_strategy$minimum_SSB_ratio < SSB_threshold, na.rm = TRUE)) {
  stop("At least one selected strategy violates the 30% SSB0 rule.")
}

cat("\n========================================\n")
cat("SELECTED REGIONAL STRATEGIES\n")
cat("========================================\n\n")

print(
  as.data.frame(
    best_strategy %>%
      select(
        region,
        decade,
        best_approach,
        optimal_F_or_Fcoef,
        minimum_SSB_ratio,
        limiting_group
      )
  ),
  row.names = FALSE
)

write.csv(
  best_strategy,
  file.path(
    tables_dir,
    "NWS_full_grid_selected_regional_strategies.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 5. LOAD AND SELECT SPATIAL CELLS
# ============================================================

full_cells <- readRDS(classified_cells_file) %>%
  as_tibble() %>%
  transmute(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    sea = as.character(sea)
  ) %>%
  filter(sea %in% study_regions) %>%
  arrange(match(sea, study_regions), lat, lon)

if (nrow(full_cells) != 22317L) {
  warning(
    "The historical-box classification previously contained 22,317 cells; ",
    "the current file contains ",
    format(nrow(full_cells), big.mark = ","),
    "."
  )
}

if (full_run) {

  # FINAL ANALYSIS: every available cell in the four study regions
  analysis_cells <- full_cells

  run_label <- "FULL_GRID_22317_CELLS_BAY_FINAL"

} else {

  analysis_cells <- full_cells %>%
    group_by(sea) %>%
    slice(round(seq(1, n(), length.out = test_cells_per_region))) %>%
    ungroup()

  run_label <- "TEST"
}

analysis_cells <- analysis_cells %>%
  mutate(cell_id = row_number()) %>%
  select(cell_id, everything())

cat("\nRun mode:", run_label, "\n")
cat("Cells:", format(nrow(analysis_cells), big.mark = ","), "\n")
cat("FEISTY simulations:", format(nrow(analysis_cells) * 3L, big.mark = ","), "\n")
cat("Parallel workers:", n_cores, "\n")

print(
  as.data.frame(
    analysis_cells %>%
      count(sea, name = "number_of_cells")
  ),
  row.names = FALSE
)


# ============================================================
# 6. PREPARE THE LATEST FFMSY GRID
# ============================================================

ffmsy_long <- readRDS(ffmsy_file)
latest_ffmsy_year <- max(ffmsy_long$year, na.rm = TRUE)

ffmsy_reference <- ffmsy_long %>%
  filter(
    year == latest_ffmsy_year,
    Lat >= 39,
    Lat <= 68,
    Lon >= -21,
    Lon <= 16
  ) %>%
  select(Lat, Lon, group, value) %>%
  group_by(Lat, Lon, group) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = group, values_from = value) %>%
  mutate(
    lon_index = round((Lon - 0.25) / 0.5),
    lat_index = round((Lat - 0.25) / 0.5),
    grid_key = paste(lon_index, lat_index, sep = ":")
  ) %>%
  as.data.frame()

rm(ffmsy_long)
invisible(gc())

missing_ffmsy_groups <- setdiff(target_groups, names(ffmsy_reference))

if (length(missing_ffmsy_groups) > 0) {
  stop(
    "Missing FFMSY groups: ",
    paste(missing_ffmsy_groups, collapse = ", ")
  )
}


# ============================================================
# 7. MATCH EACH NWS CELL TO THE NEAREST FFMSY CELL
#
# The FFMSY grid is regular at 0.5 degrees. Candidate cells are
# searched within two degrees of the rounded grid position.
# ============================================================

ffmsy_lookup <- setNames(
  seq_len(nrow(ffmsy_reference)),
  ffmsy_reference$grid_key
)

candidate_offsets <- expand.grid(
  lon_offset = -4:4,
  lat_offset = -4:4
)

nearest_ffmsy_index <- integer(nrow(analysis_cells))
nearest_ffmsy_distance <- numeric(nrow(analysis_cells))

for (i in seq_len(nrow(analysis_cells))) {

  base_lon_index <- round(
    (analysis_cells$lon[i] - 0.25) / 0.5
  )

  base_lat_index <- round(
    (analysis_cells$lat[i] - 0.25) / 0.5
  )

  candidate_keys <- paste(
    base_lon_index + candidate_offsets$lon_offset,
    base_lat_index + candidate_offsets$lat_offset,
    sep = ":"
  )

  candidate_indices <- unname(ffmsy_lookup[candidate_keys])
  candidate_indices <- candidate_indices[!is.na(candidate_indices)]

  if (length(candidate_indices) == 0) {
    stop(
      "No nearby FFMSY cell found for lon = ",
      analysis_cells$lon[i],
      ", lat = ",
      analysis_cells$lat[i]
    )
  }

  dx <- (
    ffmsy_reference$Lon[candidate_indices] - analysis_cells$lon[i]
  ) * cos(analysis_cells$lat[i] * pi / 180)

  dy <-
    ffmsy_reference$Lat[candidate_indices] - analysis_cells$lat[i]

  candidate_distance_squared <- dx^2 + dy^2
  selected_candidate <- which.min(candidate_distance_squared)

  nearest_ffmsy_index[i] <- candidate_indices[selected_candidate]
  nearest_ffmsy_distance[i] <- sqrt(
    candidate_distance_squared[selected_candidate]
  )
}

analysis_cells <- analysis_cells %>%
  mutate(
    lon_ffmsy = ffmsy_reference$Lon[nearest_ffmsy_index],
    lat_ffmsy = ffmsy_reference$Lat[nearest_ffmsy_index],
    ffmsy_demersals = ffmsy_reference$demersals[nearest_ffmsy_index],
    ffmsy_smallPel = ffmsy_reference$smallPel[nearest_ffmsy_index],
    ffmsy_largePel = ffmsy_reference$largePel[nearest_ffmsy_index],
    ffmsy_distance = nearest_ffmsy_distance
  )

if (any(!complete.cases(
  analysis_cells[
    c("ffmsy_demersals", "ffmsy_smallPel", "ffmsy_largePel")
  ]
))) {
  stop("At least one cell has an incomplete FFMSY match.")
}

cat("\nFFMSY reference year:", latest_ffmsy_year, "\n")
cat(
  "Maximum approximate matching distance:",
  round(max(analysis_cells$ffmsy_distance), 4),
  "degrees\n"
)

write.csv(
  analysis_cells,
  file.path(
    tables_dir,
    paste0("NWS_full_grid_FFMSY_match_", run_label, ".csv")
  ),
  row.names = FALSE
)


# ============================================================
# 8. PREPARE ENVIRONMENTAL TIME SERIES
# ============================================================

environment <- readRDS(environment_file) %>%
  as_tibble()

required_environment_columns <- c(
  "lon",
  "lat",
  "year",
  "lzprod",
  "lzbio",
  "dfbot",
  "depth",
  "Tb",
  "Tm",
  "Tp"
)

missing_environment_columns <- setdiff(
  required_environment_columns,
  names(environment)
)

if (length(missing_environment_columns) > 0) {
  stop(
    "Missing environmental columns: ",
    paste(missing_environment_columns, collapse = ", ")
  )
}

cell_coordinates <- analysis_cells %>%
  select(cell_id, lon, lat)

environment_study <- environment %>%
  inner_join(cell_coordinates, by = c("lon", "lat")) %>%
  select(cell_id, all_of(required_environment_columns)) %>%
  arrange(cell_id, year)

rm(environment)
invisible(gc())

year_check <- environment_study %>%
  group_by(cell_id) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    number_of_years = n_distinct(year),
    duplicated_years = n() - n_distinct(year),
    .groups = "drop"
  )

if (
  any(year_check$first_year != 2015L) ||
  any(year_check$last_year != 2100L) ||
  any(year_check$number_of_years != 86L) ||
  any(year_check$duplicated_years != 0L)
) {
  stop("The environmental time series is incomplete for at least one cell.")
}

environment_study <- as.data.table(environment_study)
setkey(environment_study, cell_id, year)

cell_information <- as.data.table(analysis_cells)
setkey(cell_information, cell_id)

strategy_information <- as.data.table(best_strategy)
setkey(strategy_information, region, start_year)


# ============================================================
# 9. FEISTY FUNCTION FOR ONE CELL AND ONE PERIOD
# ============================================================

run_one_period <- function(
    environmental_time_series,
    cell_row,
    strategy_row
) {

  approach <- strategy_row$best_approach[1]
  fishing_value <- strategy_row$optimal_F_or_Fcoef[1]

  if (approach == "constant_F") {

    fishing_smallPel <- rep(
      fishing_value,
      nrow(environmental_time_series)
    )

    fishing_demersals <- fishing_smallPel
    fishing_largePel <- fishing_smallPel

  } else if (approach == "FFMSY_temp") {

    fishing_smallPel <-
      fishing_value *
      cell_row$ffmsy_smallPel[1] *
      exp(
        temperature_slope *
          (environmental_time_series$Tp - 10)
      )

    fishing_demersals <-
      fishing_value *
      cell_row$ffmsy_demersals[1] *
      exp(
        temperature_slope *
          (environmental_time_series$Tb - 10)
      )

    fishing_largePel <-
      fishing_value *
      cell_row$ffmsy_largePel[1] *
      exp(
        temperature_slope *
          (environmental_time_series$Tp - 10)
      )

  } else {
    stop("Unknown fishing approach: ", approach)
  }

  model_output <- NULL

  invisible(
    capture.output({

      parameters <- setupTimeseries(
        p = setupVertical2(
          depth = as.numeric(environmental_time_series$depth[1])
        ),
        tStep_ts = 1,
        tSpin = 40,
        szprod_ts = as.numeric(
          environmental_time_series$lzprod *
            factorProduction *
            factorSmallLarge
        ),
        szbio_ts = as.numeric(
          environmental_time_series$lzbio *
            factorSmallLarge *
            factorBiomass
        ),
        lzprod_ts = as.numeric(
          environmental_time_series$lzprod *
            factorProduction *
            (1 - factorSmallLarge)
        ),
        lzbio_ts = as.numeric(
          environmental_time_series$lzbio *
            (1 - factorSmallLarge) *
            factorBiomass
        ),
        dfbot = as.numeric(
          environmental_time_series$dfbot *
            factorBottomflux
        ),
        Tp = as.numeric(environmental_time_series$Tp),
        Tm = as.numeric(environmental_time_series$Tm),
        Tb = as.numeric(environmental_time_series$Tb),
        Fsmp_ts = as.numeric(fishing_smallPel),
        Fdem_ts = as.numeric(fishing_demersals),
        Flgp_ts = as.numeric(fishing_largePel)
      )

      model_output <- simulateFEISTY(p = parameters)
    })
  )

  model_year_rows <- 2:(nrow(environmental_time_series) + 1L)
  group_names <- parameters$groupnames[5:(4 + parameters$nGroups)]

  yield_matrix <- as.data.frame(
    model_output$yield[model_year_rows, , drop = FALSE]
  )

  names(yield_matrix) <- group_names
  yield_matrix$year <- environmental_time_series$year

  SSB_matrix <- as.data.frame(
    model_output$SSB[model_year_rows, , drop = FALSE]
  )

  names(SSB_matrix) <- group_names
  SSB_matrix$year <- environmental_time_series$year

  period_yield <- yield_matrix %>%
    filter(
      year >= strategy_row$start_year[1],
      year <= strategy_row$end_year[1]
    ) %>%
    summarise(
      smallPel = mean(smallPel, na.rm = TRUE),
      demersals = mean(demersals, na.rm = TRUE),
      largePel = mean(largePel, na.rm = TRUE)
    )

  period_SSB <- SSB_matrix %>%
    filter(
      year >= strategy_row$start_year[1],
      year <= strategy_row$end_year[1]
    ) %>%
    summarise(
      SSB_smallPel = mean(smallPel, na.rm = TRUE),
      SSB_demersals = mean(demersals, na.rm = TRUE),
      SSB_largePel = mean(largePel, na.rm = TRUE)
    )

  bind_cols(period_yield, period_SSB) %>%
    mutate(
      total = smallPel + demersals + largePel,
      climate_period = strategy_row$climate_period[1],
      decade = strategy_row$decade[1],
      approach = approach,
      F_or_Fcoef = fishing_value,
      .before = 1
    )
}


# ============================================================
# 10. UNFISHED SSB0 FOR ONE CELL
# ============================================================

run_unfished_baseline <- function(environmental_time_series) {

  model_output <- NULL

  invisible(
    capture.output({

      parameters <- setupTimeseries(
        p = setupVertical2(
          depth = as.numeric(environmental_time_series$depth[1])
        ),
        tStep_ts = 1,
        tSpin = 40,
        szprod_ts = as.numeric(
          environmental_time_series$lzprod *
            factorProduction *
            factorSmallLarge
        ),
        szbio_ts = as.numeric(
          environmental_time_series$lzbio *
            factorSmallLarge *
            factorBiomass
        ),
        lzprod_ts = as.numeric(
          environmental_time_series$lzprod *
            factorProduction *
            (1 - factorSmallLarge)
        ),
        lzbio_ts = as.numeric(
          environmental_time_series$lzbio *
            (1 - factorSmallLarge) *
            factorBiomass
        ),
        dfbot = as.numeric(
          environmental_time_series$dfbot *
            factorBottomflux
        ),
        Tp = as.numeric(environmental_time_series$Tp),
        Tm = as.numeric(environmental_time_series$Tm),
        Tb = as.numeric(environmental_time_series$Tb),
        Fsmp_ts = rep(0, nrow(environmental_time_series)),
        Fdem_ts = rep(0, nrow(environmental_time_series)),
        Flgp_ts = rep(0, nrow(environmental_time_series))
      )

      model_output <- simulateFEISTY(p = parameters)
    })
  )

  model_year_rows <- 2:(nrow(environmental_time_series) + 1L)
  group_names <- parameters$groupnames[5:(4 + parameters$nGroups)]

  SSB_matrix <- as.data.frame(
    model_output$SSB[model_year_rows, , drop = FALSE]
  )

  names(SSB_matrix) <- group_names
  SSB_matrix$year <- environmental_time_series$year

  map_dfr(
    seq_len(nrow(period_definition)),
    function(period_index) {

      period_row <- period_definition[period_index, ]

      SSB_matrix %>%
        filter(
          year >= period_row$start_year,
          year <= period_row$end_year
        ) %>%
        summarise(
          SSB0_smallPel = mean(smallPel, na.rm = TRUE),
          SSB0_demersals = mean(demersals, na.rm = TRUE),
          SSB0_largePel = mean(largePel, na.rm = TRUE)
        ) %>%
        mutate(
          climate_period = period_row$climate_period,
          .before = 1
        )
    }
  )
}


# ============================================================
# 11. FEISTY FUNCTION FOR ONE CELL
# ============================================================

run_one_cell <- function(cell_id_value) {

  tryCatch({

    environmental_time_series <- environment_study[
      .(cell_id_value)
    ]

    cell_row <- cell_information[
      .(cell_id_value)
    ]

    regional_strategies <- strategy_information[
      .(cell_row$sea[1])
    ]

    if (nrow(environmental_time_series) != 86L) {
      stop("Expected 86 environmental years.")
    }

    if (nrow(regional_strategies) != 2L) {
      stop("Expected two period strategies for the region.")
    }

    period_results <- lapply(
      seq_len(nrow(regional_strategies)),
      function(period_index) {
        run_one_period(
          environmental_time_series,
          cell_row,
          regional_strategies[period_index]
        )
      }
    )

    unfished_baseline <- run_unfished_baseline(
      environmental_time_series
    )

    bind_rows(period_results) %>%
      left_join(
        unfished_baseline,
        by = "climate_period"
      ) %>%
      mutate(
        cell_id = cell_id_value,
        sea = cell_row$sea[1],
        lon = cell_row$lon[1],
        lat = cell_row$lat[1],
        error = NA_character_,
        .before = 1
      )

  }, error = function(error_condition) {

    cell_row <- cell_information[
      .(cell_id_value)
    ]

    tibble(
      cell_id = cell_id_value,
      sea = cell_row$sea[1],
      lon = cell_row$lon[1],
      lat = cell_row$lat[1],
      error = conditionMessage(error_condition)
    )
  })
}


# ============================================================
# 12. CHUNKED PARALLEL RUN WITH RESUME SUPPORT
# ============================================================

strategy_signature <- paste(
  "validated_v2",
  run_label,
  nrow(analysis_cells),
  paste(
    best_strategy$region,
    best_strategy$decade,
    best_strategy$best_approach,
    format(best_strategy$optimal_F_or_Fcoef, scientific = FALSE),
    collapse = "|"
  ),
  sep = "__"
)

chunk_directory <- file.path(
  rds_dir,
  paste0(
    "NWS_full_grid_MSY_yield_change_VALIDATED_",
    run_label,
    "_chunks"
  )
)

dir.create(chunk_directory, recursive = TRUE, showWarnings = FALSE)

cell_chunks <- split(
  analysis_cells$cell_id,
  ceiling(seq_along(analysis_cells$cell_id) / chunk_size)
)

cat("\nChunks:", length(cell_chunks), "\n")

for (chunk_index in seq_along(cell_chunks)) {

  chunk_file <- file.path(
    chunk_directory,
    sprintf("chunk_%05d.RDS", chunk_index)
  )

  if (file.exists(chunk_file)) {

    existing_chunk <- readRDS(chunk_file)

    if (!identical(existing_chunk$signature, strategy_signature)) {
      stop(
        "Checkpoint signature mismatch in: ",
        chunk_file,
        "\nMove the old checkpoint directory before restarting."
      )
    }

    cat(
      "Chunk",
      chunk_index,
      "/",
      length(cell_chunks),
      "already complete.\n"
    )

    next
  }

  cat(
    "Running chunk",
    chunk_index,
    "/",
    length(cell_chunks),
    "- cells",
    min(cell_chunks[[chunk_index]]),
    "to",
    max(cell_chunks[[chunk_index]]),
    "\n"
  )

  if (.Platform$OS.type == "windows" || n_cores == 1L) {

    chunk_results <- lapply(
      cell_chunks[[chunk_index]],
      run_one_cell
    )

  } else {

    chunk_results <- parallel::mclapply(
      cell_chunks[[chunk_index]],
      run_one_cell,
      mc.cores = n_cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  }

  chunk_results <- bind_rows(chunk_results)

  temporary_chunk_file <- paste0(chunk_file, ".tmp")

  saveRDS(
    list(
      signature = strategy_signature,
      results = chunk_results
    ),
    temporary_chunk_file
  )

  if (!file.rename(temporary_chunk_file, chunk_file)) {
    stop("Could not finalize checkpoint: ", chunk_file)
  }
}


# ============================================================
# 13. COMBINE AND VERIFY RESULTS
# ============================================================

chunk_files <- file.path(
  chunk_directory,
  sprintf("chunk_%05d.RDS", seq_along(cell_chunks))
)

period_results <- map_dfr(
  chunk_files,
  function(chunk_file) {

    checkpoint <- readRDS(chunk_file)

    if (!identical(checkpoint$signature, strategy_signature)) {
      stop("Checkpoint signature mismatch in: ", chunk_file)
    }

    checkpoint$results
  }
)

error_results <- period_results %>%
  filter(!is.na(error))

if (nrow(error_results) > 0) {

  error_file <- file.path(
    tables_dir,
    paste0("NWS_full_grid_MSY_errors_", run_label, ".csv")
  )

  write.csv(error_results, error_file, row.names = FALSE)

  stop(
    nrow(error_results),
    " cell simulations failed. See: ",
    error_file
  )
}

expected_period_rows <- nrow(analysis_cells) * 2L

if (nrow(period_results) != expected_period_rows) {
  stop(
    "Expected ",
    expected_period_rows,
    " period rows but found ",
    nrow(period_results),
    "."
  )
}

if (anyDuplicated(period_results[c("cell_id", "climate_period")])) {
  stop("Duplicated cell-period results were produced.")
}

period_output_file <- file.path(
  rds_dir,
  paste0("NWS_full_grid_MSY_period_yields_", run_label, ".RDS")
)

saveRDS(period_results, period_output_file)

write.csv(
  period_results,
  file.path(
    tables_dir,
    paste0("NWS_full_grid_MSY_period_yields_", run_label, ".csv")
  ),
  row.names = FALSE
)


# ============================================================
# 14. RECALCULATE THE 30% SSB0 RULE ON THE FULL GRID
# ============================================================

regional_SSB_check <- period_results %>%
  group_by(sea, climate_period) %>%
  summarise(
    SSB_smallPel = mean(SSB_smallPel, na.rm = TRUE),
    SSB0_smallPel = mean(SSB0_smallPel, na.rm = TRUE),
    SSB_demersals = mean(SSB_demersals, na.rm = TRUE),
    SSB0_demersals = mean(SSB0_demersals, na.rm = TRUE),
    SSB_largePel = mean(SSB_largePel, na.rm = TRUE),
    SSB0_largePel = mean(SSB0_largePel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ratio_smallPel = SSB_smallPel / SSB0_smallPel,
    ratio_demersals = SSB_demersals / SSB0_demersals,
    ratio_largePel = SSB_largePel / SSB0_largePel,
    minimum_SSB_ratio_full_grid = pmin(
      ratio_smallPel,
      ratio_demersals,
      ratio_largePel,
      na.rm = TRUE
    ),
    limiting_group_full_grid = case_when(
      minimum_SSB_ratio_full_grid == ratio_smallPel ~ "smallPel",
      minimum_SSB_ratio_full_grid == ratio_demersals ~ "demersals",
      minimum_SSB_ratio_full_grid == ratio_largePel ~ "largePel",
      TRUE ~ NA_character_
    ),
    valid_30pct_full_grid =
      is.finite(minimum_SSB_ratio_full_grid) &
      minimum_SSB_ratio_full_grid >= SSB_threshold
  )

ssb_check_output_file <- file.path(
  tables_dir,
  if (full_run) {
    "NWS_FINAL_FULL_GRID_22317_cells_30pct_SSB_check.csv"
  } else {
    paste0("NWS_full_grid_30pct_SSB_check_", run_label, ".csv")
  }
)

write.csv(regional_SSB_check, ssb_check_output_file, row.names = FALSE)

cat("\n========================================\n")
cat("FULL-GRID 30% SSB0 CHECK\n")
cat("========================================\n\n")

print(
  as.data.frame(
    regional_SSB_check %>%
      select(
        sea,
        climate_period,
        ratio_smallPel,
        ratio_demersals,
        ratio_largePel,
        minimum_SSB_ratio_full_grid,
        limiting_group_full_grid,
        valid_30pct_full_grid
      )
  ),
  row.names = FALSE
)

if (full_run && any(!regional_SSB_check$valid_30pct_full_grid)) {
  stop(
    "At least one selected regional strategy fails the 30% SSB0 rule ",
    "when recalculated on the full spatial grid. The fishing value must ",
    "be reduced and the affected period rerun before producing maps."
  )
}


# ============================================================
# 15. PRESENT-FUTURE CHANGES
# ============================================================

reference_results <- period_results %>%
  filter(climate_period == "Reference 2020-2029") %>%
  transmute(
    cell_id,
    sea,
    lon,
    lat,
    reference_approach = approach,
    reference_F_or_Fcoef = F_or_Fcoef,
    reference_total = total,
    reference_smallPel = smallPel,
    reference_demersals = demersals,
    reference_largePel = largePel
  )

future_results <- period_results %>%
  filter(climate_period == "Future 2090-2100") %>%
  transmute(
    cell_id,
    future_approach = approach,
    future_F_or_Fcoef = F_or_Fcoef,
    future_total = total,
    future_smallPel = smallPel,
    future_demersals = demersals,
    future_largePel = largePel
  )

safe_percent_change <- function(reference_value, future_value) {
  if_else(
    is.finite(reference_value) &
      is.finite(future_value) &
      reference_value > 1e-10,
    100 * (future_value - reference_value) / reference_value,
    NA_real_
  )
}

yield_change <- reference_results %>%
  inner_join(future_results, by = "cell_id") %>%
  mutate(
    absolute_total = future_total - reference_total,
    absolute_smallPel = future_smallPel - reference_smallPel,
    absolute_demersals = future_demersals - reference_demersals,
    absolute_largePel = future_largePel - reference_largePel,

    percent_total = safe_percent_change(reference_total, future_total),
    percent_smallPel = safe_percent_change(
      reference_smallPel,
      future_smallPel
    ),
    percent_demersals = safe_percent_change(
      reference_demersals,
      future_demersals
    ),
    percent_largePel = safe_percent_change(
      reference_largePel,
      future_largePel
    )
  )

saveRDS(
  yield_change,
  file.path(
    rds_dir,
    paste0("NWS_full_grid_MSY_yield_change_", run_label, ".RDS")
  )
)

write.csv(
  yield_change,
  file.path(
    tables_dir,
    paste0("NWS_full_grid_MSY_yield_change_", run_label, ".csv")
  ),
  row.names = FALSE
)


# ============================================================
# 16. CHANGE SUMMARY
# ============================================================

change_long <- bind_rows(
  yield_change %>%
    transmute(
      cell_id,
      sea,
      lon,
      lat,
      group = "Total yield",
      absolute_change = absolute_total,
      percent_change = percent_total
    ),
  yield_change %>%
    transmute(
      cell_id,
      sea,
      lon,
      lat,
      group = "Small pelagics",
      absolute_change = absolute_smallPel,
      percent_change = percent_smallPel
    ),
  yield_change %>%
    transmute(
      cell_id,
      sea,
      lon,
      lat,
      group = "Demersals",
      absolute_change = absolute_demersals,
      percent_change = percent_demersals
    ),
  yield_change %>%
    transmute(
      cell_id,
      sea,
      lon,
      lat,
      group = "Large pelagics",
      absolute_change = absolute_largePel,
      percent_change = percent_largePel
    )
) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Total yield",
        "Small pelagics",
        "Demersals",
        "Large pelagics"
      )
    )
  )

change_summary <- change_long %>%
  group_by(sea, group) %>%
  summarise(
    number_of_cells = n(),
    mean_absolute_change = mean(absolute_change, na.rm = TRUE),
    median_absolute_change = median(absolute_change, na.rm = TRUE),
    mean_percent_change = mean(percent_change, na.rm = TRUE),
    median_percent_change = median(percent_change, na.rm = TRUE),
    percentage_cells_decreasing = 100 * mean(
      absolute_change < 0,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

summary_output_file <- file.path(
  tables_dir,
  if (full_run) {
    "NWS_FINAL_FULL_GRID_22317_cells_yield_change_summary.csv"
  } else {
    paste0("NWS_full_grid_MSY_change_summary_", run_label, ".csv")
  }
)

write.csv(change_summary, summary_output_file, row.names = FALSE)

cat("\n========================================\n")
cat("CHANGE SUMMARY\n")
cat("========================================\n\n")
print(as.data.frame(change_summary), row.names = FALSE)


# ============================================================
# 17. FULL-RUN MAPS
# ============================================================

if (full_run) {

  world <- map_data("world")

  grid_width <- median(diff(sort(unique(analysis_cells$lon))))
  grid_height <- median(diff(sort(unique(analysis_cells$lat))))

  group_file_names <- c(
    "Total yield" = "total_yield",
    "Small pelagics" = "smallPel",
    "Demersals" = "demersals",
    "Large pelagics" = "largePel"
  )

  create_change_map <- function(
      map_data_group,
      value_column,
      legend_title,
      plot_title,
      output_file
  ) {

    values <- map_data_group[[value_column]]

    robust_limit <- as.numeric(
      quantile(
        abs(values[is.finite(values)]),
        probs = 0.98,
        na.rm = TRUE
      )
    )

    if (!is.finite(robust_limit) || robust_limit <= 0) {
      robust_limit <- 1
    }

    plot_object <- ggplot() +
      geom_tile(
        data = map_data_group,
        aes(
          x = lon,
          y = lat,
          fill = .data[[value_column]]
        ),
        width = grid_width,
        height = grid_height
      ) +
      geom_polygon(
        data = world,
        aes(x = long, y = lat, group = group),
        fill = "grey85",
        colour = "grey40",
        linewidth = 0.2
      ) +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-robust_limit, robust_limit),
        oob = scales::squish,
        na.value = "grey70",
        name = legend_title
      ) +
      coord_quickmap(
        xlim = c(-11, 11),
        ylim = c(42, 63),
        expand = FALSE
      ) +
      labs(
        title = plot_title,
        subtitle = paste0(
          "2090-2100 minus 2020-2029 under the best valid regional strategy; ",
          "colour scale capped at the 98th percentile"
        ),
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_bw(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold")
      )

    ggsave(
      file.path(figures_dir, output_file),
      plot_object,
      width = 10,
      height = 9,
      dpi = 300
    )
  }

  for (group_name in levels(change_long$group)) {

    group_map_data <- change_long %>%
      filter(group == group_name)

    group_file_name <- unname(group_file_names[group_name])

    create_change_map(
      group_map_data,
      "percent_change",
      "Change (%)",
      paste0(group_name, " - percentage yield change"),
      paste0("NWS_FINAL_FULL_", group_file_name, "_percent.png")
    )

    create_change_map(
      group_map_data,
      "absolute_change",
      expression(Delta * " yield (g " * m^{-2} * " " * year^{-1} * ")"),
      paste0(group_name, " - absolute yield change"),
      paste0("NWS_FINAL_FULL_", group_file_name, "_absolute.png")
    )
  }
}


# ============================================================
# 18. FINAL MESSAGE
# ============================================================

cat("\n========================================\n")
cat("FULL-GRID YIELD-CHANGE WORKFLOW COMPLETED\n")
cat("========================================\n")
cat("\nRun mode:", run_label, "\n")
cat("Cells completed:", format(nrow(analysis_cells), big.mark = ","), "\n")
cat("Period results:", period_output_file, "\n")

if (!full_run) {
  cat(
    "\nTEST PASSED. Set full_run <- TRUE at the top of the script ",
    "to start the checkpointed full run.\n"
  )
}
