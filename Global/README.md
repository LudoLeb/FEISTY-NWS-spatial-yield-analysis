# FEISTY global simulations

This folder contains the global continuation of the North West Shelf FEISTY
analysis. **All R code is in the `R/` folder.** The scripts are numbered in the
order in which they should normally be used, and the main steps are documented
with comments that appear in green in RStudio.

The main report is `Global_FEISTY_spatial_analysis.qmd`. It describes the full
workflow and includes the selection, validation, time-series, spatial and LME
comparison results. The rendered PDF can be read without running the code.

## R code workflow

| Script | Purpose |
|---|---|
| `R/01_global_climate_baseline_smoke_test.R` | Checks the forcing and the no-fishing setup on 12 cells. |
| `R/02_global_fixed_F_smoke_test.R` | Adds an initial constant-F sensitivity on the same cells. |
| `R/03_global_constant_F_scan_smoke_test.R` | Tests F = 0.05 to 0.60 and selects F = 0.40 under the 30% SSB0 rule. |
| `R/04_global_pilot_run.R` | Validates F = 0.40 on 134 stratified cells. |
| `R/05_global_production_5deg.R` | Runs an intermediate representative five-degree grid. |
| `R/06_global_full_grid_production.R` | Runs F = 0 and F = 0.40 on all 41,008 one-degree cells. |
| `R/07_global_full_grid_maps.R` | Creates the four global SSB and yield map figures. |
| `R/08_global_LME_observed_vs_simulated_yield.R` | Compares observed and simulated yield among LMEs. |
| `R/monitor_progress.R` | Displays progress for the five-degree run. |
| `R/monitor_full_grid_progress.R` | Displays progress for the full-grid run. |

## Analysis scope

- Monthly global forcing from January 2015 to December 2087.
- FEISTY global correction factors from the established calibration workflow.
- Full one-degree ocean grid containing 41,008 cells.
- Two scenarios: climate-only `F = 0` and constant fishing `F = 0.40 yr-1`.
- Annual SSB and yield for small pelagics, large pelagics and demersals.
- Five-year means for 2015-2019 and 2083-2087.

The same fishing pressure is applied across time, cells and commercial groups.
It is a controlled experiment, not an estimate of observed fishing pressure in
each region of the world.

## Smoke-test result

A scan from `F = 0.05` to `0.60 yr-1` was evaluated against the 30% SSB0
constraint using final-period SSB aggregated over the twelve test cells. The
largest tested value that kept all three commercial groups above the threshold
was `F = 0.40 yr-1`. Large pelagics were the limiting group: their aggregate
SSB was 31.7% of the no-fishing reference at `F = 0.40`, but 28.4% at
`F = 0.45`. Therefore `F = 0.40` is the provisional fishing scenario for the
larger pilot run; it is not yet a final global management recommendation.

## Stratified pilot result

The larger pilot contained 134 spatially stratified cells: 64 shelf cells and
70 open-ocean cells. Both `F = 0` and `F = 0.40` ran successfully from 2015 to
2087, producing 268 completed simulations with no missing output values.

In the 2083-2087 mean, weighted aggregate SSB under `F = 0.40` represented:

- 89.4% of SSB0 for small pelagics;
- 33.8% for large pelagics;
- 60.2% for demersals;
- 56.6% for the three-group total.

The pilot therefore confirms `F = 0.40` as a suitable candidate for the full
global run under the aggregate 30% SSB0 rule. Large pelagics remain the
limiting group, so the full-grid constraint must be checked before treating
this value as final.

## Full-grid result

The production run covers all 41,008 ocean cells on the available 1-degree
grid. Both `F = 0` and `F = 0.40 yr-1` were simulated from January 2015 to
December 2087. The completed output contains all expected cell, scenario,
metric and group combinations, with no missing or duplicated result keys.

For the 2083-2087 mean, aggregate SSB under `F = 0.40` represented 89.7% of
SSB0 for small pelagics, 36.8% for large pelagics, 55.5% for demersals and
57.6% for the three-group total. All groups therefore remain above the 30%
SSB0 threshold in the full-grid simulation.

Four comparable global map figures are produced with
`R/07_global_full_grid_maps.R`:

- SSB in 2015-2019, total and by functional group;
- yield in 2015-2019, total and by functional group;
- absolute SSB change between 2015-2019 and 2083-2087;
- absolute yield change between 2015-2019 and 2083-2087.

Each four-panel figure uses one common colour range so values can be compared
directly among the total, small pelagic, demersal and large pelagic panels.
The Word file `Global_FEISTY_spatial_results.docx` collects the four figures
on separate landscape pages for discussion and editing.

## Input data

The raw file `global_feisty_projection.RDS` is approximately 1 GB and is not
copied into this folder or intended for GitHub. Its local path can be supplied
through the `FEISTY_GLOBAL_FORCING` environment variable.

The LME comparison additionally requires these two small research inputs:

```text
data/raw/grid_LME_1deg.RData
data/raw/catch_Ftype_LME_year.RData
```

They are not redistributed in this repository. Alternative local paths can be
provided through `FEISTY_LME_GRID` and `FEISTY_LME_CATCH`.

## Run the complete workflow

```bash
FEISTY_GLOBAL_FORCING="/path/to/global_feisty_projection.RDS" \
  Rscript R/01_global_climate_baseline_smoke_test.R

Rscript R/02_global_fixed_F_smoke_test.R
Rscript R/03_global_constant_F_scan_smoke_test.R

FEISTY_GLOBAL_FORCING="/path/to/global_feisty_projection.RDS" \
  Rscript R/04_global_pilot_run.R

FEISTY_GLOBAL_FORCING="/path/to/global_feisty_projection.RDS" \
  Rscript R/05_global_production_5deg.R

FEISTY_GLOBAL_FORCING="/path/to/global_feisty_projection.RDS" \
  Rscript R/06_global_full_grid_production.R

Rscript R/07_global_full_grid_maps.R

FEISTY_LME_GRID="/path/to/grid_LME_1deg.RData" \
FEISTY_LME_CATCH="/path/to/catch_Ftype_LME_year.RData" \
  Rscript R/08_global_LME_observed_vs_simulated_yield.R
```

To render the report from the selected tables and figures already included in
the repository:

```bash
quarto render Global_FEISTY_spatial_analysis.qmd
```

The long production runs use checkpoints under `outputs/checkpoints/`, so an
interrupted calculation can be resumed. The selected final tables and figures
are kept in the repository; raw forcing, large intermediate objects,
checkpoints and logs are excluded.

## LME comparison

The current LME figure compares mean observed catch for 1995-2004 with mean
simulated yield for 2015-2019 under F = 0.40. Simulated cell yield is converted
from g wet weight m-2 yr-1 to tonnes yr-1 using the cell area and is then summed
within each LME. The observed total contains only the three groups represented
by FEISTY. The periods do not overlap, so this is a broad spatial comparison
rather than a strict same-year model validation.

The small derived table
`outputs/tables/global_full_grid_yield_2015_2019_for_LME.csv` contains only the
four yield groups and coordinates needed by script 08. It is included so the
LME comparison can be rerun without downloading the complete 88 MB period
summary.
