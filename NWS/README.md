# FEISTY spatial yield analysis for the North West Shelf

This repository documents the spatial FEISTY analysis of yield change between 2020–2029 and 2090–2100 for four North West Shelf study regions: Bay of Biscay, Celtic Sea, English Channel, and North Sea.

The analysis compares total yield, small pelagics, demersals, and large pelagics under regionally selected fishing strategies. Every selected strategy is checked against the management constraint that spawning-stock biomass remains at or above 30% of unfished spawning-stock biomass (SSB0).

## Main files

- `NWS_FEISTY_spatial_yield_analysis.qmd`: reproducible Quarto document with live R code.
- `NWS FEISTY spatial yield analysis 2.0.pdf`: formatted report for direct reading.
- `NWS_FINAL_maps_FULL_22317_cells.R`: full FEISTY simulation, validation, summary, and mapping workflow.
- `outputs/tables/`: selected strategies and final tabular results.
- `outputs/figures/`: the eight final spatial yield-change maps.
- `outputs/rds/NWS_full_grid_four_regions_cells.RDS`: small spatial-cell classification used by the report and full workflow.
- `INSTRUCTIONS_AVANT_PUBLICATION.md`: short instructions in French for preparing and publishing the repository.

## Render the report

The repository already contains the small derived files and final outputs needed to render the report. The two large FEISTY inputs are not needed for this step.

1. Open `FEISTY_NWS.Rproj` in RStudio.
2. Ensure that Quarto and the R packages `tidyverse`, `data.table`, and `knitr` are available.
3. Open `NWS_FEISTY_spatial_yield_analysis.qmd` and select **Render**.

The equivalent terminal command is:

```sh
quarto render NWS_FEISTY_spatial_yield_analysis.qmd
```

## Rerun the full FEISTY simulation

The complete calculation requires the FEISTY package from the research environment and these two large local inputs:

```text
data/nws_for_feisty_proj_year.RDS
data/all_FFMSY.RDS
```

These files are deliberately excluded from GitHub. Their exact names and locations must not be changed.

The repository already includes the other two inputs required by the script:

```text
outputs/rds/NWS_full_grid_four_regions_cells.RDS
outputs/tables/NWS_decadal_optimisation_best_approach.csv
```

The script defaults to a one-cell-per-region test. Run that test first. Only after it succeeds should `full_run <- TRUE` be set in `NWS_FINAL_maps_FULL_22317_cells.R`. The full 22,317-cell calculation can take many hours and writes resumable checkpoints under `outputs/rds/`.

## Repository structure

```text
FEISTY_NWS_spatial_yield_analysis/
├── README.md
├── INSTRUCTIONS_AVANT_PUBLICATION.md
├── FEISTY_NWS.Rproj
├── _quarto.yml
├── NWS_FEISTY_spatial_yield_analysis.qmd
├── NWS FEISTY spatial yield analysis 2.0.pdf
├── NWS_FINAL_maps_FULL_22317_cells.R
├── data/
│   ├── README.md
│   ├── nws_for_feisty_proj_year.RDS       # local only; not on GitHub
│   └── all_FFMSY.RDS                      # local only; not on GitHub
└── outputs/
    ├── rds/
    │   └── NWS_full_grid_four_regions_cells.RDS
    ├── tables/
    │   ├── NWS_decadal_optimisation_best_approach.csv
    │   ├── NWS_FINAL_FULL_GRID_22317_cells_30pct_SSB_check.csv
    │   └── NWS_FINAL_FULL_GRID_22317_cells_yield_change_summary.csv
    └── figures/
        ├── NWS_FINAL_FULL_total_yield_percent.png
        ├── NWS_FINAL_FULL_total_yield_absolute.png
        ├── NWS_FINAL_FULL_smallPel_percent.png
        ├── NWS_FINAL_FULL_smallPel_absolute.png
        ├── NWS_FINAL_FULL_demersals_percent.png
        ├── NWS_FINAL_FULL_demersals_absolute.png
        ├── NWS_FINAL_FULL_largePel_percent.png
        └── NWS_FINAL_FULL_largePel_absolute.png
```

## Reproducibility boundary

The repository contains the analysis code, small spatial classification, selected strategies, final validation and summary tables, and final maps. It does not redistribute the two large FEISTY input datasets or the FEISTY package. Reproducing the published report is therefore possible from the repository alone; recomputing all model results additionally requires those local research inputs and the FEISTY software environment.
