# FEISTY spatial yield analyses

This repository contains two related FEISTY analyses. They are kept in
separate folders so that the regional North West Shelf work and the global
simulations can be understood and rerun independently.

## Analyses

### [NWS analysis](NWS/)

Spatial yield analysis for the Bay of Biscay, Celtic Sea, English Channel and
North Sea. The folder contains the Quarto report, the full R workflow, selected
tables and the final maps.

The main R script is:

```text
NWS/NWS_FINAL_maps_FULL_22317_cells.R
```

### [Global analysis](Global/)

Global FEISTY simulations on the one-degree ocean grid from 2015 to 2087. The
folder contains the complete sequence of R scripts, from the first smoke test
to the full 41,008-cell production run, global maps and the LME comparison.

All global R code is in:

```text
Global/R/
```

The scripts are numbered in their normal execution order. Each script starts
with a description of its objective, inputs, outputs and main calculation
steps. Comments beginning with `#` appear in green in RStudio.

## Repository structure

```text
.
├── README.md
├── NWS/
│   ├── README.md
│   ├── NWS_FEISTY_spatial_yield_analysis.qmd
│   ├── NWS_FINAL_maps_FULL_22317_cells.R
│   ├── data/
│   └── outputs/
└── Global/
    ├── README.md
    ├── FEISTY_Global.Rproj
    ├── R/
    │   ├── 01_global_climate_baseline_smoke_test.R
    │   ├── 02_global_fixed_F_smoke_test.R
    │   ├── 03_global_constant_F_scan_smoke_test.R
    │   ├── 04_global_pilot_run.R
    │   ├── 05_global_production_5deg.R
    │   ├── 06_global_full_grid_production.R
    │   ├── 07_global_full_grid_maps.R
    │   └── 08_global_LME_observed_vs_simulated_yield.R
    ├── data/
    └── outputs/
```

## Data availability

The code, selected derived tables and final figures are included. Large raw
forcing files, full checkpoints and other intermediate model objects are not
stored on GitHub. The README inside each analysis folder lists the exact local
input filenames and explains how to supply their paths.

The global climate forcing currently available ends in December 2087, so the
global projection cannot be extended to 2100 without a new forcing dataset.
