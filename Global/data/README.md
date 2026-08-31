# Data

The raw global projection is not included because of its size. The required
input is:

- `global_feisty_projection.RDS`: monthly SSP3-7.0 forcing with `lon`, `lat`,
  `time`, `bprodin`, `zbio`, `zprod`, `Tm`, `Tp`, `Tb` and `depth`.

In the file currently available for this project, the forcing covers January
2015 to December 2087. The supervisor confirmed that no later version is
available.

`derived/global_smoke_test_forcing_12_cells.rds` is generated from the raw
file by the smoke-test script and contains only twelve representative cells:
six shelf locations and six open-ocean locations.

The LME comparison in `R/08_global_LME_observed_vs_simulated_yield.R` also
expects these two local files:

- `raw/grid_LME_1deg.RData`: one-degree cells with LME number, name and area;
- `raw/catch_Ftype_LME_year.RData`: observed catch by LME, functional type and
  year from 1961 to 2004.

These research inputs are not committed to GitHub. They can be stored under
`data/raw/`, or their locations can be supplied with the environment variables
`FEISTY_LME_GRID` and `FEISTY_LME_CATCH`.
