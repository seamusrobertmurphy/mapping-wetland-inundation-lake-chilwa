# GEE SAR-Landsat Pipeline: Operational Workplan

**Date:** 5 July 2026
**Script:** `05.scripts/gee_sar_landsat_pipeline.R`
**Purpose:** regenerate the Sentinel-1 and Landsat inundation time series end to end under your own Earth Engine credentials. I wrote the script; you run it, since I cannot authenticate your Earth Engine session.

## Prerequisites

You need R (4.x) with the packages `reticulate`, `sf`, and `jsonlite`, and a Python environment with the Earth Engine API (`earthengine-api`) installed and visible to reticulate. You need a Google Earth Engine account with access to a cloud project; the script defaults to `murphys-deforisk` but reads the `EE_PROJECT` environment variable if set. You must have completed a one-time `ee$Authenticate()` (or `earthengine authenticate` at the command line) under that account before the first run.

## One-time setup

Install the Earth Engine Python API into the environment reticulate will use, for example `pip install earthengine-api` in that environment. In R, confirm reticulate points at it: set `RETICULATE_PYTHON` to the interpreter path if the default is wrong. Authenticate once from R with `library(reticulate); ee <- import("ee"); ee$Authenticate()`, which opens a browser and stores a token. After that, `ee$Initialize(project = "your-project")` should succeed without prompting.

## Configuration

All settings sit in the `cfg` list at the top of the script: the Earth Engine project, an optional explicit Python path, the lake point that seeds the AOI, the HydroSHEDS basin level, the Sentinel-1 and Landsat year windows, the cloud-cover cap, and the output directory. Change these there rather than in the body. The defaults reproduce the manuscript's configuration.

## Run order and what happens

Run the whole script with `Rscript 05.scripts/gee_sar_landsat_pipeline.R` from the repository root, or step through it in RStudio. It initialises Earth Engine, sources the AOI from HydroSHEDS around the lake point, then in sequence builds the Sentinel-1 monthly water-fraction series, the Landsat annual index series, and the spectral-mixture annual fraction series, printing scene counts and the adaptive water threshold as it goes.

## Expected outputs

Three CSVs land in the configured output directory (`03.outputs/` by default): `s1_monthly_timeseries.csv` (date, year, month, scene count, mean VV, mean VH, water fraction), `landsat_annual_indices.csv` (year, scene count, the five basin-mean indices), and `sma_annual_fractions.csv` (year, open-water, flooded-vegetation, dry-vegetation, and soil fractions). These are the series the manuscript's Results figures are drawn from.

## Run-time and quota notes

The reductions run server-side on Earth Engine, but the script pulls one `getInfo()` per month and per year, so expect the Sentinel-1 loop (120 months) and the Landsat and SMA loops (41 years each) to take on the order of several minutes to tens of minutes depending on Earth Engine load. If you hit a `computation timed out` or `Too many concurrent aggregations` error, rerun; the loops are idempotent and simply overwrite the CSVs. For very large exports, the manuscript notes that selected composites can be exported to Drive or Assets for interferometric work in SNAP, which this script does not do.

## Checkpoints

After the AOI step, confirm the printed Sentinel-1 scene count is in the expected range (the manuscript reports 913 for the study window). After the Landsat step, confirm the harmonised scene count is plausible (the manuscript reports 1,335). If either is far off, check the AOI and the year windows before trusting the series. Spot-check the CSVs: Sentinel-1 water fraction should oscillate seasonally with wet-season peaks well above dry-season troughs, and the SMA open-water fraction should collapse in the 1995 and 2012 recession years.

## Known limitations carried from the notebook

The endmember reflectance values in the spectral-mixture step are reasonable defaults and should be replaced with locally sampled endmembers before final results. The adaptive Sentinel-1 threshold is calibrated on a single reference scene with a fixed fallback; review it against the printed value. Neither of these blocks a first run, but both should be revisited before the numbers are treated as final.
