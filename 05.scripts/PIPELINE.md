# The Earth Engine pipeline, preserved

Written 4 September 2026. Every script that assembled, resampled, processed,
evaluated and normalised the multi-source raster time series, in the order they
run, with what each one decides and what it leaves behind. This is the record of
the Earth Engine stage, kept so the analysis can be rerun or audited without
reconstructing it from memory.

Everything runs against `python3` (3.12, MacPorts) with `earthengine-api`,
`rasterio`, `geopandas`, `scipy`, `numpy` and `matplotlib`, authenticated as
`seamusrobertmurphy@gmail.com` under the Earth Engine project
`murphys-deforisk`. No virtual environment is needed and none should be created.

## Order of execution

### 1. Establish what exists

| Script | Writes | Decides |
|:---|:---|:---|
| `survey_optical_archive.py` | `optical_archive_survey.csv` | Which sensors see this basin and when, across eleven collections at three cloud screens. Settled that 1988 is not empty, that no Multispectral Scanner reaches here, and that nothing precedes 1984. |

Run this first and before any claim of absence. Filenames are not evidence of
coverage and a catalogue is not a footprint.

### 2. Extract observations

| Script | Writes | Decides |
|:---|:---|:---|
| `extract_optical_water_16day.py` | `optical_water_16day.csv` | Water area actually seen in each sixteen-day window, per instrument family, with the clear fraction behind it. Steps with no clear observation return nothing. |
| `extract_avhrr_index.py` | `avhrr_index_16day.csv` | Kept for the record. Advanced Very High Resolution Radiometer was tested as the pre-2000 dense stream and rejected at a best correlation of 0.33. Do not reinstate without a new idea. |
| `asar_process.py` | `03.outputs/TIF/asar/*_sigma0.tif` | Envisat products to calibrated backscatter on the optical grid, without SNAP. GDAL reads the format directly. |
| `asar_series_join.py` | `asar_series_join.csv` | The radar joined to the optical model date by date. Seven of sixteen dates have a swath covering the lake. |
| `palsar_flooded_vegetation.py` | `palsar_flooded_vegetation.csv` | L-band double bounce for 2007 to 2010 and 2015 to 2024. |
| `recover_split_members.py` | extracted `.N1` products | Recovers members from a concatenated split zip whose local headers sit past the first slice, which `unzip` cannot read. |

### 3. Fuse and normalise

| Script | Writes | Decides |
|:---|:---|:---|
| `build_inundation_16day.py` | `chilwa_inundation_16day.csv`, its validation and provenance | The record itself. Calibrates the coarse stream onto the fine one, sets the clear-cover floor, and runs a local linear trend with an annual cycle through a Kalman filter and a backward smoothing pass, fitting both the step size and the observation error scale by maximum likelihood. |

Three constants in that script are measurements, not choices, and each has its
reason in the file. The clear-cover floor of eighty per cent, because below it
the scale-up error runs to hundreds of square kilometres. The Landsat retrieval
error of forty square kilometres, refined by the fitted scale of 0.377. And the
radar ratio of 0.11 that carries the vegetated column between radar dates.

### 4. Evaluate

| Script | Writes | Decides |
|:---|:---|:---|
| `eval_timeseries.py` | `eval_timeseries.json` | Whether the record is finished. Sixteen criteria, exit zero only when all pass. This file is the contract and nothing else decides completeness. |
| `diagnose_cube_consistency.py` | `cube_consistency_tests.csv` | Whether the fusion is internally consistent: innovation whiteness, interval calibration, per-sensor pull, handover steps, and leave-one-sensor-out. |
| `eda_cube_statistics.py` | `eda_table1` to `table3`, and the markdown twins | Inventory, distributions and bias against the Landsat reference, and what each era can support. |

`eval_timeseries_annual.py` is the superseded annual contract, kept only so the
decision to abandon an annual product can be re-read.

### 5. Products and aids

| Script | Writes | Purpose |
|:---|:---|:---|
| `plot_16day_record.py` | `fig08_inundation_16day.png` | The record drawn, with its credible band and the observation count beneath it. |
| `export_basemap_s2.py` | `basemap/chilwa_s2_10m_*.tif` | Sentinel-2 at ten metres, two seasons, nine years, gap free. Stretch measured from the image. |
| `export_basemap_esri.py` | `basemap/chilwa_esri_2.3m_z16.tif` | Sub-five-metre positioning aid. No calibrated radiometry; never enters a measurement. |
| `export_ground_truth_aid.py` | `ground_truth_aid/*.tif` | Per-date imagery for relocating photographs, in a near-date and a gap-filled version. |
| `export_sits_cube_sample.py` | `02.inputs/CUBE/chilwa_2012_30m/` | The handover. Pixels on a sixteen-day grid, named for the sits local-cube reader. |
| `export_harmonic_coefficients.py` | `harmonic_2012/harmonic_coefficients.tif` | The fitted optical surface, so the model can be read on any date offline. |
| `harmonic_fit_2012.py` | `mndwi_harmonic_2012-02-22.tif` | The settled 2012 water map, validated against Water Observations from Space. |

### 6. Not for reuse

`harmonic_open_water_series.py` and `sma_harmonic_flooded_veg.py` built the
annual product that was discarded. They remain because the adaptive window logic
in the first is sound and would be needed again if a per-pixel record is ever
built, but neither feeds anything current.

## Two Earth Engine limits every script here works around

A request carrying more than roughly sixteen concurrent aggregations is refused,
so reductions are folded into single multi-band calls and steps are sent in
chunks of eight.

An interactive download recomputes its whole expression for every tile, which
exceeds the user memory limit on anything built from hundreds of scenes. The
pattern that works is to export once to an asset and then pull tiles from the
stored pixels, where no computation happens. `export_basemap_s2.py` and
`export_harmonic_coefficients.py` both show it.

## What is deliberately not here

Advanced Very High Resolution Radiometer, tested and rejected. Sentinel-2 in the
forty-one-year record, held out so a Landsat series does not gain a sensor step
where it is densest. Planet NICFI and HISTARFM, neither of which covers this
basin. The evidence for each is in
`03.outputs/TABLES/table4_candidate_datasets.md`.
