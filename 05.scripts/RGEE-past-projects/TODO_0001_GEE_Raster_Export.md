# TODO 0001 — GEE Raster Export and Offline Map Conversion

**Date:** 2026-04-22
**Branch:** `todo/gee-raster-export` (feature branch off `main`, PR on completion)
**Status:** Approved 2026-04-22 — design decisions resolved in §6 below

---

## 1. Objective

Decouple the rendered book from live Google Earth Engine tile sessions. Every interactive map in the book currently stacks one or more `tmap::tm_tiles(ee_tile_url(...))` layers. Those tile URLs carry a session token that GEE invalidates within weeks, so rendered HTML pages go blank for the reader. Export each underlying GEE layer once as a Cloud-Optimised GeoTIFF, commit it under `data/rendered/`, and rewrite the map-building chunks to read from disk.

Prose, equations, tables, and static figures are out of scope. Only the GEE-dependent map widgets change.

## 2. Scope

### 2.1 Files in scope

| File | `ee_tile_url` calls | Interactive map widgets |
|---|---|---|
| `index.qmd` (preface) | 15 | 4 |
| `01-burn-area/index.qmd` | 12 | 5 |
| `02-fuel-stratification/index.qmd` | 18 | 7 |
| `03-fuel-consumption/index.qmd` | 3 | 2 |
| `04-emission-factors/index.qmd` | 3 | up to 2 |
| **Total** | **51** | **~20 widgets** |

The same underlying layer (e.g. the 25-year MCD64A1 burn composite, or the IGBP 17-class raster) is reused across chapters, so the number of unique rasters to export is ~12–15, not 51.

### 2.2 Unique rasters to export

| ID | Source | Derivation | Spatial resolution | Notes |
|----|--------|------------|--------------------|-------|
| `mcd64a1_burn_composite_2000_2025` | `MODIS/061/MCD64A1` | max over 2000-01 to 2025-12 | 500 m | binary burn/no-burn |
| `mcd64a1_burn_edges_2000_2025` | derived from above | morphological edge | 500 m | yellow edge overlay |
| `mcd64a1_annual_2001_2024` | `MODIS/061/MCD64A1` | one band per year | 500 m | multi-band stack |
| `mcd64a1_burn_frequency` | derived from annual stack | pixelwise sum | 500 m | integer 0–24 |
| `mcd64a1_filter_pre` | `MODIS/061/MCD64A1` | no uncertainty filter | 500 m | for filter-comparison map |
| `mcd64a1_filter_post` | `MODIS/061/MCD64A1` | <20 % uncertainty filter | 500 m | for filter-comparison map |
| `mod14a1_active_fires_2020` | `MODIS/061/MOD14A1` | confidence ≥ 30 | 1 km | active-fire/burned-area comparison |
| `mcd12q1_igbp_2020` | `MODIS/061/MCD12Q1` | LC_Type1 | 500 m | 17-class IGBP |
| `mcd12q1_fire_categories_2020` | derived from MCD12Q1 | reclassified to 3 IPCC fire categories | 500 m | |
| `forest_strata_honduras` | `WORLDCLIM/V1/BIO` bio01, bio12 + MCD12Q1 forest mask | tropical moist / dry / temperate / boreal | 500 m | |
| `ipcc_climate_zones_honduras` | `users/philipaudebert/.../ipcc_climate_1985-2015_corrigenda_raster` | clip to AOI | source native | |
| `ipcc_soils_honduras` | `projects/murphys-deforisk/assets/ipcc_soils` | clip to AOI | source native | |
| `landsat_dnbr_2020` | `LANDSAT/LC08/C02/T1_L2` pre/post fire pair | dNBR = (preNBR − postNBR) | 30 m | validation only |
| `hansen_forest_height_2020` | `UMD_hansen_global_forest_change_2024_v1_12` | treecover2000 thresholded | 30 m | validation |
| `esa_worldcover_2020` | `ESA_WorldCover_v100_2020` | LC class | 10 m | validation |

All rasters clipped to the Honduras GAUL level-0 boundary and reprojected to **EPSG:32616** (UTM zone 16N). Resolution preserved at source where possible; 10 m ESA WorldCover resampled to 30 m to keep file size manageable.

### 2.3 Out of scope

- Any change to equations, tables, or narrative prose.
- Any change to static ggplot figures or Reactable tables already cached in `_freeze/`.
- Any new analytical content or methodological change.
- FAOSTAT data, already local in `data/raw/honduras_faostat_fires.csv`.

## 3. IPCC section referenced

No IPCC equation is implemented in this task. The exports preserve the inputs to Eq. 2.27 (Vol. 4 Ch. 2). No parameter values change.

## 4. Inputs required

- Authenticated GEE account with read access to:
  - `MODIS/061/*` (public)
  - `FAO/GAUL/2015/*` (public)
  - `WORLDCLIM/V1/BIO` (public)
  - `UMD_hansen_global_forest_change_2024_v1_12` (public)
  - `ESA_WorldCover_v100_2020` (public)
  - `users/philipaudebert/IPCC/ipcc_climate_1985-2015_corrigenda_raster` (public-shared)
  - `projects/murphys-deforisk/assets/ipcc_soils` (requires user's own access)
- Python `earthengine-api` in the MacPorts Python at `/opt/local/bin/python`.
- `gcloud` authentication for large exports to Google Cloud Storage; or direct `ee.batch.Export.image.toDrive` if Drive is acceptable.
- R packages already in the pipeline: `terra`, `sf`, `tmap`, `leaflet`, `mapview`.

## 5. Expected outputs

### 5.1 Directory structure

```
data/
├── raw/                           (unchanged — immutable source evidence)
├── rendered/                      (NEW — committed derived rasters)
│   ├── metadata/
│   │   └── <layer_id>_<YYYYMMDD>.json
│   ├── mcd64a1_burn_composite_2000_2025.tif
│   ├── mcd64a1_burn_edges_2000_2025.tif
│   ├── mcd64a1_annual_2001_2024.tif
│   ├── ...
│   └── esa_worldcover_2020_hnd.tif
└── raw/metadata/                  (unchanged)
```

Each `<layer_id>_<YYYYMMDD>.json` records: source collection, asset ID, date range, region (GeoJSON bbox), CRS, resolution, export timestamp, and the GEE reducer or visualisation expression used. Per `CLAUDE.md` §5 Step 1.

### 5.2 Code artefacts

- `gee/export_rasters.py` — single Python script that iterates the layer list in §2.2 and writes each COG + its metadata JSON. Idempotent: skip an output that already exists unless `--force` is passed.
- `gee/layer_manifest.yaml` — machine-readable registry of the 15 exports, consumed by the Python script and by the R chunks when loading.
- `R/local_rasters.R` — small helper that maps a layer ID to its on-disk path and returns a `terra::SpatRaster`, optionally re-projected for display.

### 5.3 Rewritten chunks

Each of the ~20 interactive map chunks is rewritten to:

1. Load the relevant COG(s) via `terra::rast(here::here("data/rendered/<id>.tif"))`.
2. Build the map with `leaflet::leaflet() |> addRasterImage(...)` or equivalent, matching the original colour palette and legend.
3. Wrap the original GEE extraction chunk with `#| eval: false` so the methodology stays visible in the source but does not execute on render.

### 5.4 Rendered artefacts

- `_book/` HTML re-rendered end to end with no GEE session active, confirming independence.
- PDF target also renders (per `CLAUDE.md` §3 Quarto section): static `tmap` fallback used in PDF, interactive `leaflet` in HTML, conditional on `knitr::is_latex_output()`.

## 6. Design decisions (resolved 2026-04-22)

1. **Two-tier persistence: GEE asset + local COG.** Source imagery is materialised as a persistent asset under `projects/<user>/assets/trees-wildfire-replications/<layer_id>` using `ee.batch.Export.image.toAsset`. From that asset, a COG is downloaded to `data/rendered/` for the book to load on render. The asset is the authoritative backup; the COG is the render-time consumer. This matters because GEE tile URLs — even those generated from user-owned assets — carry an expiring session token, so loading from tile URL alone would not fix the expiry problem.

2. **Download method.** No Google Drive. `rgee::ee_as_raster(via = "getDownloadURL")` is the primary path (synchronous, hits the `computePixels` endpoint, capped at 32 MB uncompressed per request). The 32 MB ceiling covers every 500 m and 1 km layer over the Honduras AOI comfortably. Validation rasters at native 30 m (Landsat dNBR, Hansen) and 10 m (ESA WorldCover) exceed this. For those, the asset is exported at native resolution (canonical copy), then a **display COG** is downloaded resampled to 100 m for book rendering — this is sufficient for the qualitative visual checks they support and avoids both Drive and release attachments.

3. **Size budget and overflow to GitHub release.** Expected compressed COG sizes (Honduras, EPSG:32616):
   - 500 m single-band: <1 MB
   - 500 m 24-band annual stack: ~15 MB
   - 1 km MOD14A1: <0.5 MB
   - 100 m resampled Landsat dNBR: ~4 MB
   - 100 m resampled Hansen: ~2 MB
   - 100 m resampled ESA WorldCover: ~2 MB

   No layer is expected to exceed 50 MB after compression. **If any layer does exceed 50 MB**, it moves to a GitHub release on `seamusrobertmurphy/TREES-wildfire-replications` tagged `gee-exports-v1`, fetched at render time with `piggyback::pb_download()` into a local cache directory (gitignored). The book loads from the downloaded local file the same way it loads committed COGs — the render chunk is indifferent to which. I will report any overflow cases to you before uploading to the release.

4. **IPCC soils asset ownership.** Confirmed: `projects/murphys-deforisk/assets/ipcc_soils` is your own asset. Export will run under your authenticated GEE account. No substitution needed.

5. **Peatland screening gate.** `CLAUDE.md` §6.2 requires peatland screening before SOC calculation. This TODO touches no SOC work, so the gate does not apply. The FEATURE file will note that the IPCC soils export does **not** imply peatland clearance for downstream SOC work.

6. **AOI definition.** `FAO/GAUL/2015/level0` filtered by `ADM0_NAME == "Honduras"` remains authoritative. Exported once to `data/rendered/aoi_honduras.gpkg` and reused by every raster download and every render chunk.

7. **Colour palettes.** Preserve every existing palette and breakpoint exactly. Stored in `R/palettes.R` so the render-time map builder uses the same values as the current `ee_tile_url(..., visParams)` calls.

8. **Working tree.** User will rebase/tidy `main` independently. I branch off current `main` for the work but will not touch the 53 uncommitted working-tree changes. My branch commits only the new files and the edited `.qmd` chunks.

## 7. Execution plan (phased)

**Phase A — Infrastructure (no GEE calls yet).**
- Create `todo/gee-raster-export` branch off `main`.
- Create `data/rendered/`, `data/rendered/metadata/`, `gee/`, `R/palettes.R`, `R/local_rasters.R`.
- Write `gee/layer_manifest.yaml` with the 15 entries.
- Write `gee/export_to_asset.py` (source → `projects/<user>/assets/trees-wildfire-replications/<id>`) and `gee/download_cogs.R` (asset → `data/rendered/<id>.tif`).
- Create `INDEX.md` at repo root (currently missing, spec requires it).
- Commit infrastructure. Pause for user review.

**Phase B — Export to assets.**
- Authenticate GEE as user.
- Run `gee/export_to_asset.py` to submit one `ee.batch.Export.image.toAsset` task per layer. Monitor with `ee.batch.Task.list()` until all complete.
- Commit the metadata JSON for each asset (timestamp, source, reducer expression).

**Phase C — Download to COG.**
- Run `gee/download_cogs.R` against the now-persistent assets.
- For each layer: if compressed size ≤50 MB, commit to `data/rendered/`; if >50 MB, upload to the GitHub release `gee-exports-v1` and record the piggyback path in `gee/layer_manifest.yaml`.
- Report actual sizes back to user before any release upload.

**Phase D — Chunk rewrites.**
- Rewrite preface (`index.qmd`) map chunks first as proof of concept.
- Render HTML + PDF and visually inspect each map matches the previous live version.
- On confirmation, repeat for chapters 1–4 in order.

**Phase E — Verification and handover.**
- Render the full book with the GEE Python kernel disabled (set `reticulate::use_python` to a venv without `earthengine-api`) to prove independence.
- Write `FEATURE_0001_GEE_Raster_Export.md`.
- Update `INDEX.md`.
- Open PR for user review.

## 8. Approval checkpoints

- **Phase A → B:** user reviews infrastructure commits before asset exports begin.
- **Phase C:** user reviews actual COG sizes before any GitHub release upload.
- **Phase D:** user inspects preface renders before chapters 1–4 are rewritten.
- **Phase E:** user reviews PR.

## 9. Known open risks

- **GEE asset project path.** User's default GEE project is required. I will inspect the user's `~/.config/earthengine/credentials` and `ee.data.getProjectConfig()` to determine the correct prefix. If the user's default project does not permit asset writes, the export halts and I ask before proceeding.
- **Palette fidelity in leaflet vs. tmap.** Leaflet's `colorNumeric` and `colorFactor` do not always match GEE's `visParams$palette` exactly for categorical rasters. A side-by-side test chunk renders each map (GEE live vs. local) during Phase D to catch drift before the book-wide rewrite.
- **Asset export concurrency.** GEE limits concurrent export tasks. If the queue grows beyond 5–10 tasks, I will submit in batches rather than all at once, and will not advance to Phase C until every task reports `COMPLETED`.
- **Uncommitted working-tree changes.** User will tidy `main` independently. I branch from current `main` and commit only new files and edited chunks.

---

*End of TODO 0001 — beginning Phase A.*
