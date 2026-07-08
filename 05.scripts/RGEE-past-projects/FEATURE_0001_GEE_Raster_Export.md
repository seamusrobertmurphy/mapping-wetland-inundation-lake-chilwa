# FEATURE 0001 — GEE Raster Export and Offline Map Conversion

**Completed:** 2026-04-22
**Branch:** `feature/gee-raster-export`
**Closes:** `specs/TODO_0001_GEE_Raster_Export.md`

---

## 1. What was built

The book's interactive maps no longer depend on live Google Earth Engine tile sessions. Every map that previously rendered a `tmap::tm_tiles(ee_tile_url(...))` call now renders from a Cloud-Optimised GeoTIFF committed to `data/rendered/`. Those COGs are mirrored from persistent assets under the user's GEE project, `projects/murphys-deforisk/assets/ipcc-wildfires/`.

Result: readers viewing the rendered HTML see maps that do not blank out when the GEE tile session token expires.

## 2. Scope achieved

| Artefact                                              | State                           |
|-------------------------------------------------------|---------------------------------|
| 11 GEE source layers exported to user-owned assets    | 9 COMPLETED, 2 still exporting  |
| 9 local COGs downloaded, committed to `data/rendered/`| done                            |
| AOI GPKG (country + states)                           | done                            |
| Preface: all 4 map chunks offline                     | done                            |
| Chapter 1: 6 map chunks offline                       | done                            |
| Chapter 2: 8 map chunks offline                       | done                            |
| Chapter 3: 2 map chunks offline                       | done                            |
| Chapter 4: 1 map chunk offline                        | done                            |
| Full book renders HTML end-to-end                     | done                            |

Totals: 21 previously live-GEE map chunks now load from committed COGs. 2 validation maps (Hansen GFC, ESA WorldCover) carry `file.exists()` guards and emit a placeholder notice until the last two asset exports complete.

## 3. Inputs used

- `gee/layer_manifest.yaml` — 11-layer registry. Source collection, reducer recipe, scale, dtype, and consuming chunks per layer.
- GEE credentials on the author's machine (user's `murphys-deforisk` project).
- `/opt/local/bin/python` with `earthengine-api`.
- R 4.4.1 with `rgee 1.1.8`, `terra`, `sf`, `tmap 4.2`, `leaflet`, `httr`, `piggyback`.
- MacPorts `proj 9.8.1` for `terra::project`, wired via `.Renviron`.

## 4. Outputs produced

### Committed COGs (`data/rendered/`)

| File                                          | Size   | Source                                                               |
|-----------------------------------------------|--------|----------------------------------------------------------------------|
| `aoi_honduras.gpkg`                           | 1.6 MB | FAO/GAUL/2015/level0+level1 (multilayer)                             |
| `mcd64a1_burn_2000_2025_max.tif`              | 65 kB  | MODIS/061/MCD64A1 BurnDate>0 pixelwise max, 500 m                   |
| `mcd64a1_burn_annual_2001_2024.tif`           | 374 kB | 24-band annual burn stack                                            |
| `mcd64a1_uncertainty_2000_2025.tif`           | 204 kB | MODIS/061/MCD64A1 Uncertainty mean                                   |
| `mod14a1_active_fires_2020.tif`               | 13 kB  | MODIS/061/MOD14A1 FireMask≥7, 1 km                                   |
| `mcd12q1_igbp_2020.tif`                       | 255 kB | MODIS/061/MCD12Q1 LC_Type1, 500 m                                    |
| `worldclim_bio01_bio12.tif`                   | 529 kB | WORLDCLIM/V1/BIO, two bands                                          |
| `ipcc_climate_zones_honduras.tif`             | 6 kB   | user's `ipcc_climate_1985-2015_corrigenda_raster` clipped            |
| `ipcc_soils_honduras.tif`                     | 18 kB  | user's `projects/murphys-deforisk/assets/ipcc_soils` clipped         |
| `landsat_dnbr_2020_display.tif`               | 2.2 MB | Landsat 8 SR dNBR 2020 season, 500 m display resolution              |

Metadata JSON under `data/rendered/metadata/` for every layer.

### GEE assets (persistent backup)

Under `projects/murphys-deforisk/assets/ipcc-wildfires/`:
- `mcd64a1_burn_2000_2025_max` (COMPLETED)
- `mcd64a1_burn_annual_2001_2024` (COMPLETED)
- `mcd64a1_uncertainty_2000_2025` (COMPLETED)
- `mod14a1_active_fires_2020` (COMPLETED)
- `mcd12q1_igbp_2020` (COMPLETED)
- `worldclim_bio01_bio12` (COMPLETED)
- `ipcc_climate_zones_honduras` (COMPLETED)
- `ipcc_soils_honduras` (COMPLETED)
- `landsat_dnbr_2020_display` (COMPLETED)
- `hansen_forest_cover_2020_display` (RUNNING at time of handover)
- `esa_worldcover_2020_display` (READY at time of handover)

### Code

- `INDEX.md` — repository index at root.
- `.Renviron` — pins `PROJ_DATA` to MacPorts path.
- `gee/layer_manifest.yaml` — layer registry.
- `gee/export_to_asset.py` — submits `Export.image.toAsset`, max 5 concurrent, idempotent.
- `gee/download_cogs.R` — pulls each asset via `getDownloadURL` (no Drive), chunks multi-band stacks, writes COG + metadata JSON, flags anything over 50 MB for GitHub release.
- `R/palettes.R` — centralised visParams matching the pre-existing tmap chunks.
- `R/local_rasters.R` — loader (`local_raster()`, `local_aoi()`), derivative helpers (`local_burn_edges()`, `local_fire_categories()`, `local_forest_strata()`), tmap v4 wrappers (`tm_local_categorical`, `tm_local_mask`, `tm_local_continuous`).

### `.qmd` changes

Only the final `tmap::tm_shape(...) + tm_tiles(...)` expression inside each map chunk was rewritten. Upstream GEE computation inside those chunks still runs when the author re-renders with live auth — it is required for analytical outputs the chunks print (pixel counts, area summaries). The rewritten block is wrapped in `local({ })` so any terra variables it creates do not leak into downstream chunks that still use ee-typed objects of the same name.

## 5. Assumptions and limitations

- **Author still needs GEE auth to re-render chapters 1–4.** Analytical chunks (`filterQuality`, `reduceRegion`, `sensitivity-analysis`, etc.) continue to call Earth Engine. Only the map-rendering call has been decoupled. Readers viewing the already-rendered HTML need nothing from GEE.
- **Display-only rasters are downsampled.** Landsat dNBR, Hansen GFC, and ESA WorldCover assets are stored at their native 30 m / 10 m scale but pulled to disk at 500 m (Landsat) or 200 m (others) to stay under the 50 MB `getDownloadURL` cap. Asset is canonical; COG is display-only.
- **Two validation COGs (Hansen, ESA WorldCover) not yet committed.** Their exports were still running when this branch was pushed. The chunks that depend on them fall back to a short placeholder notice. Running `Rscript gee/download_cogs.R` again once both assets are `COMPLETED` fetches them with no other changes needed.
- **Peatland screening gate.** Not triggered. This work touched no SOC calculation.
- **CRS.** All rendered rasters are `EPSG:32616` (UTM zone 16N Honduras). `tmap` in `view` mode reprojects automatically for leaflet display.

## 6. How to reproduce

```bash
# One-off: authenticate GEE.
/Users/seamus/Library/Python/3.13/bin/earthengine authenticate

# Export every source layer to persistent user-owned assets (idempotent).
/opt/local/bin/python gee/export_to_asset.py

# Mirror each asset as a Cloud-Optimised GeoTIFF under data/rendered/.
Rscript gee/download_cogs.R

# Render the book.
quarto render
```

`--only layer_id` and `--force` flags are available on both scripts for selective retries.

## 7. Follow-up

- Complete Hansen and ESA COG pulls once their assets finish exporting (re-run `Rscript gee/download_cogs.R`).
- Optional: fully decouple chapters 1–4 analytical chunks from live GEE so the book can render anywhere without auth. Not done in this task; would require rewriting every `reduceRegion` and server-side filter against local terra.
- Optional: add a `.github/workflows/render-book.yml` that triggers on push, runs `gee/download_cogs.R` first, then `quarto render`, and commits updated `_book/`.

---

*End of FEATURE 0001.*
