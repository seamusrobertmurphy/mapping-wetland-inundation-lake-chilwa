# rgee Workflow Design for the Lake Chilwa Analysis

**Date:** 5 July 2026
**Goal:** control the whole Earth Engine analysis from R via rgee, at full speed, with credentials hidden in `.Renviron`.

## Short answer

rgee is an R wrapper over the same Python `earthengine-api` that the manuscript already calls, so switching to it costs almost no speed and buys R-native syntax, `sf`/`terra`/`stars` interoperability, an interactive `Map`, and managed export tasks. The one real change you must make is authentication: Earth Engine uses an OAuth credential from `ee_Authenticate()`, not the `GOOGLE_API_KEY` you have in `.Renviron`. Re-authenticate once and everything runs from R.

## 1. Authentication model (the important part)

Earth Engine compute authorises through one of two mechanisms: an OAuth user credential created by `rgee::ee_Authenticate()` (interactive, browser, stored locally), or a service-account JSON. It is always tied to a Google Cloud **project** that has the Earth Engine API enabled. A `GOOGLE_API_KEY` authenticates other Google APIs (Maps, Places) and is **not accepted** by Earth Engine, so the key in your `.Renviron` plays no part here and can stay for its other uses.

Because your usual EE credentials are absent, the one-time path is: enable the Earth Engine API on a Cloud project, put its ID in `.Renviron` as `EE_PROJECT`, then in R run `rgee::ee_Authenticate()` (browser flow) followed by `ee_Initialize(project = Sys.getenv("EE_PROJECT"))`. After that, initialisation is silent and no secret ever enters the repository.

## 2. Environment variables (.Renviron)

Your existing `.Renviron` already sets `RETICULATE_PYTHON=/opt/local/bin/python`, which is exactly what rgee needs to find the Python that holds `earthengine-api`. Add one line, `EE_PROJECT=your-project-id`. Keep `TAUDEM_PATH`, `PROJ_LIB`, `R_LIBS_USER`, and `R_C_STACK_LIMIT` as they are; they serve the terrain and spatial stack, not EE. `.Renviron` is git-ignored, so none of this is committed.

## 3. Package setup

Install rgee (`install.packages("rgee")`), then let it verify the Python side with `ee_check()`; if `earthengine-api` or `numpy` are missing from the `RETICULATE_PYTHON` interpreter, `ee_install_upgrade()` or `reticulate::py_install("earthengine-api", pip = TRUE)` adds them. The manuscript's `environment-setup` chunk now loads rgee and calls `ee_Initialize(project = Sys.getenv("EE_PROJECT"))`; the `ee` object and rgee helpers become available to every downstream chunk, and the existing `ee$ImageCollection(...)` calls keep working unchanged.

## 4. What rgee replaces in the current notebook

The manuscript carries hand-rolled helpers that rgee supersedes, and adopting them is where the ergonomic gain is:

The custom `ee_to_sf()` (a `geojsonsf` round-trip through `getInfo`) becomes `ee_as_sf(x)`, which also supports `via = "drive"` or `"gcs"` for geometries too large to pull through `getInfo`. The `ee_tile_url()` plus `tmap::tm_tiles` hack for showing EE layers becomes `Map$addLayer(image, vis, "name")` and `Map$addLayers(...)`, a Leaflet map driven entirely from R. Pulling rasters to disk, which the notebook does not yet do cleanly, becomes `ee_as_rast()` / `ee_as_stars()` (small extents) or an export task for large ones. Table and image exports become `ee_table_to_drive()`, `ee_image_to_drive()`, and `ee_imagecollection_to_drive()`, monitored from R with `ee_monitoring()`.

## 5. Speed: keep it server-side

rgee does not slow Earth Engine; `getInfo()` round-trips and client-side loops do. The manuscript's monthly and annual series are currently built with a `getInfo()` inside an R `lapply` over 96 or 41 iterations, which makes 96 or 41 separate server calls. The faster pattern is to compute the whole reduction server-side and pull it once: map the per-image reducer over the `ImageCollection`, then a single `ee_extract()` or one `getInfo()` on the resulting `FeatureCollection`. For per-geometry statistics use `ee_extract(image, region, fun)` or `reduceRegions` rather than a loop of `reduceRegion`. For anything large, prefer an export task over a synchronous pull, set an explicit `scale` and `maxPixels`, and use `tileScale` or `bestEffort` on heavy `reduceRegion` calls. With those habits, the R-first workflow runs at the same speed as a Python one.

## 6. Migration plan for the manuscript (low-risk order)

The init is already migrated. The safe next steps, in order, are: swap `ee_to_sf()` call sites for `ee_as_sf()`; replace the `ee_tile_url` and `tm_tiles` map blocks with `Map$addLayer` (this also removes the fragile tile-URL helper); convert the Sentinel-1 monthly and Landsat annual `lapply`+`getInfo` loops to a single mapped reducer plus one `ee_extract`, which is both cleaner and much faster; and route the CSV and any raster exports through `ee_table_to_drive` / `ee_image_to_drive`. Each step is independent and testable on its own, and none changes the science, only the plumbing. I can carry these out on the master under the usual backup-and-log protocol whenever you want; they are held until you have authenticated and confirmed the smoke test in `00-gee-setup.qmd` runs.

## 7. Your immediate to-do

Add `EE_PROJECT` to `.Renviron` and restart R; run `ee_Authenticate()` once; run `05.scripts/00-gee-setup.qmd` and confirm the Sentinel-1 and Landsat scene counts print. Tell me it works and I will begin the rgee migration of the manuscript chunks.
