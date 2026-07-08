#!/usr/bin/env Rscript
# =============================================================================
# Lake Chilwa SAR + Landsat time-series pipeline (Google Earth Engine)
# Consolidated, parameterised, runnable version of the manuscript's GEE chunks.
# Built 2026-07-05 for the Journal of Hydrology submission.
#
# WHAT THIS DOES
#   1. Sources the AOI from HydroSHEDS level-6 basins around a point in the lake.
#   2. Sentinel-1 C-band: speckle filter, ratio bands, adaptive water threshold,
#      monthly composites, basin-mean time series -> CSV.
#   3. Landsat C2 (L5/L7/L8): cloud mask, harmonise, five water indices,
#      annual composites, basin-mean time series -> CSV.
#   4. Spectral mixture analysis: annual open-water and flooded-vegetation
#      fractions -> CSV.
#
# WHAT THIS DOES NOT DO
#   It does not authenticate for you. You must have run ee$Authenticate() once
#   under your own Google Earth Engine account, and the project below must be one
#   you can access. It does not build the socio-hydrological model (next phase).
#
# HOW TO RUN
#   Rscript 05.scripts/gee_sar_landsat_pipeline.R
#   or step through interactively in RStudio. See the operational workplan:
#   _staging/gee_pipeline_operational_workplan.md
# =============================================================================

# ------------------------- CONFIG (edit here) --------------------------------
cfg <- list(
  ee_project   = Sys.getenv("EE_PROJECT", unset = "murphys-deforisk"),
  python       = Sys.getenv("RETICULATE_PYTHON", unset = ""),  # "" = let reticulate choose
  lon          = 35.55,     # point inside Lake Chilwa (lon)
  lat          = -15.35,    # point inside Lake Chilwa (lat)
  hybas_level  = "hybas_6", # HydroSHEDS basin level for the AOI
  s1_years     = 2015:2024, # Sentinel-1 monthly-composite window
  ls_years     = 1984:2024, # Landsat annual-composite window
  cloud_max    = 30,        # Landsat scene cloud-cover cap (percent)
  scale_export = 30L,       # reduceRegion scale (m)
  out_dir      = "03.outputs" # relative to repo root
)
# -----------------------------------------------------------------------------

suppressMessages({
  library(reticulate)
  library(sf)
})
if (nzchar(cfg$python)) reticulate::use_python(cfg$python, required = TRUE)

message("Initialising Earth Engine (project: ", cfg$ee_project, ") ...")
ee <- reticulate::import("ee")
ee$Initialize(project = cfg$ee_project)

if (!dir.exists(cfg$out_dir)) dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

ee_to_sf <- function(fc) {
  sf::st_read(jsonlite::toJSON(fc$getInfo(), auto_unbox = TRUE), quiet = TRUE)
}

# --------------------------- 1. AOI ------------------------------------------
message("Sourcing AOI from HydroSHEDS ", cfg$hybas_level, " ...")
poi     <- ee$Geometry$Point(c(cfg$lon, cfg$lat))
basins  <- ee$FeatureCollection(paste0("WWF/HydroSHEDS/v1/Basins/", cfg$hybas_level))$filterBounds(poi)
aoi_ee  <- basins$geometry()

# --------------------------- 2. Sentinel-1 -----------------------------------
message("Building Sentinel-1 monthly water-fraction series ...")
s1 <- ee$ImageCollection("COPERNICUS/S1_GRD")$
  filterBounds(aoi_ee)$
  filter(ee$Filter$eq("instrumentMode", "IW"))$
  filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV"))$
  filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VH"))$
  filter(ee$Filter$eq("orbitProperties_pass", "DESCENDING"))$
  select(c("VV", "VH", "angle"))
message("  Sentinel-1 scenes: ", s1$size()$getInfo())

apply_speckle <- function(image) {
  vv <- image$select("VV")$focal_mean(radius = 3.5, kernelType = "square", units = "pixels")$rename("VV_f")
  vh <- image$select("VH")$focal_mean(radius = 3.5, kernelType = "square", units = "pixels")$rename("VH_f")
  image$addBands(vv)$addBands(vh)$copyProperties(image, list("system:time_start"))
}
s1f <- s1$map(apply_speckle)

# Adaptive VV water threshold (15th percentile), with a -15 dB fallback
sample_vv   <- ee$Image(s1f$first())$select("VV_f")
p15         <- sample_vv$reduceRegion(reducer = ee$Reducer$percentile(list(15L)),
                 geometry = aoi_ee, scale = cfg$scale_export, maxPixels = 1e9)$getInfo()
thr_vv      <- p15$VV_f
if (is.null(thr_vv) || thr_vv > -10 || thr_vv < -25) thr_vv <- -15
message("  Sentinel-1 water threshold (VV): ", round(thr_vv, 2), " dB")

classify_water <- function(image) {
  w <- image$select("VV_f")$lt(thr_vv)$rename("water")
  image$addBands(w)$copyProperties(image, list("system:time_start"))
}
s1w <- s1f$map(classify_water)

ym <- expand.grid(year = as.integer(cfg$s1_years), month = 1:12)
s1_rows <- lapply(seq_len(nrow(ym)), function(i) {
  y <- as.integer(ym$year[i]); m <- as.integer(ym$month[i])
  start <- ee$Date$fromYMD(y, m, 1L); end <- start$advance(1L, "month")
  mc <- s1w$filterDate(start, end)
  n  <- mc$size()$getInfo()
  if (n == 0L) return(NULL)
  stats <- mc$mean()$select(c("VV_f", "VH_f", "water"))$reduceRegion(
    reducer = ee$Reducer$mean(), geometry = aoi_ee, scale = cfg$scale_export, maxPixels = 1e9)$getInfo()
  data.frame(date = sprintf("%04d-%02d-01", y, m), year = y, month = m, n_scenes = n,
             VV_mean = stats$VV_f, VH_mean = stats$VH_f, water_frac = stats$water)
})
s1_df <- do.call(rbind, Filter(Negate(is.null), s1_rows))
write.csv(s1_df, file.path(cfg$out_dir, "s1_monthly_timeseries.csv"), row.names = FALSE)
message("  Wrote ", file.path(cfg$out_dir, "s1_monthly_timeseries.csv"))

# --------------------------- 3. Landsat --------------------------------------
message("Building Landsat annual index series ...")
l5b <- list(from = c("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7"),
            to   = c("blue","green","red","nir","swir1","swir2"))
l8b <- list(from = c("SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7"),
            to   = c("blue","green","red","nir","swir1","swir2"))
scale_l <- function(image) {
  o <- image$select("SR_B.*")$multiply(0.0000275)$add(-0.2)
  image$addBands(o, NULL, TRUE)$copyProperties(image, list("system:time_start"))
}
mask_l <- function(image) {
  qa <- image$select("QA_PIXEL")
  m  <- qa$bitwiseAnd(bitwShiftL(1L, 3L))$eq(0L)$And(qa$bitwiseAnd(bitwShiftL(1L, 4L))$eq(0L))
  image$updateMask(m)$copyProperties(image, list("system:time_start"))
}
harm <- function(image, from, to) image$select(from, to)$copyProperties(image, list("system:time_start"))
build_col <- function(id, bands) {
  ee$ImageCollection(id)$filterBounds(aoi_ee)$filter(ee$Filter$lt("CLOUD_COVER", cfg$cloud_max))$
    map(mask_l)$map(scale_l)$map(function(img) harm(img, bands$from, bands$to))
}
landsat <- build_col("LANDSAT/LT05/C02/T1_L2", l5b)$
  merge(build_col("LANDSAT/LE07/C02/T1_L2", l5b))$
  merge(build_col("LANDSAT/LC08/C02/T1_L2", l8b))$sort("system:time_start")
message("  Landsat scenes: ", landsat$size()$getInfo())

indices <- function(image) {
  g <- image$select("green"); n <- image$select("nir"); s1b <- image$select("swir1")
  r <- image$select("red"); b <- image$select("blue"); s2b <- image$select("swir2")
  ndwi  <- g$subtract(n)$divide(g$add(n))$rename("NDWI")
  mndwi <- g$subtract(s1b)$divide(g$add(s1b))$rename("MNDWI")
  aweish<- b$add(g$multiply(2.5))$subtract(n$add(s1b)$multiply(1.5))$subtract(s2b$multiply(0.25))$rename("AWEIsh")
  wri   <- g$add(r)$divide(n$add(s1b))$rename("WRI")
  ndpi  <- s1b$subtract(g)$divide(s1b$add(g))$rename("NDPI")
  image$addBands(c(ndwi, mndwi, aweish, wri, ndpi))$copyProperties(image, list("system:time_start"))
}
landsat_idx <- landsat$map(indices)

idx_bands <- c("NDWI","MNDWI","AWEIsh","WRI","NDPI","blue","green","red","nir","swir1","swir2")
ls_rows <- lapply(cfg$ls_years, function(y) {
  start <- ee$Date$fromYMD(as.integer(y), 1L, 1L); end <- ee$Date$fromYMD(as.integer(y), 12L, 31L)
  yc <- landsat_idx$filterDate(start, end)
  n  <- yc$size()$getInfo()
  if (n == 0L) return(NULL)
  comp  <- yc$select(idx_bands)$median()
  stats <- comp$reduceRegion(reducer = ee$Reducer$mean(), geometry = aoi_ee,
             scale = cfg$scale_export, maxPixels = 1e9)$getInfo()
  data.frame(year = y, n_scenes = n,
             NDWI = stats$NDWI, MNDWI = stats$MNDWI, AWEIsh = stats$AWEIsh,
             WRI = stats$WRI, NDPI = stats$NDPI)
})
ls_df <- do.call(rbind, Filter(Negate(is.null), ls_rows))
write.csv(ls_df, file.path(cfg$out_dir, "landsat_annual_indices.csv"), row.names = FALSE)
message("  Wrote ", file.path(cfg$out_dir, "landsat_annual_indices.csv"))

# --------------------------- 4. Spectral mixture analysis --------------------
message("Building spectral-mixture annual fraction series ...")
# Endmembers (surface reflectance) for a linear unmixing; adjust to local samples.
em <- list(
  water        = c(0.03, 0.05, 0.04, 0.02, 0.01, 0.01),
  flooded_veg  = c(0.04, 0.07, 0.05, 0.28, 0.12, 0.05),
  dry_veg      = c(0.06, 0.09, 0.08, 0.32, 0.22, 0.12),
  soil         = c(0.14, 0.18, 0.24, 0.30, 0.38, 0.34)
)
em_list <- ee$List(unname(lapply(em, as.list)))
sma_rows <- lapply(cfg$ls_years, function(y) {
  start <- ee$Date$fromYMD(as.integer(y), 1L, 1L); end <- ee$Date$fromYMD(as.integer(y), 12L, 31L)
  yc <- landsat$filterDate(start, end)
  if (yc$size()$getInfo() == 0L) return(NULL)
  comp <- yc$select(c("blue","green","red","nir","swir1","swir2"))$median()
  frac <- comp$unmix(em_list, TRUE, TRUE)$rename(c("f_water","f_floodveg","f_dryveg","f_soil"))
  stats <- frac$reduceRegion(reducer = ee$Reducer$mean(), geometry = aoi_ee,
             scale = cfg$scale_export, maxPixels = 1e9)$getInfo()
  data.frame(year = y, water_frac = stats$f_water, flooded_veg_frac = stats$f_floodveg,
             dry_veg_frac = stats$f_dryveg, soil_frac = stats$f_soil)
})
sma_df <- do.call(rbind, Filter(Negate(is.null), sma_rows))
write.csv(sma_df, file.path(cfg$out_dir, "sma_annual_fractions.csv"), row.names = FALSE)
message("  Wrote ", file.path(cfg$out_dir, "sma_annual_fractions.csv"))

message("Done. Three CSVs written to ", cfg$out_dir, ": ",
        "s1_monthly_timeseries.csv, landsat_annual_indices.csv, sma_annual_fractions.csv")
