#!/usr/bin/env Rscript
# Download GEE assets as local Cloud-Optimised GeoTIFFs.
#
# Reads gee/layer_manifest.yaml, pulls each asset's pixels via the computePixels
# endpoint (getDownloadURL), and writes a COG to data/rendered/<id>.tif.
# Metadata JSON written alongside.
#
# Usage:
#   Rscript gee/download_cogs.R [--only ID ...] [--force]
#
# Design decisions (see specs/TODO_0001_GEE_Raster_Export.md §6):
#   - Uses rgee::ee_as_raster(via = "getDownloadURL"); no Google Drive.
#   - Target CRS EPSG:32616 (UTM 16N). Book is indifferent to this once loaded.
#   - COGs compressed with DEFLATE + PREDICTOR=2, 256x256 blocks.
#   - Size reported per layer; anything >50 MB flagged for GitHub release.

suppressPackageStartupMessages({
  library(yaml)
  library(reticulate)
  reticulate::use_python("/opt/local/bin/python", required = TRUE)
  library(rgee)
  library(terra)
  library(sf)
  library(jsonlite)
  library(httr)
})

repo_root   <- here::here()
manifest_fp <- file.path(repo_root, "gee", "layer_manifest.yaml")
rendered    <- file.path(repo_root, "data", "rendered")
metadata    <- file.path(rendered, "metadata")

dir.create(rendered, showWarnings = FALSE, recursive = TRUE)
dir.create(metadata, showWarnings = FALSE, recursive = TRUE)

args   <- commandArgs(trailingOnly = TRUE)
force  <- "--force" %in% args
only   <- if ("--only" %in% args) {
  idx <- which(args == "--only")
  args[(idx + 1):length(args)]
} else {
  NULL
}

manifest <- yaml::read_yaml(manifest_fp)

rgee::ee_Initialize(project = manifest$asset_project)

# AOI as ee$Geometry for clipping at download time.
aoi_ee <- ee$FeatureCollection(manifest$aoi$source)$
  filter(ee$Filter$eq("ADM0_NAME", "Honduras"))$
  geometry()

# Helper matching the preface's ee_to_sf — avoids the geojsonio dependency
# that rgee::ee_as_sf pulls in.
ee_fc_to_sf <- function(fc) {
  info <- fc$getInfo()
  geojsonsf::geojson_sf(jsonlite::toJSON(info, auto_unbox = TRUE))
}

# One-off: export AOI as committed GPKG if missing. Exports both country
# (level0) and state (level1) layers into a single multilayer GPKG so the
# book's `local_aoi()` can pull either.
aoi_gpkg <- file.path(repo_root, manifest$aoi$committed_path)
if (!file.exists(aoi_gpkg) || force) {
  message("Writing AOI GPKG...")
  country_fc <- ee$FeatureCollection(manifest$aoi$source)$
    filter(ee$Filter$eq("ADM0_NAME", "Honduras"))
  country <- sf::st_transform(ee_fc_to_sf(country_fc), 32616)
  states  <- ee$FeatureCollection(sub("level0", "level1", manifest$aoi$source))$
    filter(ee$Filter$eq("ADM0_NAME", "Honduras"))
  states_sf <- sf::st_transform(ee_fc_to_sf(states), 32616)
  sf::st_write(country,   aoi_gpkg, layer = "country", delete_dsn = TRUE, quiet = TRUE)
  sf::st_write(states_sf, aoi_gpkg, layer = "states",  append = TRUE,    quiet = TRUE)
}

cog_options <- c(
  "COMPRESS=DEFLATE",
  "PREDICTOR=2",
  "TILED=YES",
  "BLOCKXSIZE=256",
  "BLOCKYSIZE=256",
  "COPY_SRC_OVERVIEWS=YES"
)

# Manifest dtype strings → terra::writeRaster datatype codes.
terra_dtype <- function(dtype) {
  switch(tolower(dtype %||% "uint8"),
    uint8   = "INT1U",
    int8    = "INT1S",
    uint16  = "INT2U",
    int16   = "INT2S",
    uint32  = "INT4U",
    int32   = "INT4S",
    float32 = "FLT4S",
    float64 = "FLT8S",
    "INT1U"
  )
}

asset_ready <- function(asset_id) {
  tryCatch({
    rgee::ee$data$getAsset(asset_id)
    TRUE
  }, error = function(e) FALSE)
}

download_layer <- function(spec) {
  id       <- spec$id
  asset_id <- paste0(manifest$asset_root, "/", id)
  out_fp   <- file.path(rendered, paste0(id, ".tif"))

  if (file.exists(out_fp) && !force) {
    message(sprintf("[skip] %s — COG already exists", id))
    return(invisible(NULL))
  }

  if (!asset_ready(asset_id)) {
    message(sprintf("[wait] %s — asset not ready yet, skipping for now", id))
    return(invisible(NULL))
  }

  message(sprintf("[pull] %s <- %s", id, asset_id))
  img <- ee$Image(asset_id)

  download_scale <- if (isTRUE(spec$display_only)) max(spec$scale_m, 500) else spec$scale_m

  fetch_tif <- function(ee_img) {
    url <- ee_img$getDownloadURL(list(
      region = aoi_ee,
      scale  = download_scale,
      crs    = "EPSG:32616",
      format = "GEO_TIFF"
    ))
    tmp <- tempfile(fileext = ".tif")
    resp <- httr::GET(url, httr::write_disk(tmp, overwrite = TRUE), httr::timeout(600))
    if (httr::http_error(resp)) {
      stop(sprintf("HTTP %d fetching %s", httr::status_code(resp), id), call. = FALSE)
    }
    terra::rast(tmp)
  }

  # Assets with many bands blow past the 50 MB per-request limit. Pull in
  # groups of 6 bands and stack on disk. Reset names from the original
  # GEE band list so years like "burn_2020" survive the chunked round-trip.
  band_names <- img$bandNames()$getInfo()
  if (length(band_names) > 6) {
    chunks <- split(band_names, ceiling(seq_along(band_names) / 6))
    parts <- lapply(chunks, function(bs) fetch_tif(img$select(bs)))
    r <- Reduce(c, parts)
  } else {
    r <- fetch_tif(img)
  }
  names(r) <- band_names
  terra::writeRaster(
    r,
    filename  = out_fp,
    overwrite = TRUE,
    filetype  = "COG",
    gdal      = cog_options,
    datatype  = terra_dtype(spec$dtype)
  )

  size_mb <- file.info(out_fp)$size / 1e6
  message(sprintf("  wrote %s (%.1f MB)", basename(out_fp), size_mb))

  meta <- list(
    id              = id,
    asset_id        = asset_id,
    local_cog       = basename(out_fp),
    size_mb         = round(size_mb, 2),
    over_50mb       = size_mb > 50,
    target_crs      = "EPSG:32616",
    scale_m         = spec$scale_m,
    display_only    = isTRUE(spec$display_only),
    download_utc    = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    rgee_version    = as.character(packageVersion("rgee"))
  )
  jsonlite::write_json(
    meta,
    file.path(metadata, sprintf("%s_%s.json", id, format(Sys.Date(), "%Y%m%d"))),
    auto_unbox = TRUE, pretty = TRUE
  )

  if (size_mb > 50) {
    warning(sprintf(
      "%s exceeds 50 MB (%.1f MB). Move to GitHub release 'gee-exports-v1' with piggyback.",
      id, size_mb
    ))
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

layers <- manifest$layers
if (!is.null(only)) {
  layers <- Filter(function(s) s$id %in% only, layers)
}

oversized <- character()
for (spec in layers) {
  size <- tryCatch(download_layer(spec), error = function(e) {
    message(sprintf("[fail] %s: %s", spec$id, conditionMessage(e)))
    NA_real_
  })
}

message("\nDone. Inspect data/rendered/ for COGs and data/rendered/metadata/ for per-layer JSON.")
