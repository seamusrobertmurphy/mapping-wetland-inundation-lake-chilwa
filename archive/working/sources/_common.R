# Shared setup sourced by all computational pages.
# Loads packages, configures knitr, authenticates GEE, and derives AOI.

easypackages::packages(
  "bslib",
  "cols4all", "covr", "cowplot",
  "dendextend", "digest", "DiagrammeR", "dtwclust", "downlit",
  "e1071", "exactextractr", "elevatr",
  "FNN", "future", "forestdata",
  "gdalcubes", "gdalUtilities", "geojsonsf", "geos", "ggplot2", "ggstats",
  "ggspatial", "ggmap", "ggplotify", "ggpubr", "ggrepel", "giscoR",
  "hdf5r", "httr", "httr2", "htmltools",
  "jsonlite",
  "kohonen",
  "leaflet.providers", "leafem", "libgeos", "luz", "lwgeom", "leaflet", "leafgl",
  "mapedit", "mapview", "maptiles", "methods", "mgcv",
  "ncdf4", "nnet",
  "openxlsx", "parallel", "plotly",
  "randomForest", "rasterVis", "raster", "Rcpp", "RcppArmadillo",
  "RcppCensSpatial", "rayshader", "RcppEigen", "RcppParallel",
  "RColorBrewer", "reactable", "rgl", "rsconnect", "RStoolbox", "rts",
  "s2", "sf", "scales", "sits", "spdep", "stars", "stringr", "supercells",
  "terra", "testthat", "tidyverse", "tidyterra", "tools",
  "tmap", "tmaptools", "terrainr",
  "webshot2",
  "xgboost",
  prompt = FALSE)

sf::sf_use_s2(use_s2 = FALSE)

knitr::opts_knit$set(root.dir = here::here())
options(repos = c(CRAN = "https://cloud.r-project.org"),
        htmltools.dir.version = FALSE,
        htmltools.preserve.raw = FALSE)

knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  error = FALSE, comment = NA,
  tidy.opts = list(width.cutoff = 60))

# GEE
reticulate::use_python("/opt/local/bin/python", required = TRUE)
ee <- reticulate::import("ee")
ee$Initialize(project = "murphys-deforisk")

# Helpers
ee_tile_url <- function(ee_image, vis_params) {
  ee_image$getMapId(vis_params)$tile_fetcher$url_format
}

ee_to_sf <- function(fc) {
  geojsonsf::geojson_sf(
    jsonlite::toJSON(fc$getInfo(), auto_unbox = TRUE))
}

# AOI
lon <- 35.55
lat <- -15.35
poi <- ee$Geometry$Point(c(lon, lat))
poi_sf <- st_sfc(st_point(c(lon, lat)), crs = 4326)

basin_l6 <- ee$FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_6")$filterBounds(poi)
basin_l8 <- ee$FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_8")$filterBounds(poi)
basin_l6_sf <- ee_to_sf(basin_l6)
basin_l8_sf <- ee_to_sf(basin_l8)

aoi_sf <- st_union(basin_l6_sf)
aoi_ee <- ee$FeatureCollection(basin_l6)$geometry()
aoi_bbox <- st_bbox(aoi_sf)

aoi_country_ee <- ee$FeatureCollection("FAO/GAUL/2015/level0")$filter(ee$Filter$eq("ADM0_NAME", "Malawi"))
aoi_country_sf <- ee_to_sf(aoi_country_ee)
