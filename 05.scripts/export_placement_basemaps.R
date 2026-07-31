# Static RGB basemaps for the manual ground-truth placement tool.
#
# The field photographs run from December 2011 to November 2014. Sentinel-2 does
# not reach them: the first Sentinel-2 scene over this basin is 4 September 2015,
# and there is no Sentinel-2 at all before 2015. Landsat 5 stopped imaging the
# basin on 18 October 2011 and Landsat 7 is excluded for scan-line-corrector loss,
# so 2012 has no optical coverage either. The nearest usable composites bracketing
# the fieldwork are therefore Landsat 5 for 2011 and Landsat 8 for 2013.
#
# Writes PNGs plus their bounds to 03.outputs/PNG/placement/.

suppressMessages({library(rgee); library(sf); library(jsonlite)})
ee_py <- Sys.getenv("EARTHENGINE_PYTHON", unset = Sys.getenv("RETICULATE_PYTHON"))
Sys.setenv(RETICULATE_PYTHON = ee_py); reticulate::use_python(ee_py, required = TRUE)
ee_Initialize(user = "seamusrobertmurphy", drive = FALSE,
              project = "murphys-deforisk", quiet = TRUE)

basin <- st_read(here::here("03.outputs","JSON","chilwa_basin.geojson"), quiet = TRUE)
bb    <- as.numeric(st_bbox(basin))
aoi   <- ee$Geometry$Rectangle(bb[1], bb[2], bb[3], bb[4])
out   <- here::here("03.outputs","PNG","placement")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# Collection 2 Level 2 surface reflectance: scale 0.0000275, offset -0.2.
sr <- function(img) img$select("SR_B.")$multiply(0.0000275)$add(-0.2)

grab <- function(id, s, e, rgb, label) {
  ic <- ee$ImageCollection(id)$filterBounds(aoi)$filterDate(s, e)$
    filter(ee$Filter$lt("CLOUD_COVER", 25))
  n <- ic$size()$getInfo()
  if (n == 0) { cat(sprintf("  %s: no scenes\n", label)); return(invisible(NULL)) }
  comp <- ic$map(ee_utils_pyfunc(function(i) sr(ee$Image(i))))$median()$clip(aoi)
  url <- comp$visualize(bands = rgb, min = 0.02, max = 0.32, gamma = 1.15)$
    getThumbURL(list(region = aoi, dimensions = 2048, format = "png"))
  f <- file.path(out, paste0(label, ".png"))
  download.file(url, f, quiet = TRUE, mode = "wb")
  cat(sprintf("  %-22s n=%-3d -> %s\n", label, n, basename(f)))
  invisible(NULL)
}

cat("exporting placement basemaps\n")
grab("LANDSAT/LT05/C02/T1_L2","2011-01-01","2011-12-31", c("SR_B3","SR_B2","SR_B1"), "landsat5_2011_rgb")
grab("LANDSAT/LC08/C02/T1_L2","2013-01-01","2013-12-31", c("SR_B4","SR_B3","SR_B2"), "landsat8_2013_rgb")

write_json(list(
  crs    = "EPSG:4326",
  bounds = list(west = bb[1], south = bb[2], east = bb[3], north = bb[4]),
  note   = paste("Sentinel-2 does not cover the 2011-2014 photographs; first S2 scene",
                 "over this basin is 2015-09-04. Landsat 5 ends 2011-10-18 and 2012 has",
                 "no optical coverage, so these bracket the fieldwork.")),
  file.path(out, "bounds.json"), auto_unbox = TRUE, pretty = TRUE)
cat("bounds.json written\n")
