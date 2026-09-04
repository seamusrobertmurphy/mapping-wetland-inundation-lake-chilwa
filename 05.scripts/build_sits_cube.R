#!/usr/bin/env Rscript
# Read the downloaded Lake Chilwa sample as a sits data cube.
#
# This is the first step of the local stage. Everything before it ran in Earth
# Engine and returned numbers; from here the pixels are on this machine and the
# analysis is reproducible without a network.
#
# The cube is the 2012 fieldwork year at sixteen-day steps, thirty metres, over
# the lake and a ten-kilometre belt, with four bands. Steps holding no clear
# observation are absent rather than interpolated, so the cube's own timeline
# records where the archive is empty. That is deliberate: a cube that silently
# fills its gaps cannot be told from one that measured them.
#
# Run: Rscript 05.scripts/build_sits_cube.R

suppressPackageStartupMessages({
  library(sits); library(terra); library(sf)
})

root     <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."), mustWork = FALSE)
if (!dir.exists(root)) root <- getwd()
data_dir <- file.path(root, "02.inputs", "CUBE", "chilwa_2012_30m")
if (!dir.exists(data_dir)) data_dir <- "02.inputs/CUBE/chilwa_2012_30m"
stopifnot(dir.exists(data_dir))

files <- list.files(data_dir, pattern = "\\.tif$", full.names = FALSE)
files <- files[!startsWith(files, "._")]
cat(sprintf("%d files in %s\n", length(files), data_dir))
if (length(files) == 0) stop("nothing downloaded yet; run export_sits_cube_sample.py first")

# Filenames are CHILWA_LANDSAT_lake_<BAND>_<YYYY-MM-DD>.tif, five fields split
# on the underscore, which is exactly what parse_info describes.
cube <- sits_cube(
  source     = "BDC",
  collection = "LANDSAT-OLI-16D",
  data_dir   = data_dir,
  parse_info = c("X1", "X2", "tile", "band", "date"),
  delim      = "_",
  bands      = c("GREEN", "NIR", "SWIR1", "MNDWI"),
  progress   = FALSE
)

cat("\n--- cube as sits sees it ---\n")
print(cube)
cat(sprintf("\ntiles      : %d\n", nrow(cube)))
cat(sprintf("bands      : %s\n", paste(sits_bands(cube), collapse = ", ")))
tl <- sits_timeline(cube)
cat(sprintf("timeline   : %d steps, %s to %s\n", length(tl), min(tl), max(tl)))
cat(sprintf("spacing    : %s days\n",
            paste(unique(as.numeric(diff(tl))), collapse = ", ")))
bb <- sits_bbox(cube)
cat(sprintf("bbox       : %.4f to %.4f east, %.4f to %.4f north (%s)\n",
            bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]], bb[["crs"]]))

# Read one band back through terra to confirm the values survived the handover.
r <- terra::rast(file.path(data_dir, files[grep("MNDWI", files)][1]))
v <- terra::values(r); v <- v[!is.na(v) & v != -32768] / 10000
cat(sprintf("\nfirst MNDWI raster: %d by %d, %.1f%% of pixels carry data\n",
            terra::ncol(r), terra::nrow(r),
            100 * length(v) / (terra::ncol(r) * terra::nrow(r))))
cat(sprintf("MNDWI range %.3f to %.3f, median %.3f, water above zero on %.1f%%\n",
            min(v), max(v), median(v), 100 * mean(v > 0)))

saveRDS(cube, file.path(root, "03.outputs", "chilwa_sits_cube_2012.rds"))
cat("\nwrote 03.outputs/chilwa_sits_cube_2012.rds\n")
cat("next: sits_regularize() to force a strict sixteen-day grid, then\n")
cat("      sits_uncertainty() once a classification exists to give the\n")
cat("      per-pixel uncertainty raster this study needs.\n")
