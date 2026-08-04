#!/usr/bin/env Rscript
# Mosaic the tiled high-resolution RGB downloads and write quick-look previews.
#
# Earth Engine caps a single download at 50 MB, so export_highres_groundtruth.py
# fetches each product as a tile grid. This joins them, compresses with DEFLATE
# so nothing approaches GitHub's 100 MB per-file limit, and writes a PNG preview
# of each so the imagery can be checked without loading the rasters.

suppressPackageStartupMessages(library(terra))

root <- getwd()
while (!dir.exists(file.path(root, "03.outputs")) && dirname(root) != root) {
  root <- dirname(root)
}
d <- file.path(root, "03.outputs", "TIF", "RGB_highres")
ql <- file.path(d, "quicklook")
dir.create(ql, showWarnings = FALSE, recursive = TRUE)

tifs <- list.files(d, pattern = "_t[0-9]{2}[.]tif$", full.names = TRUE)
if (length(tifs) == 0) stop("no tiles found in ", d)

products <- unique(sub("_t[0-9]{2}[.]tif$", "", basename(tifs)))
cat("products:", paste(products, collapse = ", "), "\n\n")

for (p in products) {
  fs <- tifs[startsWith(basename(tifs), paste0(p, "_t"))]
  cat(p, ":", length(fs), "tiles\n")
  rs <- lapply(fs, rast)
  m <- if (length(rs) == 1) rs[[1]] else do.call(mosaic, c(rs, list(fun = "mean")))
  names(m) <- c("red", "green", "blue")[seq_len(nlyr(m))]
  out <- file.path(d, paste0(p, ".tif"))
  writeRaster(m, out, overwrite = TRUE, datatype = "INT2U",
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9", "TILED=YES"))
  mb <- round(file.info(out)$size / 1048576, 1)
  cat("   wrote", basename(out), paste0(mb, " MB"), paste0(ncol(m), " x ", nrow(m)), "\n")

  png(file.path(ql, paste0(p, ".png")), width = 900, height = 850)
  try(plotRGB(m, r = 1, g = 2, b = 3, stretch = "lin", mar = c(0, 0, 2, 0), main = p),
      silent = TRUE)
  dev.off()

  file.remove(fs)
}

cat("\ntiles removed; mosaics and previews in", d, "\n")
