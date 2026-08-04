#!/usr/bin/env Rscript
# Mosaic the tiled ASTER downloads and measure how much of the LAKE each covers.
#
# The measurement is the point. ASTER's 60 km swath means a scene can intersect
# the basin while missing the water entirely: the cloud-free 12 May 2011 pass
# reaches 0.1 per cent of lake pixels. A mosaic is only useful here if it
# actually contains the lake, so the coverage figure is reported before the
# file is offered for use, and mosaics that miss the water are deleted rather
# than left to mislead.

suppressPackageStartupMessages({library(terra); library(sf)})

root <- getwd()
while (!dir.exists(file.path(root, "03.outputs")) && dirname(root) != root) root <- dirname(root)
d <- file.path(root, "03.outputs", "TIF", "RGB_highres")
ql <- file.path(d, "quicklook"); dir.create(ql, showWarnings = FALSE, recursive = TRUE)

lake <- st_read(file.path(root, "02.inputs", "SHP", "lakes_site.shp"), quiet = TRUE) |>
  st_transform(4326) |> st_union() |> vect()

tifs <- list.files(d, pattern = "^ASTER15m_.*_t[0-9]{2}[.]tif$", full.names = TRUE)
prods <- unique(sub("_t[0-9]{2}[.]tif$", "", basename(tifs)))
if (length(prods) == 0) stop("no ASTER tiles found")

cat(sprintf("%-24s %8s %10s %8s\n", "product", "tiles", "lake cov", "size"))
for (p in prods) {
  fs <- tifs[startsWith(basename(tifs), paste0(p, "_t"))]
  fs <- fs[file.info(fs)$size > 100000]      # drop empty off-swath tiles
  if (length(fs) == 0) { cat(sprintf("%-24s %8s %10s\n", p, 0, "no data")); next }
  rs <- lapply(fs, rast)
  m <- if (length(rs) == 1) rs[[1]] else do.call(mosaic, c(rs, list(fun = "max")))
  names(m) <- c("nir", "red", "green")

  inl <- try(mask(crop(m, lake), lake), silent = TRUE)
  cov <- if (inherits(inl, "try-error")) 0 else {
    v <- values(inl[[1]]); 100 * sum(!is.na(v)) / length(v)
  }

  out <- file.path(d, paste0(p, ".tif"))
  writeRaster(m, out, overwrite = TRUE, datatype = "INT1U",
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"))
  mb <- file.info(out)$size / 1048576
  cat(sprintf("%-24s %8d %9.1f%% %7.1f MB\n", p, length(fs), cov, mb))

  if (cov < 20) {
    file.remove(out)
    cat("    removed: covers under 20% of the lake, not usable here\n")
  } else {
    png(file.path(ql, paste0(p, ".png")), width = 850, height = 850)
    try(plotRGB(m, r = 1, g = 2, b = 3, stretch = "lin", mar = c(0, 0, 2, 0),
                main = paste(p, sprintf("(lake %.0f%%)", cov))), silent = TRUE)
    try(plot(lake, add = TRUE, border = "cyan", lwd = 2), silent = TRUE)
    dev.off()
  }
  file.remove(fs)
}
cat("\ndone\n")
