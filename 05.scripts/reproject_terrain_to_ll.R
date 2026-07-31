# Geographic (EPSG:4326) copies of the 30 m terrain grids, for the manuscript
# map chunks, which crop using the AOI extent in degrees. Routing itself stays
# in EPSG:32736 so slope, distance, and contributing area remain metric.
if (Sys.getenv("PROJ_DATA") == "") {
  .pd <- system.file("proj", package = "sf")
  if (file.exists(file.path(.pd, "proj.db"))) Sys.setenv(PROJ_DATA=.pd, PROJ_LIB=.pd)
}
suppressPackageStartupMessages(library(terra)); terraOptions(memfrac = 0.55, progress = 0)
out <- "03.outputs/DEM/extracted"
pairs <- list(c("srtm30m_utm.tif","srtm30m_ll.tif","bilinear"),
              c("dinf30m_FlowDirectionang.tif","dinf30m_FlowDirectionang_ll.tif","near"),
              c("dinf30m_areasca.tif","dinf30m_areasca_ll.tif","bilinear"))
for (p in pairs) {
  fi <- file.path(out, p[1]); fo <- file.path(out, p[2])
  if (!file.exists(fi)) { message("MISSING ", fi); next }
  terra::writeRaster(terra::project(terra::rast(fi), "EPSG:4326", method = p[3]),
                     fo, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE","TILED=YES"))
  x <- terra::rast(fo)
  message(sprintf("  %-34s %d x %d  %.6f deg (~%.0f m)", p[2], ncol(x), nrow(x),
                  terra::res(x)[1], terra::res(x)[1]*111320))
}
