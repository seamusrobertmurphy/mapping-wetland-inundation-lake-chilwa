# Fetch 1 arc-second (30 m) SRTM over the Lake Chilwa basin.
# Replaces output_SRTM15Plus.tif, which is a 15 arc-second OCEAN BATHYMETRY
# product (Tozer et al. 2019) carrying land at ~450 m. Wrong product, wrong
# scale, and the reason the derived basin came out at 17,318 km2.
if (Sys.getenv("PROJ_DATA") == "") {
  .pd <- system.file("proj", package = "sf")
  if (file.exists(file.path(.pd, "proj.db"))) Sys.setenv(PROJ_DATA=.pd, PROJ_LIB=.pd)
}
suppressPackageStartupMessages({library(rgee); library(sf); library(terra)})
ee_py <- Sys.getenv("EARTHENGINE_PYTHON", unset = Sys.getenv("RETICULATE_PYTHON"))
Sys.setenv(RETICULATE_PYTHON = ee_py); reticulate::use_python(ee_py, required = TRUE)
ee_Initialize(user="seamusrobertmurphy", drive=FALSE, project="murphys-deforisk", quiet=TRUE)

basin <- st_read("03.outputs/JSON/chilwa_basin.geojson", quiet = TRUE)
bb <- as.numeric(st_bbox(basin)); pad <- 0.05
w<-bb[1]-pad; s<-bb[2]-pad; e<-bb[3]+pad; n<-bb[4]+pad
message(sprintf("bbox %.3f %.3f %.3f %.3f  (%.2f x %.2f deg)", w,s,e,n, e-w, n-s))

srtm <- ee$Image("USGS/SRTMGL1_003")$select("elevation")
out  <- "03.outputs/DEM/extracted"; dir.create(out, showWarnings=FALSE, recursive=TRUE)
work <- file.path(tempdir(), "srtm"); dir.create(work, showWarnings=FALSE)

# 3 x 3 tiles keeps each request well inside the getDownloadURL size limit.
NT <- 3L; xs <- seq(w, e, length.out=NT+1); ys <- seq(s, n, length.out=NT+1)
tiles <- character(0)
for (i in 1:NT) for (j in 1:NT) {
  g <- ee$Geometry$Rectangle(list(xs[i], ys[j], xs[i+1], ys[j+1]))
  url <- srtm$getDownloadURL(list(scale=30, crs="EPSG:4326", region=g,
                                  format="GEO_TIFF"))
  f <- file.path(work, sprintf("t_%d_%d.tif", i, j))
  download.file(url, f, quiet=TRUE, mode="wb")
  tiles <- c(tiles, f)
  message(sprintf("  tile %d/%d  %5.1f MB", length(tiles), NT*NT, file.size(f)/1e6))
}
m <- do.call(terra::merge, lapply(tiles, terra::rast))
dst <- file.path(out, "srtm30m_chilwa.tif")
terra::writeRaster(m, dst, overwrite=TRUE, datatype="INT2S",
                   gdal=c("COMPRESS=DEFLATE","TILED=YES"))
r <- terra::rast(dst)
message(sprintf("\nWROTE %s\n  %d x %d cells (%.1f million)\n  pixel %.6f deg (~%.0f m)\n  elevation %.0f to %.0f m\n  file %.1f MB",
        dst, ncol(r), nrow(r), ncell(r)/1e6, terra::res(r)[1],
        terra::res(r)[1]*111320, min(terra::values(r),na.rm=TRUE),
        max(terra::values(r),na.rm=TRUE), file.size(dst)/1e6))
