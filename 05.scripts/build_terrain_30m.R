#!/usr/bin/env Rscript
# D-infinity terrain grids at 1 arc-second (30 m) for the Lake Chilwa basin.
# Replaces the 15 arc-second SRTM15Plus products (~450 m), which are an OCEAN
# BATHYMETRY dataset and the wrong input for a terrestrial catchment.
# Conditioning is breaching only, never filling (DECISIONS.md).
# Outputs carry a 30m_ prefix so the legacy tar archives cannot overwrite them.
if (Sys.getenv("PROJ_DATA") == "") {
  .pd <- system.file("proj", package = "sf")
  if (file.exists(file.path(.pd, "proj.db"))) Sys.setenv(PROJ_DATA=.pd, PROJ_LIB=.pd)
}
suppressPackageStartupMessages({library(whitebox); library(terra); library(sf)})
wbt_init(); terraOptions(memfrac = 0.55, progress = 0)

out  <- "03.outputs/DEM/extracted"
src  <- file.path(out, "srtm30m_chilwa.tif")
stopifnot(file.exists(src))

p_utm <- file.path(out, "srtm30m_utm.tif")
if (!file.exists(p_utm)) {
  message("[1/5] reprojecting to EPSG:32736 at 30 m")
  terra::writeRaster(terra::project(terra::rast(src), "EPSG:32736", res = 30,
                     method = "bilinear"), p_utm, overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE","TILED=YES"))
}
r <- terra::rast(p_utm)
message(sprintf("      %d x %d = %.1f M cells at %.0f m", ncol(r), nrow(r),
                ncell(r)/1e6, terra::res(r)[1]))

p_pit <- file.path(out, "srtm30m_pits.tif")
p_brc <- file.path(out, "srtm30m_breached.tif")
message("[2/5] breaching single-cell pits")
if (!file.exists(p_pit)) wbt_breach_single_cell_pits(dem = p_utm, output = p_pit)
message("[3/5] least-cost breaching (fill = FALSE)")
if (!file.exists(p_brc))
  wbt_breach_depressions_least_cost(dem = p_pit, output = p_brc, dist = 200,
                                    min_dist = TRUE, flat_increment = 0.001,
                                    fill = FALSE)

message("[4/5] D-infinity pointer")
wbt_d_inf_pointer(dem = p_brc, output = file.path(out, "dinf30m_FlowDirectionang.tif"))
message("[5/5] D-infinity specific contributing area")
wbt_d_inf_flow_accumulation(input = p_brc,
  output = file.path(out, "dinf30m_areasca.tif"),
  out_type = "Specific Contributing Area")

for (f in c("srtm30m_utm.tif","dinf30m_FlowDirectionang.tif","dinf30m_areasca.tif")) {
  x <- terra::rast(file.path(out, f))
  message(sprintf("  WROTE %-32s %d x %d  %.0f m", f, ncol(x), nrow(x), terra::res(x)[1]))
}
message("done")
