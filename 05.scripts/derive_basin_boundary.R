#!/usr/bin/env Rscript
# Derive the Lake Chilwa basin from the terrain.
#
# Method of record (DECISIONS.md): least-cost depression breaching only, never
# filling; D-infinity flow routing (Tarboton 1997) for the drainage structure,
# with the MFD/FD8 variant (Freeman 1991; Quinn 1991) as the named alternative.
#
# The basin is endorheic, so the divide is the boundary of the region draining
# internally to the terminal depression rather than to the grid edge. The
# terminal sink is located on the unconditioned surface, before any breaching
# can carve an artificial outlet through it.
#
# Writes: 03.outputs/SHP/chilwa_basin_derived.shp
#         03.outputs/JSON/chilwa_basin_provenance.json
#         03.outputs/CSV/basin_derivation_checks.csv

# GDAL on this machine does not find proj.db unless PROJ_DATA points at the
# copy shipped with sf. Set before any spatial package is attached.
if (Sys.getenv("PROJ_DATA") == "") {
  .pd <- system.file("proj", package = "sf")
  if (file.exists(file.path(.pd, "proj.db"))) {
    Sys.setenv(PROJ_DATA = .pd, PROJ_LIB = .pd)
  }
}

suppressPackageStartupMessages({
  library(whitebox); library(terra); library(sf); library(jsonlite)
})

wbt_init()
sf_use_s2(FALSE)  # matches the manuscript's area convention

root <- getwd()
while (!dir.exists(file.path(root, "03.outputs")) && dirname(root) != root) root <- dirname(root)
stopifnot(dir.exists(file.path(root, "03.outputs")))

dem_src  <- file.path(root, "03.outputs", "DEM", "extracted", "output_SRTM15Plus.tif")
work     <- file.path(tempdir(), "chilwa_basin")
dir.create(work, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dem_src))

# ---- 0. Project to metres so slopes, areas, and thresholds are physical -----
message("0. reprojecting to EPSG:32736")
dem_ll <- terra::rast(dem_src)
res_m  <- 450  # native SRTM15+ is ~450 m at this latitude
dem_m  <- terra::project(dem_ll, "EPSG:32736", res = res_m, method = "bilinear")
p_raw  <- file.path(work, "dem_00_raw.tif")
terra::writeRaster(dem_m, p_raw, overwrite = TRUE)

# ---- 1. Terminal point of the internal drainage --------------------------
# Located on the unconditioned surface, before breaching can carve through it.
# For a closed basin the terminal point is the lowest cell of the lake floor.
message("1. locating the terminal point on the unconditioned surface")
raw_r <- terra::rast(p_raw)
min_cell <- terra::where.min(raw_r)
terminal_xy <- terra::xyFromCell(raw_r, min_cell[1, "cell"])
terminal_elev <- min_cell[1, "value"]
message("   terminal cell at ", round(terminal_xy[1]), ", ", round(terminal_xy[2]),
        "  elevation ", round(terminal_elev, 1), " m")

# ---- 2. Conditioning: breaching only, never filling -------------------------
message("2. breaching (single-cell pits, then least-cost, fill = FALSE)")
p_pits <- file.path(work, "dem_01_pits.tif")
p_brch <- file.path(work, "dem_02_breached.tif")
wbt_breach_single_cell_pits(dem = p_raw, output = p_pits)
wbt_breach_depressions_least_cost(
  dem = p_pits, output = p_brch,
  dist = 100, min_dist = TRUE, flat_increment = 0.001, fill = FALSE)

brch_r <- terra::rast(p_brch)
diff_r <- raw_r - brch_r
raised <- sum(terra::values(diff_r) < -1e-6, na.rm = TRUE)
message("   cells raised by conditioning (must be 0): ", raised)
message("   max breach depth (m): ", round(max(terra::values(diff_r), na.rm = TRUE), 3))

# ---- 3. D-infinity drainage structure (method of record) --------------------
message("3. D-infinity pointer and accumulation")
p_dinf_ptr <- file.path(work, "dinfFlowDirectionang.tif")
p_dinf_sca <- file.path(work, "Dinfareasca.tif")
p_dinf_ca  <- file.path(work, "Dinfareaca.tif")
wbt_d_inf_pointer(dem = p_brch, output = p_dinf_ptr)
wbt_d_inf_flow_accumulation(input = p_brch, output = p_dinf_sca, out_type = "Specific Contributing Area")
wbt_d_inf_flow_accumulation(input = p_brch, output = p_dinf_ca,  out_type = "Catchment Area")

# Named alternative, for the dispersion comparison
message("   MFD/FD8 alternative")
p_mdinf <- file.path(work, "MDinfarea.tif")
wbt_md_inf_flow_accumulation(dem = p_brch, output = p_mdinf, out_type = "catchment area")

# ---- 4. The divide: the internally drained region ---------------------------
# Endorheic basin, so the divide is not a pour-point trace to an outlet. Every
# cell either drains internally to the terminal depression or leaves by the grid
# edge. Basins() partitions the grid by catchment; the one containing the
# terminal cell is the internally drained region. The D8 pointer is used for the
# tracing step only and is not a reported product.
message("4. delineating the internally drained region")
p_d8ptr  <- file.path(work, "d8_pointer.tif")
p_basins <- file.path(work, "edge_basins.tif")
wbt_d8_pointer(dem = p_brch, output = p_d8ptr)
wbt_basins(d8_pntr = p_d8ptr, output = p_basins)

basins_r <- terra::rast(p_basins)
basin_id_at_terminal <- terra::extract(basins_r, matrix(terminal_xy, ncol = 2))[1, 1]
message("   basin id containing the terminal cell: ", basin_id_at_terminal)

mask_r <- terra::ifel(basins_r == basin_id_at_terminal, 1, NA)
poly   <- sf::st_as_sf(terra::as.polygons(mask_r, dissolve = TRUE))
poly   <- sf::st_make_valid(sf::st_union(poly))

area_km2_utm <- as.numeric(sf::st_area(poly)) / 1e6
poly_ll <- sf::st_transform(sf::st_sf(geometry = sf::st_sfc(poly, crs = 32736)), 4326)
area_km2_ll  <- as.numeric(sf::st_area(poly_ll)) / 1e6

message("   derived area (EPSG:32736): ", round(area_km2_utm, 1), " km2")
message("   derived area (EPSG:4326, s2 off): ", round(area_km2_ll, 1), " km2")

# ---- 5. Comparison against the inherited HydroBASINS polygon ----------------
inherited_p <- file.path(root, "03.outputs", "SHP", "chilwa_basin.shp")
cmp <- list()
if (file.exists(inherited_p)) {
  inh <- sf::st_transform(sf::read_sf(inherited_p), 32736)
  inh <- sf::st_make_valid(sf::st_union(inh))
  der <- sf::st_sfc(poly, crs = 32736)
  inter <- as.numeric(sf::st_area(sf::st_intersection(der, inh))) / 1e6
  uni   <- as.numeric(sf::st_area(sf::st_union(der, inh))) / 1e6
  cmp <- list(
    inherited_area_km2 = round(as.numeric(sf::st_area(inh)) / 1e6, 1),
    intersection_km2   = round(inter, 1),
    union_km2          = round(uni, 1),
    jaccard            = round(inter / uni, 4),
    area_difference_km2 = round(area_km2_utm - as.numeric(sf::st_area(inh)) / 1e6, 1))
  message("   Jaccard vs inherited HydroBASINS: ", cmp$jaccard,
          "   area difference: ", cmp$area_difference_km2, " km2")
}

# ---- 6. Write outputs with a self-derived schema ----------------------------
message("6. writing outputs")
out_sf <- sf::st_sf(
  basin_id    = 1L,
  dem_src     = "SRTM15+ (03.outputs/DEM)",
  dem_res_m   = res_m,
  conditioning = "BreachSingleCellPits + BreachDepressionsLeastCost(fill=FALSE)",
  route_struct = "Dinf",
  route_alt    = "MDInf (FD8)",
  divide_rule  = "internally drained region containing the terminal cell",
  terminal_x = round(terminal_xy[1], 1),
  terminal_y = round(terminal_xy[2], 1),
  area_km2    = round(area_km2_utm, 2),
  crs_calc    = "EPSG:32736",
  wbt_ver     = as.character(whitebox::wbt_version()[1]),
  geometry    = sf::st_sfc(poly, crs = 32736))

dir.create(file.path(root, "03.outputs", "JSON"), showWarnings = FALSE, recursive = TRUE)
sf::st_write(sf::st_transform(out_sf, 4326),
             file.path(root, "03.outputs", "SHP", "chilwa_basin_derived.shp"),
             delete_dsn = TRUE, quiet = TRUE)

prov <- c(list(
  derived_utc        = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  dem_source         = dem_src,
  dem_res_m          = res_m,
  conditioning       = "breaching only, fill = FALSE",
  cells_raised       = raised,
  max_breach_depth_m = round(max(terra::values(diff_r), na.rm = TRUE), 3),
  routing_structure  = "D-infinity (Tarboton 1997)",
  routing_alternative = "MDInf / FD8 (Freeman 1991; Quinn 1991)",
  divide_rule        = "internally drained region containing the terminal cell",
  terminal_cell_xy   = round(as.numeric(terminal_xy), 1),
  terminal_elev_m    = round(as.numeric(terminal_elev), 2),
  area_km2_utm       = round(area_km2_utm, 2),
  area_km2_wgs84_s2off = round(area_km2_ll, 2),
  hydrobasins_fields_present = FALSE), cmp)

write_json(prov, file.path(root, "03.outputs", "JSON", "chilwa_basin_provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)

checks <- data.frame(
  check = c("cells raised by conditioning", "max breach depth (m)",
            "derived area km2 (EPSG:32736)", "derived area km2 (EPSG:4326, s2 off)",
            "jaccard vs inherited"),
  value = c(raised, round(max(terra::values(diff_r), na.rm = TRUE), 3),
            round(area_km2_utm, 2), round(area_km2_ll, 2),
            if (length(cmp)) cmp$jaccard else NA))
write.csv(checks, file.path(root, "03.outputs", "CSV", "basin_derivation_checks.csv"), row.names = FALSE)

message("done")
