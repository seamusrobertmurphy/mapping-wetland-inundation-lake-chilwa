# Basin-mean climate series for the Lake Chilwa basin, 1984 to 2024.
#
# The window matches the optical and radar record: Landsat from 1984, Sentinel-1
# from 2015, ALOS PALSAR across 2007-2010 and 2015-2020. Three products are read
# so that no single reanalysis or satellite estimate carries the record alone.
#
#   CHIRPS pentad        precipitation, station-blended, 0.05 deg, 1981 onward
#   TerraClimate         precipitation, PET, min and max temperature, PDSI, 1958-2024
#   ERA5-Land monthly    precipitation and 2 m air temperature, 1950 onward
#   CSIC SPEI v2.11      standardised precipitation-evapotranspiration index, 1901-2024
#
# SPEI rather than SPI: the basin is endorheic, water leaves only by evaporation
# and seepage, so the evaporative term belongs in the drought index.
#
# Writes 03.outputs/CSV/climate_monthly.csv and climate_annual.csv.

suppressMessages({library(rgee); library(sf); library(dplyr); library(tidyr)})

YR_START <- 1984L
YR_END   <- 2024L
OUT_DIR  <- here::here("03.outputs", "CSV")

ee_py <- Sys.getenv("EARTHENGINE_PYTHON", unset = Sys.getenv("RETICULATE_PYTHON"))
Sys.setenv(RETICULATE_PYTHON = ee_py)
reticulate::use_python(ee_py, required = TRUE)
ee_Initialize(user = "seamusrobertmurphy", drive = FALSE,
              project = "murphys-deforisk", quiet = TRUE)

basin <- st_read(here::here("03.outputs", "JSON", "chilwa_basin.geojson"), quiet = TRUE)
aoi   <- sf_as_ee(st_geometry(basin))

start <- ee$Date$fromYMD(YR_START, 1L, 1L)
end   <- ee$Date$fromYMD(YR_END + 1L, 1L, 1L)
n_mon <- (YR_END - YR_START + 1L) * 12L
months <- ee$List$sequence(0, n_mon - 1)

# Basin mean of one image, returned as a dated feature. Scale is the product's
# native grid; reduceRegion at a coarser scale than the pixel would alias.
# The date is passed in rather than read off the image: copyProperties and set
# both return an Element, and calling Image.date on that fails through reticulate.
basin_mean <- function(img, scale, date_str) {
  v <- img$reduceRegion(reducer = ee$Reducer$mean(), geometry = aoi,
                        scale = scale, maxPixels = 1e9, bestEffort = TRUE)
  ee$Feature(NULL, v)$set("date", date_str)
}

# Month label of an image already known to carry system:time_start.
month_of <- function(img) ee$Date(img$get("system:time_start"))$format("YYYY-MM")

# Pull a monthly FeatureCollection down in one call and tidy it.
fetch <- function(fc, value_names) {
  raw <- fc$getInfo()$features
  out <- lapply(raw, function(f) {
    p <- f$properties
    row <- list(date = p$date)
    for (nm in value_names) row[[nm]] <- if (is.null(p[[nm]])) NA_real_ else p[[nm]]
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  bind_rows(out)
}

# ---- CHIRPS: sum the pentads falling in each month --------------------------
chirps <- ee$ImageCollection("UCSB-CHG/CHIRPS/PENTAD")$select("precipitation")
chirps_monthly <- ee$FeatureCollection(months$map(ee_utils_pyfunc(function(m) {
  s <- start$advance(m, "month")
  img <- chirps$filterDate(s, s$advance(1, "month"))$sum()$rename("pr_chirps")
  basin_mean(img, 5566, s$format("YYYY-MM"))
})))
chirps_df <- fetch(chirps_monthly, "pr_chirps")

# ---- TerraClimate: already monthly; apply the documented scale factors -------
# pr mm (unscaled), pet 0.1 mm, tmmn/tmmx 0.1 degC, pdsi 0.01.
tc <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")$filterDate(start, end)
tc_monthly <- tc$map(ee_utils_pyfunc(function(img) {
  img <- ee$Image(img)
  scaled <- img$select("pr")$rename("pr_terraclim")
  scaled <- scaled$addBands(img$select("pet")$multiply(0.1)$rename("pet_terraclim"))
  scaled <- scaled$addBands(img$select("tmmx")$multiply(0.1)$rename("tmax_terraclim"))
  scaled <- scaled$addBands(img$select("tmmn")$multiply(0.1)$rename("tmin_terraclim"))
  scaled <- scaled$addBands(img$select("pdsi")$multiply(0.01)$rename("pdsi_terraclim"))
  basin_mean(scaled, 4638, month_of(img))
}))
tc_df <- fetch(tc_monthly, c("pr_terraclim", "pet_terraclim",
                             "tmax_terraclim", "tmin_terraclim", "pdsi_terraclim"))

# ---- ERA5-Land monthly: precipitation in m, air temperature in K ------------
era <- ee$ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR")$filterDate(start, end)
era_monthly <- era$map(ee_utils_pyfunc(function(img) {
  img <- ee$Image(img)
  b <- img$select("total_precipitation_sum")$multiply(1000)$rename("pr_era5")
  b <- b$addBands(img$select("temperature_2m")$subtract(273.15)$rename("tmean_era5"))
  basin_mean(b, 11132, month_of(img))
}))
era_df <- fetch(era_monthly, c("pr_era5", "tmean_era5"))

# ---- SPEI: 12-month accumulation, the wet/dry state of the preceding year ---
spei <- ee$ImageCollection("CSIC/SPEI/2_11")$select("SPEI_12_month")$filterDate(start, end)
spei_monthly <- spei$map(ee_utils_pyfunc(function(img) {
  img <- ee$Image(img)
  basin_mean(img$rename("spei12"), 55660, month_of(img))
}))
spei_df <- fetch(spei_monthly, "spei12")

# ---- Join, derive, write ----------------------------------------------------
monthly <- Reduce(function(a, b) full_join(a, b, by = "date"),
                  list(chirps_df, tc_df, era_df, spei_df)) |>
  mutate(year  = as.integer(substr(date, 1, 4)),
         month = as.integer(substr(date, 6, 7)),
         tmean_terraclim = (tmax_terraclim + tmin_terraclim) / 2,
         # Ensemble precipitation: mean of the three products, which is a plain
         # average and not the triple-collocation weighting of Arash (2024).
         pr_ensemble    = rowMeans(cbind(pr_chirps, pr_terraclim, pr_era5), na.rm = TRUE),
         tmean_ensemble = rowMeans(cbind(tmean_terraclim, tmean_era5), na.rm = TRUE),
         # Climatic water balance: the endorheic basin loses water only upward.
         p_minus_pet    = pr_ensemble - pet_terraclim) |>
  arrange(year, month) |>
  select(date, year, month, pr_chirps, pr_terraclim, pr_era5, pr_ensemble,
         pet_terraclim, p_minus_pet, tmax_terraclim, tmin_terraclim,
         tmean_terraclim, tmean_era5, tmean_ensemble, pdsi_terraclim, spei12)

annual <- monthly |>
  group_by(year) |>
  summarise(
    pr_chirps      = sum(pr_chirps, na.rm = TRUE),
    pr_terraclim   = sum(pr_terraclim, na.rm = TRUE),
    pr_era5        = sum(pr_era5, na.rm = TRUE),
    pr_ensemble    = sum(pr_ensemble, na.rm = TRUE),
    pet_terraclim  = sum(pet_terraclim, na.rm = TRUE),
    p_minus_pet    = sum(p_minus_pet, na.rm = TRUE),
    tmax_terraclim = mean(tmax_terraclim, na.rm = TRUE),
    tmin_terraclim = mean(tmin_terraclim, na.rm = TRUE),
    tmean_ensemble = mean(tmean_ensemble, na.rm = TRUE),
    pdsi_terraclim = mean(pdsi_terraclim, na.rm = TRUE),
    # SPEI-12 in October closes the hydrological year that feeds the recession.
    spei12_oct     = spei12[match(10L, month)],
    n_months       = sum(!is.na(pr_ensemble)),
    .groups = "drop")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(monthly, file.path(OUT_DIR, "climate_monthly.csv"), row.names = FALSE)
write.csv(annual,  file.path(OUT_DIR, "climate_annual.csv"),  row.names = FALSE)

cat(sprintf("climate_monthly.csv  %d rows  %s to %s\n",
            nrow(monthly), min(monthly$date), max(monthly$date)))
cat(sprintf("climate_annual.csv   %d rows  %d to %d\n",
            nrow(annual), min(annual$year), max(annual$year)))
cat(sprintf("mean annual precipitation (ensemble) %.0f mm; mean PET %.0f mm; mean air temperature %.1f C\n",
            mean(annual$pr_ensemble), mean(annual$pet_terraclim), mean(annual$tmean_ensemble)))
