#!/usr/bin/env Rscript
# =============================================================================
# Spectral mixture analysis and endmember modelling for Lake Chilwa (Landsat)
# Reference implementation, literature-grounded. Built 2026-07-05.
#
# WHY THIS EXISTS / WHAT IT IMPROVES over tasks/02.d-methods-landsat.qmd:
#   1. Endmembers are IMAGE-DERIVED from pure pixels (threshold-guided, per
#      Halabisky et al. 2016), not hard-coded guesses. Hard-coded reflectance
#      vectors that do not match the scene bias every fraction.
#   2. The built/urban endmember is REMOVED. Chilwa settlements are Typha grass
#      and bamboo, sub-pixel at 30 m, spectrally vegetation, and collinear with
#      soil; a built endmember is poorly constrained (see recommendation memo).
#   3. Flooded vegetation is not forced as an optical endmember. In a 6-band
#      linear SMA it reads as a water+vegetation mixture; SAR (Sentinel-1
#      double-bounce) is the physically correct tool. Two schemes are provided:
#      SCHEME_3 (water, vegetation, soil) with SAR carrying flooded vegetation,
#      and SCHEME_4 (water, emergent veg, dry veg, soil) if you want the
#      wet-to-dry split optically (validate the flooded fraction carefully).
#   4. An RMS RESIDUAL image is computed and summarised, the standard SMA check
#      (Halabisky 2016): if the residual carries spatial structure inside the
#      target, an endmember is missing or wrong.
#   5. Unmixing is constrained (sum-to-one + non-negativity); a relaxed variant
#      (Halabisky) is noted for keeping pure-water / no-water extremes.
#
# Literature: Halabisky et al. 2016; Mattson 2023; Zhang et al. 2014;
# Zhang, Li & Zheng 2016; Zeng et al. 2014. Band limit: 6 Landsat optical bands
# allow ~<=6 endmembers for a determined solution; the stable region is 3-4.
#
# RUN: Rscript 05.scripts/sma_endmember_modelling.R   (needs an authenticated
# Earth Engine session; see gee_pipeline_operational_workplan.md). Cannot run
# here; you execute under your own credentials.
# =============================================================================

# ------------------------- CONFIG --------------------------------------------
cfg <- list(
  ee_project = Sys.getenv("EE_PROJECT", unset = "murphys-deforisk"),
  scheme     = "SCHEME_3",          # "SCHEME_3" or "SCHEME_4"
  ref_start  = "2020-06-01",         # reference composite window (dry season)
  ref_end    = "2020-10-31",
  bands      = c("blue","green","red","nir","swir1","swir2"),
  scale      = 30L,
  out_dir    = "03.outputs"
)
BANDS <- cfg$bands

suppressMessages({ library(reticulate); library(sf) })
ee <- reticulate::import("ee"); ee$Initialize(project = cfg$ee_project)
if (!dir.exists(cfg$out_dir)) dir.create(cfg$out_dir, recursive = TRUE)

# Assumes aoi_ee / aoi_sf exist (source _common.R first). Minimal fallback:
if (!exists("aoi_ee")) {
  poi <- ee$Geometry$Point(c(35.55, -15.35))
  aoi_ee <- ee$FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_6")$
    filterBounds(poi)$geometry()
}

# --------------------- reference composite + indices -------------------------
scale_l <- function(image) image$select("SR_B.*")$multiply(0.0000275)$add(-0.2)
# (Assumes a harmonised, scaled, cloud-masked 'landsat' collection is available
# from the main pipeline; here we build a quick L8 reference for endmember work.)
l8 <- ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(aoi_ee)$
  filterDate(cfg$ref_start, cfg$ref_end)$filter(ee$Filter$lt("CLOUD_COVER", 30))
harmonise_l8 <- function(image) {
  o <- image$select(c("SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7"),
                    BANDS)$multiply(0.0000275)$add(-0.2)
  o$copyProperties(image, list("system:time_start"))
}
ref <- l8$map(harmonise_l8)$median()$clip(aoi_ee)
mndwi <- ref$normalizedDifference(c("green","swir1"))$rename("MNDWI")
ndvi  <- ref$normalizedDifference(c("nir","red"))$rename("NDVI")
ref_i <- ref$addBands(mndwi)$addBands(ndvi)

# --------------------- image-derived endmembers ------------------------------
# Threshold-guided pure-pixel masks, then take the MEDIAN spectrum per class as
# the endmember (robust to outliers). Tune thresholds to your scene.
pure_water <- ref_i$updateMask(ref_i$select("MNDWI")$gt(0.2))
pure_veg   <- ref_i$updateMask(ref_i$select("NDVI")$gt(0.35)$
                               And(ref_i$select("MNDWI")$lt(0)))
pure_soil  <- ref_i$updateMask(ref_i$select("NDVI")$lt(0.12)$
                               And(ref_i$select("MNDWI")$lt(0)))
# For SCHEME_4, split vegetation into emergent (wetter, near swamp) and dry.
pure_emerg <- ref_i$updateMask(ref_i$select("NDVI")$gt(0.35)$
                               And(ref_i$select("MNDWI")$gt(-0.1))$
                               And(ref_i$select("MNDWI")$lt(0.2)))
pure_dry   <- ref_i$updateMask(ref_i$select("NDVI")$gt(0.35)$
                               And(ref_i$select("MNDWI")$lt(-0.1)))

median_spectrum <- function(masked) {
  v <- masked$select(BANDS)$reduceRegion(
    reducer = ee$Reducer$median(), geometry = aoi_ee,
    scale = cfg$scale, maxPixels = 1e9)$getInfo()
  as.numeric(v[BANDS])
}

if (cfg$scheme == "SCHEME_3") {
  em_names <- c("water_frac","veg_frac","soil_frac")
  em <- list(median_spectrum(pure_water),
             median_spectrum(pure_veg),
             median_spectrum(pure_soil))
} else {                                  # SCHEME_4
  em_names <- c("water_frac","emerg_veg_frac","dry_veg_frac","soil_frac")
  em <- list(median_spectrum(pure_water),
             median_spectrum(pure_emerg),
             median_spectrum(pure_dry),
             median_spectrum(pure_soil))
}
cat("Image-derived endmembers (", cfg$scheme, "):\n"); names(em) <- em_names
for (n in em_names) cat(sprintf("  %-16s %s\n", n, paste(round(em[[n]],4), collapse=", ")))

# Guard against near-collinear endmembers (unstable fractions)
em_mat <- do.call(rbind, em)
cc <- cor(t(em_mat))
if (any(abs(cc[lower.tri(cc)]) > 0.999))
  warning("Two endmembers are near-collinear; fractions may be unstable. Revisit thresholds or drop a class.")

# --------------------- constrained unmixing ----------------------------------
fractions <- ref$select(BANDS)$unmix(
  endmembers = lapply(em, as.list),
  sumToOne = TRUE, nonNegative = TRUE      # relax nonNegative for Halabisky-style
)$rename(em_names)

# --------------------- RMS residual (the SMA check) --------------------------
# Reconstruct each band from fractions x endmembers, then RMSE across bands.
recon <- ee$Image(0)
band_imgs <- lapply(seq_along(BANDS), function(bi) {
  acc <- ee$Image(0)
  for (k in seq_along(em_names))
    acc <- acc$add(fractions$select(em_names[k])$multiply(em[[k]][bi]))
  acc$rename(paste0("recon_", BANDS[bi]))
})
recon <- ee$Image$cat(band_imgs)
resid <- ee$Image$cat(lapply(seq_along(BANDS), function(bi)
  ref$select(BANDS[bi])$subtract(recon$select(paste0("recon_", BANDS[bi])))$
    pow(2)$rename(paste0("e2_", BANDS[bi]))))
rmse <- resid$reduce(ee$Reducer$mean())$sqrt()$rename("rmse")

rmse_stats <- rmse$reduceRegion(
  reducer = ee$Reducer$mean()$combine(ee$Reducer$stdDev(), sharedInputs = TRUE),
  geometry = aoi_ee, scale = cfg$scale, maxPixels = 1e9)$getInfo()
cat(sprintf("\nRMS residual: mean %.4f, sd %.4f\n",
            rmse_stats$rmse_mean, rmse_stats$rmse_stdDev))
cat("Interpretation: a low, spatially unstructured RMSE means the endmember set\n",
    "accounts for the scene. Map 'rmse' and inspect for structure inside the\n",
    "wetland; structure means a missing or wrong endmember (Halabisky 2016).\n")

# --------------------- basin-mean fractions + export -------------------------
frac_means <- fractions$reduceRegion(
  reducer = ee$Reducer$mean(), geometry = aoi_ee,
  scale = cfg$scale, maxPixels = 1e9)$getInfo()
cat("\nBasin-mean fractions:\n"); print(unlist(frac_means))
write.csv(data.frame(class = names(unlist(frac_means)),
                     fraction = as.numeric(unlist(frac_means))),
          file.path(cfg$out_dir, paste0("sma_endmembers_", cfg$scheme, ".csv")),
          row.names = FALSE)

# =============================================================================
# LOCAL R ALTERNATIVE (no Earth Engine), using RStoolbox multiple-endmember SMA.
# Use when you have a local Landsat GeoTIFF stack 'ras' with the 6 bands.
# -----------------------------------------------------------------------------
# library(RStoolbox); library(terra)
# ras <- rast("landsat_2020_6band.tif"); names(ras) <- BANDS
# em_df <- as.data.frame(em_mat); colnames(em_df) <- BANDS
# em_df$class <- em_names
# sma <- mesma(img = ras, em = em_df, method = "NNLS")  # non-negative LS
# # sma has one fraction layer per class plus an RMSE layer; check the RMSE:
# summary(sma[["RMSE"]])
# writeRaster(sma, file.path(cfg$out_dir, "sma_fractions_local.tif"), overwrite = TRUE)
# =============================================================================

cat("\nDone. Endmember scheme:", cfg$scheme,
    "| built/urban endmember removed | RMS residual reported.\n")
