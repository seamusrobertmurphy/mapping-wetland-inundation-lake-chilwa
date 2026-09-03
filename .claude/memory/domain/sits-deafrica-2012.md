# What sits and Digital Earth Africa add for 2012

2026-09-02. `sits` supplies no imagery that Earth Engine lacks, because its Landsat holdings are the
same USGS Collection 2 archive repackaged as cloud-optimised GeoTIFFs. Checked against the live
Digital Earth Africa STAC at `https://explorer.digitalearth.africa/stac/`, searching the basin
polygon from `03.outputs/JSON/chilwa_basin.geojson`. Matters because no route through `sits` produces
a February 2012 observation that does not exist in the archive.

2026-09-02. Going through `sits` gives less than going to Digital Earth Africa directly. The package
config at `inst/extdata/sources/config_source_deafrica.yml` in the clone at `/Volumes/PortableSSD/
Github/sits` fronts 15 collections, namely `ls5_sr`, `ls7_sr`, `ls8_sr`, `ls9_sr`, the Sentinel-1 and
Sentinel-2 collections, CHIRPS rainfall, an NDVI anomaly, ALOS PALSAR mosaic, Copernicus DEM, and the
`gm_ls8_ls9_annual`, `gm_s2_annual`, `gm_s2_rolling` and `gm_s2_semiannual` geomedians. The live
catalogue holds 106. Matters because the three collections worth having here are among the 91 it does
not front.

2026-09-02. The three worth having are `wofs_ls`, a per-scene Landsat water classification with 47
distinct dates over the basin in 2012 and a companion `wofs_ls_summary_annual`; `gm_ls5_ls7_annual`,
a ready-made cloud-free annual geomedian at 30 m that does include 2012; and `fc_ls`, fractional
cover. Matters most for WOfS, which is an independent water product from a different team and a
different algorithm on the same dates, so it is the external check the harmonic fit in
[[harmonic-fit-2012]] currently lacks.

2026-09-02. Earth Engine's `LANDSAT/LE07/C02/T1_L2` is not complete over this basin. Digital Earth
Africa lists three scenes it lacks, namely path 166 row 071 on 12 April 2012, and path 167 rows 070
and 071 on 15 December 2012. Measured by pixels rather than footprints, path 166 row 071 reaches 698
km2 of the basin, 8.0 per cent, while path 167 row 071 reaches 7,789 km2, 89.0 per cent, and row 070
reaches 2,762 km2. So the April scene adds little and the December pair adds a whole date outside the
fieldwork window. Matters because Earth Engine has been treated here as the archive of record.

2026-09-02. Digital Earth Africa's STAC geometries overstate coverage, so filter on pixels and never
on their footprints. A search on the true basin polygon returns path 168 row 070 for 22 February 2012,
`LE71680702012053ASN00` at 6 per cent cloud, which looks like an observation on the exact date of the
43-photograph session. It is not, because path 168 puts zero pixels inside the basin, confirmed by
masking rather than by geometry. Earth Engine's own `system:footprint` intersection is equally
unreliable here, returning zero for path 166 scenes that `filterBounds` had selected.

2026-09-02. The cost of moving the whole workflow off Earth Engine, measured rather than guessed.
Over the basin with cloud below 30 per cent the record is 1,303 Landsat scenes, being 240 from
Landsat 5 for 1984 to 2011, 627 from Landsat 7 for 1999 to 2024 and 436 from Landsat 8 for 2013 to
2024. The basin is 8,783 km2, so 9.8 million pixels at 30 m. Six bands as int16 gives 152.6 GB for
the full per-scene stack, 28.5 GB for the 2010 to 2014 harmonic window, and 4.8 GB for 41 annual
composites. Matters because the annual composite is where the cost falls off a cliff, so that is
where the boundary between Earth Engine and a local cube belongs.

2026-09-02. `sits` 1.5.4 is installed against R 4.4.1 and carries every modelling tool this study
needs, namely `sits_mixture_model` with an `rmse_band` argument for the RMS-residual check,
`sits_uncertainty`, `sits_classify`, `sits_smooth`, `sits_segment`, `sits_accuracy`, `sits_som` and
`sits_regularize`. Landsat 5, 7 and 8 reach it from four sources, being Digital Earth Africa, Microsoft
Planetary Computer, USGS and AWS, and MOD09A1 reaches it as `modis-09A1-061` on the Planetary
Computer. `sits_cube` also accepts `data_dir` with `parse_info`, so a cube can be built from local
GeoTIFFs with no API at all. Matters because archive access was assumed to be the blocker for
migration and it is not; only Envisat ASAR is genuinely absent.

2026-09-02. The manuscript at `01.manuscript/Manuscript_2026-08-03.qmd` holds 47 code chunks with 56
references to Earth Engine, and the Earth Engine work is reduction rather than modelling, being
`filterDate`, `filterBounds`, `map`, `median`, `reduceRegion` and `sampleRegions`, with only two calls
each to `unmix` and `classify`. Matters because the part that would move to `sits` is small and the
part that would have to be replaced by a 152 GB download is large.
