#!/usr/bin/env Rscript
# =============================================================================
# Build ground-truth training/validation points for the Lake Chilwa classifier.
# Joins a differential-GPS point export to cover-class labels, validates the
# points against the self-derived basin, and writes
#   03.outputs/SHP/ground_truth_points.shp  (+ 03.outputs/ground_truth_points.csv)
# which the manuscript's 2.1.1 `ground-truth-assembly` chunk and Sections 2.2.G /
# 2.2.H consume. Built 2026-07-06.
#
# WHY: the field photographs carry a date but no GPS coordinate (0 of 245), so the
# location comes from the DGPS survey and participatory mapping, and the class
# label comes from the photographic/field evidence. This script performs that join.
#
# INPUTS (edit `cfg`):
#   cfg$dgps   DGPS export. A CSV with longitude/latitude columns, OR a points
#              .gpx / .shp / .geojson / .kml. Coordinate, id, date, and (if
#              present) class columns are auto-detected from common names.
#   cfg$labels OPTIONAL separate CSV mapping a join key to a class, used only when
#              the DGPS export has no class/label column of its own.
#   cfg$key    the column common to the DGPS export and the labels file
#              (e.g. "site_id", "waypoint", "photo").
#
# OUTPUT fields: class (1 water, 2 flooded veg, 3 dry veg, 4 bare soil),
#   class_name, site_id, date, photo, source, lon, lat.
#
# If cfg$dgps does not exist, the script runs a SELF-TEST on synthetic in-basin
# points so you can see the schema and confirm it runs before your export is ready.
#
# RUN: Rscript 05.scripts/build_ground_truth_points.R
# =============================================================================

suppressMessages({library(sf); library(dplyr); library(stringr)})
sf::sf_use_s2(FALSE)

cfg <- list(
  dgps    = here::here("02.inputs", "GPS", "dgps_points.csv"),   # <-- your export
  labels  = here::here("02.inputs", "GPS", "point_labels.csv"),  # <-- optional
  key     = "site_id",
  basin   = here::here("03.outputs", "SHP", "chilwa_basin.shp"),
  out_shp = here::here("03.outputs", "SHP", "ground_truth_points.shp"),
  out_csv = here::here("03.outputs", "ground_truth_points.csv"),
  crs     = 4326
)

# ---- class dictionary: the four manuscript classes ---------------------------
CLASS_NAME <- c("1" = "open water", "2" = "flooded vegetation",
                "3" = "dry vegetation", "4" = "bare soil")

# Map free text or a 1-4 code to the integer class. Extend the patterns freely.
label_to_code <- function(x) {
  s <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    str_detect(s, "^[1-4]$")                                  ~ suppressWarnings(as.integer(s)),
    str_detect(s, "open.?water|^water$|lake|channel|lagoon|pool") ~ 1L,
    str_detect(s, "flood|emerg|typha|marsh|swamp|reed|macrophyte") ~ 2L,
    str_detect(s, "dry.?veg|grass|crop|shrub|woodland|savanna|vegetation") ~ 3L,
    str_detect(s, "bare|soil|lake.?bed|sand|mud|exposed|salt|substrate") ~ 4L,
    TRUE ~ NA_integer_)
}

# ---- column auto-detection ---------------------------------------------------
pick <- function(nms, cands) {
  hit <- nms[tolower(nms) %in% cands]; if (length(hit)) hit[1] else NA_character_
}
LON <- c("lon","long","longitude","x","point_x","gpslongitude","lng","easting")
LAT <- c("lat","latitude","y","point_y","gpslatitude","northing")
IDC <- c("site_id","id","name","waypoint","wpt","ident","point","photo")
DTC <- c("date","datetime","time","timestamp")
CLC <- c("class","label","cover","type","landcover","land_cover","classname","class_name")

read_points <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("gpx","shp","geojson","json","kml")) {
    g  <- sf::st_read(path, quiet = TRUE) |> sf::st_zm(drop = TRUE) |> sf::st_transform(cfg$crs)
    co <- sf::st_coordinates(g)
    df <- sf::st_drop_geometry(g); df$lon <- co[, 1]; df$lat <- co[, 2]; df
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

# ---- load or synthesise the DGPS export --------------------------------------
demo <- !file.exists(cfg$dgps)
if (demo) {
  message("No DGPS export at ", cfg$dgps, " -- running SELF-TEST on synthetic points.")
  bb <- sf::st_bbox(sf::st_read(cfg$basin, quiet = TRUE))
  set.seed(1); n <- 80
  dgps <- data.frame(
    site_id = sprintf("WP%03d", seq_len(n)),
    lon  = runif(n, bb["xmin"], bb["xmax"]),
    lat  = runif(n, bb["ymin"], bb["ymax"]),
    date = as.character(as.Date("2013-01-01") + sample(0:400, n, TRUE)),
    cover = sample(c("open water", "flooded Typha", "dry grassland", "exposed lakebed"),
                   n, TRUE),
    photo = sprintf("CIMG%04d.JPG", sample(1000:2000, n)),
    stringsAsFactors = FALSE)
} else {
  dgps <- read_points(cfg$dgps)
}

nms   <- names(dgps)
lon_c <- pick(nms, LON); lat_c <- pick(nms, LAT)
id_c  <- pick(nms, IDC); dt_c  <- pick(nms, DTC); cl_c <- pick(nms, CLC)
ph_c  <- if ("photo" %in% tolower(nms)) nms[tolower(nms) == "photo"][1] else NA_character_
if (is.na(lon_c) || is.na(lat_c))
  stop("Could not find longitude/latitude columns. Columns present: ",
       paste(nms, collapse = ", "))

gt <- data.frame(
  site_id   = if (!is.na(id_c)) as.character(dgps[[id_c]]) else sprintf("PT%04d", seq_len(nrow(dgps))),
  lon       = suppressWarnings(as.numeric(dgps[[lon_c]])),
  lat       = suppressWarnings(as.numeric(dgps[[lat_c]])),
  date      = if (!is.na(dt_c)) as.character(dgps[[dt_c]]) else NA_character_,
  photo     = if (!is.na(ph_c)) as.character(dgps[[ph_c]]) else NA_character_,
  raw_class = if (!is.na(cl_c)) as.character(dgps[[cl_c]]) else NA_character_,
  stringsAsFactors = FALSE)

# ---- attach labels from a separate file if the export carried none ------------
if (all(is.na(gt$raw_class)) && file.exists(cfg$labels)) {
  lab  <- read.csv(cfg$labels, stringsAsFactors = FALSE, check.names = FALSE)
  lnms <- names(lab)
  lkey <- pick(lnms, tolower(cfg$key)); lcl <- pick(lnms, CLC)
  if (is.na(lkey) || is.na(lcl))
    stop("labels file needs a '", cfg$key, "' key column and a class column.")
  gt$raw_class <- lab[[lcl]][match(gt$site_id, as.character(lab[[lkey]]))]
}

gt$class      <- label_to_code(gt$raw_class)
gt$class_name <- unname(CLASS_NAME[as.character(gt$class)])
gt$source     <- if (demo) "synthetic-selftest" else basename(cfg$dgps)

# ---- validate: finite coords, known class, inside the basin -------------------
n_in        <- nrow(gt)
n_bad_coord <- sum(!(is.finite(gt$lon) & is.finite(gt$lat)))
n_no_class  <- sum(is.na(gt$class))
gt <- gt[is.finite(gt$lon) & is.finite(gt$lat) & !is.na(gt$class), ]

basin  <- sf::st_read(cfg$basin, quiet = TRUE) |> sf::st_transform(cfg$crs)
pts    <- sf::st_as_sf(gt, coords = c("lon", "lat"), crs = cfg$crs, remove = FALSE)
inside <- lengths(sf::st_intersects(pts, sf::st_union(basin))) > 0
n_outside <- sum(!inside)
pts <- pts[inside, ]

# ---- write -------------------------------------------------------------------
dir.create(dirname(cfg$out_shp), recursive = TRUE, showWarnings = FALSE)
keep <- c("class", "class_name", "site_id", "date", "photo", "source", "lon", "lat")
sf::st_write(pts[, keep], cfg$out_shp, delete_dsn = TRUE, quiet = TRUE)
write.csv(sf::st_drop_geometry(pts[, keep]), cfg$out_csv, row.names = FALSE)

# ---- report ------------------------------------------------------------------
cat(sprintf("Ground-truth points: %d written of %d input\n", nrow(pts), n_in))
cat(sprintf("  dropped: %d bad coordinate, %d unlabelled/unmatched, %d outside basin\n",
            n_bad_coord, n_no_class, n_outside))
cat("Class balance:\n"); print(table(factor(pts$class_name, levels = CLASS_NAME)))
cat(sprintf("Wrote: %s\n       %s\n", cfg$out_shp, cfg$out_csv))
if (demo) cat("\nSELF-TEST mode. Point cfg$dgps at your real export and re-run.\n")
