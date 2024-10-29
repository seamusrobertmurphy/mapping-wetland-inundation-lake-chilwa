tiles = c("36LYH", "36LYJ", "36LZH", "36LZJ"),
aoi = read_sf("~/Repos/mapping-wetland-inundation-lake-chilwa/inputs/chilwa_watershed_4326.shp") 

cube <- sits_cube(
  source = "MPC",
  collection = "LANDSAT-C2-L2",
  bands = c("BLUE", "GREEN", "RED", "NIR08", "SWIR16", "SWIR22", "CLOUD"),
  start_date = "2012-07-01",
  end_date = "2014-07-01",
  roi = aoi
)

cube_local = sits_cube_copy(
  cube,
  roi = aoi,
  res = 30,
  n_tries = 5,
  output_dir = "/Volumes/TOSHIBA_EXT/chilwa/data/raw_cube/MPC",
  progress = T
)


"B01", "B02", "B03", "B04", "B05", "B06", "B07"





#Derive NDWI
cube_202407_spectral = sits_apply(
  data = cube_202407_spectral, 
  NDBR = (B08 - B12) / (B08 + B12), 
  output_dir = dir_reg_test,
  memsize = 6,
  multicores = 4,
  progress = T
)

#Derive NDMI
cube_202407_spectral = sits_apply(
  data = cube_202407_spectral, 
  NDMI = (B08 - B11) / (B08 + B11), 
  output_dir = dir_reg_test,
  memsize = 6,
  multicores = 8,
  progress = T
)

#Merging mosaics
ndvi = list.files(dir_reg_test, 
                  pattern = 'NDVI', full.names = T, all.files = FALSE)|>
  lapply(terra::rast)|>
  sprc() |>
  mosaic()
terra::mask(ndvi, vect(aoi))
aoi = sf::st_transform(aoi, crs(ndvi))
ndvi = terra::crop(ndvi, vect(aoi), mask=T)
ndvi = ndvi * 0.0001
writeRaster(ndvi, file.path(dir_reg_test, "NDVI_1994-08-01.tif", overwrite=T)
            
            
            
            
            
            cube_reg <- sits_cube(
              source = "MPC",
              collection = "LANDSAT-C2-L2",
              bands = c("BLUE", "GREEN", "RED", "NIR08", "SWIR16", "SWIR22", "NDVI"),
              start_date = "1994-07-01",
              end_date = "1996-07-01",
              roi = aoi,
              data_dir = dir_reg_test,
              parse_info = c("satellite", "sensor", "tile", "band", "date"),
              delim = "_"
            )
            
            
s1_vh = list.files("./cubes/2024_reg",pattern = 'VH', full.names = T) |>
  lapply(terra::rast) |>
  sprc() |>
  mosaic()
            