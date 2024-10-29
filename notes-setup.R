requirements <- c(
  "sf", "terra", "sits", "tmap", "RStoolbox", "mapview", "mapedit",
  "raster", "tidyverse", "rasterVis", "jsonlite", "parallel",
  "stringr", "httr", "gdalUtilities", "maptiles", "rgeos", "sp")

setup <- function(requirements){
  missing.packages <- requirements[!(requirements %in% installed.packages()[,"Package"])];
  if(length(missing.packages)) {install.packages(
    missing.packages, repos = "https://cloud.r-project.org/"); }
  for(package_name in requirements){library(
    package_name,character.only=TRUE,quietly=TRUE);
  }
}

install.packages("~/Repos/configs/rgdal_1.6-6.tar.gz", 
                 repos = NULL, 
                 type = "source", 
                 configure.args = c(
                   "--with-proj-include=/usr/local/include",
                   "--with-proj-lib=/usr/local/lib"))

setup(requirements)

install.packages("rgdal", type = "source",
                 configure.args = c("--with-proj-include=/usr/local/include",
                                    "--with-proj-lib=/usr/local/lib"))

install.packages("~/Repos/configs/rgdal_1.6-6.tar.gz", 
                 repos = NULL, 
                 type = "source", 
                 configure.args = c(
                   "--with-proj-include=/usr/local/include",
                   "--with-proj-lib=/usr/local/lib"))


library(gdalcubes)

gdalcubes_options(parallel=8)

files = list.files("L8_Amazon", recursive = TRUE, 
                   full.names = TRUE, pattern = ".tif") 
length(files)