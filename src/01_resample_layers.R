library(tidyverse)
library(sf)
library(terra)
library(glue)

sites <- c("avalanche_peak", "chittenden", "cub_creek")

for (site in sites) {
  data_dir <- file.path(glue("../data/input/{site}/"))

  reference <- rast(file.path(data_dir, "dem/dem_nad83.tif"))


  ### I don't think 1980_dayl_na.nc4 is needed anymore. It was used with Penman calculations using Daymet climate data.
  ## dayl <- rast(file.path("../data/1980_dayl_na.nc4"))
  ## dayl2 <- project(dayl, reference, method = "near")
  ## dayl2 <- resample(dayl2, reference, method = "near")
  ## dayl2 <- crop(dayl2, reference)

  ## writeCDF(dayl2, file.path(data_dir, glue("1980_dayl_resampled.nc4")), compression = 9, overwrite = TRUE)

  t50 <- rast(file.path("../data/merged_jennings2.tif"))
  t50 <- project(subset(t50, 1), reference, method = "near")
  t50 <- resample(t50, reference, method = "near")
  t50 <- crop(t50, reference)

  writeRaster(t50, file.path(data_dir, "jennings_t50_coefficients.tif"), overwrite=TRUE)


  # soils <- st_read("/home/steve/OneDrive/core_areas/data/StephenHuysman_GRTE_WBP_ModelingAreas/soil/grte_modelingareas_ssurgo.gpkg") %>% vect

  soils <- st_read(file.path(data_dir, "soil/ssurgo.gpkg")) %>% vect()

  ### The curve number generator ssurgo download tool in QGIS returns strings for soil AWS.
  ### R loads these as factors, and when multiplied by 10, you get the index of the string * 10 instead of
  ### the AWS * 10.  This burned me and I lost a lot of time rerunning data.  Make sure that the
  ### soil AWS values are numeric and then multiply by ten to convert cm to mm.
  soils$aws025wta <- as.numeric(soils$aws025wta) * 10 
  soils$aws050wta <- as.numeric(soils$aws050wta) * 10
  soils$aws0100wta <- as.numeric(soils$aws0100wta) * 10
  soils$aws0150wta <- as.numeric(soils$aws0150wta) * 10


  soils_025 <- soils %>%
    rasterize(reference, field = "aws025wta") 

  writeRaster(soils_025, file.path(data_dir, "soil/soil_whc_025.tif"), overwrite=TRUE)

  soils_100 <- soils %>%
    rasterize(reference, field = "aws0100wta")

  writeRaster(soils_100, file.path(data_dir, "soil/soil_whc_100.tif"), overwrite=TRUE)
}
