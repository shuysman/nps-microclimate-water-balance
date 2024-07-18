library(tidyverse)
library(sf)
library(terra)

data_dir <- file.path("../data/StephenHuysman_GRTE_WBP_ModelingAreas/surprise/")

reference <- rast(file.path(data_dir, "dem/surprise_USGS1m_nad83.tif"))

dayl <- rast(file.path(data_dir, "1980_dayl_na.nc4"))
dayl2 <- project(dayl, reference, method = "near")
dayl2 <- resample(dayl2, reference, method = "near")
dayl2 <- crop(dayl2, reference)

writeCDF(dayl2, file.path(data_dir, "1980_dayl_surprise.nc4"), compression = 9, overwrite=TRUE)

t50 <- rast(file.path(data_dir, "merged_jennings2.tif"))
t50 <- project(subset(t50, 1), reference, method = "near")
t50 <- resample(t50, reference, method = "near")
t50 <- crop(t50, reference)

writeRaster(t50, file.path(data_dir, "jennings_t50_coefficients_surprise.tif"), overwrite=TRUE)


soils <- st_read("/home/steve/OneDrive/core_areas/data/StephenHuysman_GRTE_WBP_ModelingAreas/soil/grte_modelingareas_ssurgo.gpkg") %>% vect

soils_025 <- soils %>%
    rasterize(reference, field = "aws025wta") * 10 ## cm to mm

writeRaster(soils_025, file.path(data_dir, "soil/soil_whc_025.tif"), overwrite=TRUE)

soils_100 <- soils %>%
    rasterize(reference, field = "aws0100wta") * 10 ## cm to mm

writeRaster(soils_100, file.path(data_dir, "soil/soil_whc_100.tif"), overwrite=TRUE)
