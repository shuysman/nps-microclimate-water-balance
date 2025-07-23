library(optparse)
library(tidyverse)
library(tidyterra)
library(sf)
library(terra)
library(glue)
library(janitor)
library(FedData)

option_list <- list(
  make_option(c("-n", "--name"), type="character", callback=janitor::clean_names,
              help="Site name"),
  make_option(c("-s", "--shape_file"), type="character",
              help="Input shapefile"),
  make_option(c("-d", "--dem"), type="character", ### TODO: Maybe we can automate downloading the USGS 1m DEM?
              help="Input USGS 1m DEM (uncropped)"),
  make_option(c("-c", "--target_crs"), type="character", default=crs("EPSG:26912"),
              help="Target CRS for output [default EPSG:26912]")
)

## Testing parameters
name <- "test"
shape_file <- file.path("../data/input/test/shapefile/sample.shp")
dem <- file.path("../data/input/test/dem/USGS_1m.tif")
target_crs <- crs("EPSG:26912")

data_dir <- file.path(glue("../data/input/{name}/"))

## Coerce polygon to output CRS in case it is something else
site_poly <- vect(shape_file) %>% project(target_crs)

## Crop and write out 1m DEM
reference <- rast(dem) %>% crop(site_poly) %>% writeRaster(file.path(data_dir, "dem/dem_nad83.tif"), overwrite=TRUE)

## Create and write out slope, aspect, and hillshade
## Slope in radians, aspect in degrees
## Aspect is used to run fold_aspect() which takes aspects
## in degrees and returns folded aspect in radians
## Hillshade is only used for visualizations
slope <- reference %>%
  terrain(v = "slope", unit = "radians")
writeRaster(slope, file.path(data_dir, "dem/slope_nad83.tif"), overwrite=TRUE)

aspect <- reference %>%
  terrain(v="aspect",unit="degrees")
writeRaster(aspect, file.path(data_dir, "dem/aspect_nad83.tif"), overwrite=TRUE)
## We also need aspect in radians for terra::shade
aspect_rad <- reference %>%
  terrain(v="aspect",unit="radians")

hillshade <- shade(slope = slope, aspect = aspect_rad)
writeRaster(hillshade, file.path(data_dir, "dem/hillshade_nad83.tif"), overwrite=TRUE)


t50 <- rast(file.path("../data/merged_jennings2.tif"))
t50 <- project(subset(t50, 1), reference, method = "near")
t50 <- resample(t50, reference, method = "near")
t50 <- crop(t50, reference)

writeRaster(t50, file.path(data_dir, "jennings_t50_coefficients.tif"), overwrite=TRUE)


## Download and rasterize SSURGO AWS data
soils <- get_ssurgo(reference, label = name)

soils_aws <- soils$spatial %>%
  mutate(MUKEY = as.numeric(MUKEY)) %>%
  left_join(soils$tabular$muaggatt, by = join_by("MUKEY" == "mukey")) %>%
  mutate(aws025wta = aws025wta * 10) %>% ## cm to mm
  st_transform(target_crs)

ggplot(soils_aws) +
  geom_sf(aes(fill = aws025wta)) +
  geom_sf(data = site_poly, fill = NA, lwd = 2) +
  scale_fill_viridis_c(name = "AWS (0–25 cm)\n(mm water)", na.value = "grey80") +
  labs(
    title = "Soil Water Holding Capacity to 25 cm",
    caption = "Data: SSURGO via FedData"
  ) +
  theme_minimal()

soils_025 <- rasterize(vect(soils_aws), reference, fun = "max", field = "aws025wta")

writeRaster(soils_025, file.path(data_dir, "soil/soil_whc_025.tif"), overwrite=TRUE)
