library(terra)
library(climateR)
library(glue)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

site_name <- args[1]

site_data_path <- file.path(glue("../data/input/{site_name}/"))

site_data <- read_csv("sites.csv") %>% filter(site == site_name)

lon <- site_data$lon
lat <- site_data$lat

point <- data.frame(lon = lon, lat = lat) %>%
  vect(geom = c("lon", "lat"), crs = "EPSG:4326")

gridmet_vars <- c("tmmn", "tmmx", "pr")
gridmet_data <- getGridMET(point,
  varname = gridmet_vars,
  startDate = ymd("1979-01-01"),
  endDate = ymd("2024-12-31"),
  verbose = TRUE
)

models <- c(
  "CanESM2",
  "HadGEM2-CC365",
  "MRI-CGCM3"
)
scenarios <- c("rcp45", "rcp85")
maca_vars <- c("tasmin", "tasmax", "pr")
maca_data <- getMACA(
  point,
  varname = maca_vars,
  timeRes = "day",
  model = models,
  scenario = scenarios,
  startDate = ymd("2006-01-01"),
  endDate = ymd("2099-12-31"),
  verbose = TRUE
)

write_csv(gridmet_data, file = file.path(site_data_path, glue("gridmet_1979_2023.csv")))
write_csv(maca_data, file = file.path(site_data_path, glue("macav2metdata_2006_2099.csv")))
