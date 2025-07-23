library(terra)
library(climateR)
library(glue)
library(tidyverse)

# https://gist.github.com/debruine/e096f8142d3b383dc55f195bba5c7b0a
#' Check readline input
#'
#' @param prompt the prompt for readline
#' @param type what type of check to perform, one of c("numeric", "character", "length", "minlength", "maxlength", "exact", "grep")
#' @param compare the comparator for exact and (min|max)length types
#' @param ... other arguments to pass to grep
#'
#' @return the validated result of readline
#' @export
#'
#' @examples
#' readline_check("Type a number: ", "numeric")
#' readline_check("Type two characters: ", "length", 2)
#' readline_check("Type at least 3 characters: ", "minlength", 3)
#' readline_check("Type no more than 4 characters: ", "maxlength", 4)
#' readline_check("Type a letter and a number: ", "grep", pattern = "^[a-zA-Z]\\d$")
#'
readline_check <- function(prompt, type = c("numeric", "length", "minlength", "maxlength", "exact", "grep"), compare = NULL, ...) {
  input <- readline(prompt)
  if (type == "numeric") {
    check <- suppressWarnings(!is.na(as.numeric(input)))
  } else if (type == "exact") {
    check <- input == compare
  } else if (type == "length") {
    check <- nchar(input) == compare
  } else if (type == "minlength") {
    check <- nchar(input) >= compare
  } else if (type == "maxlength") {
    check <- nchar(input) <= compare
  } else if (type == "grep") {
    check <- grep(x = input, ...) %>% length() > 0
  } else {
    check <- FALSE # default false if type is wrong?
  }
  if (!check) {
    readline_check(prompt, type[1], compare, ...)
  } else {
    input
  }
}


metdata_elev <- rast("../data/metdata_elevationdata.nc")

if (interactive()) {
  site_name <- readline_check("Enter site name. It should be in machine readable format (e.g., 'Avalanche Peak' -> 'avalanche_peak': ", "maxlength", 16)
  site_data_path <- file.path(glue("../data/input/{site_name}"))

  dem <- rast(file.path(site_data_path, "dem/dem_nad83.tif"))

  shape_files <- list.files(file.path(site_data_path, "shapefile"), pattern = "*.shp", full.names = TRUE)
  len_shape_files <- length(shape_files)
  if (len_shape_files > 1) {
    abort(glue("Found {len_shape_files} shapefiles. Only one shapefile can be in data/input/{name}/shapefile. Please remove extra files. {shape_files}"))
  } else {
    message(glue("Found len_shape_files} shapefile. {shape_files}"))
  }

  site_poly <- vect(shape_files) %>% project(dem)

  mean_elevation <- terra::extract(dem, site_poly, fun = "mean")$USGS_1m

  message(glue("Site poly mean elevation in USGS 1m DEM: {mean_elevation}"))

  elev_data <- terra::extract(metdata_elev, site_poly, xy = TRUE)
  message(glue("metdata elevation: "))
  print(elev_data)

  nrow_elev_data <- nrow(elev_data)
  if (nrow_elev_data > 1) {
    abort("TODO: implement method to select metdata gridcell pixel with close elevation")
  }

  site_elev <- elev_data$elevation

  centroid <- centroids(project(site_poly, crs("EPSG:4269")), inside = TRUE) %>% geom(df = TRUE)
  centroid <- centroid[1, ]
  lat <- centroid$y
  lon <- centroid$x

  ## ggplot() +
  ##   geom_spatraster(data = crop(metdata_elev, terra::buffer(site_poly, 7000))) +
  ##   geom_sf(data = site_poly)

  site_data <- data.frame(site = site_name, lon = lon, lat = lat, metdata_elev = site_elev)
  message(glue("Site data: \n"))
  print(site_data)

  sites_file <- file.path("sites.csv")

  if (!file.exists(sites_file)) {
    append <- readline_check("sites.csv not found. Write to new sites.csv? (y/n): ", type = "grep", pattern = "[yn]")
    if (append == "y") {
      write_csv(site_data, file = sites_file, append = FALSE, col_names = TRUE)
    }
  } else {
    append <- readline_check("sites.csv found. Append to sites.csv? (y/n): ", type = "grep", pattern = "[yn]")
    if (append == "y") {
      write_csv(site_data, file = sites_file, append = TRUE)
    }
  }
} else {
  ### Batch Mode
  ### Assumes that polygons are checked for overlap with raster grid cells
  args <- commandArgs(trailingOnly = TRUE)

  site_name <- args[1]

  site_data_path <- file.path(glue("../data/input/{site_name}/"))

  site_data <- read_csv("sites.csv") %>% filter(site == site_name)

  lon <- site_data$lon
  lat <- site_data$lat
}

point <- data.frame(lon = lon, lat = lat) %>%
  vect(geom = c("lon", "lat"), crs = "EPSG:4326")

message("Fetching climate data... This may take a while.")

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
