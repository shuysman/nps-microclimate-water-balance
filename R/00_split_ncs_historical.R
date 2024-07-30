library(tidyverse)
library(glue)
library(terra)
library(parallel)

terraOptions(verbose = TRUE)

site <- "holly_lake_small"

reference <- rast(glue("../data/input/{site}/dem/dem_nad83.tif"))

in_dir <- file.path("/media/smithers/shuysman/data/MACA/gye/historical/daily/")
## in_dir <- file.path("~/data/MACA/gye/forecasts/daily/")
out_dir <- file.path(glue("~/out/{site}/daily-split/"))
## out_dir <- file.path("~/out/daily-split/")

if (!file.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

process_gcm <- function(options) {
  variable <- options[1]
  gcm <- options[2]
  scenario <- options[3]

  filename <- file.path(in_dir,
                        glue("{variable}_1979_CurrentYear_daily_gye.nc"))

  r <- rast(filename)

  for (n in 1:nlyr(r)) {
    lyr <- subset(r, n)

    first_year <- 1979
    second_year <- 2023

    yday <- yday(time(lyr))
    year <- year(time(lyr))

    lyr <- project(lyr, reference, method = "near")
    ##lyr <- crop(lyr, reference)
    lyr <- resample(lyr, reference, method = "near")
    lyr <- crop(lyr, reference)

    writeRaster(lyr,
                filename = file.path(out_dir,
                                     glue("{n}_macav2metdata_{variable}_{gcm}_r1i1p1_{scenario}_{first_year}_{second_year}_GYE_daily_reprojected_with_extent{yday}_resampled.tif")),
                overwrite = TRUE)

    rm(lyr)
    gc()
  }
}

models <- c("historical")
scenarios <- c("gridmet")
## variables <- c("tmmx", "tmmn", "pr")
variables <- c("tmmx", "tmmn", "pr")

options <- expand.grid(variables = variables, model = models, scenario = scenarios) %>%
  t() %>%
  data.frame()  ### How to use (mc)lapply with expand.grid https://gist.github.com/ihrke/c54e8d7c82cd91ae51d1d62f8d29b936

mclapply(options,
         FUN = process_gcm,
         mc.cores = parallel::detectCores() - 2)

#MACA File Format = {YEAR_CHUNK}_macav2metdata_{PARAM}_{MODEL_PART1}_{MODEL_PART2}_{SCENARIO}_{FIRST_YEAR_OF_CHUNK}_{LAST_YEAR_OF_CHUNK}_CONUS_dail_reprojected_with_extent{DAYNUMBER}_resampled.tif


## models <- c("inmcm4")
## scenarios <- c("rcp85")
## variables <- c("tmmx")
