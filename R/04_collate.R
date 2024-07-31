library(terra)
library(glue)
library(hash)
library(reticulate)
library(tidyverse)
library(parallel)

args = commandArgs(trailingOnly=TRUE)

site <- args[1]
period <- args[2]

terraOptions(verbose = TRUE,
             memfrac = 0.9)

### Output from start_wb_v_1_5.py script is a bunch of npz files,
### one for each day/year/model/scenario.  Import python numpy library
### to read in these files as arrays then work with r spatial libraries
### to georeference and collate them into netCDFs.
use_condaenv(condaenv = "nps-wb", conda = "auto", required = NULL)
np <- import("numpy")

var_units <- hash(
  ##' soil_water' = "mm * 10",
  ##' PET' = "mm * 10",
  "AET" = "mm * 10",
  "Deficit" = "mm * 10"
  ##' runoff' = "mm * 10",
  ##' agdd' = "GDD",
  ##' accumswe' = "mm * 10",
  ##' rain' = "mm * 10"
)


if (period == "historical") {
  years <- 1979:2022 ## 2023 data is incomplete
  models <- c("historical")
  scenarios <- c("gridmet")
} else {
  years <- 2006:2099
  models <- c(
    ## "bcc-csm1-1-m",
    ## "bcc-csm1-1",
    ## "BNU-ESM",
    "CanESM2",
    ## "CNRM-CM5",
    ## "CSIRO-Mk3-6-0",
    ## "GFDL-ESM2G",
    ## "GFDL-ESM2M",
    "HadGEM2-CC365",
    ## "HadGEM2-ES365",
    ## "inmcm4",
    ## "IPSL-CM5A-LR",
    ## "IPSL-CM5A-MR",
    ## "IPSL-CM5B-LR",
    ## "MIROC5",
    ## "MIROC-ESM-CHEM",
    ## "MIROC-ESM",
    "MRI-CGCM3"
    ## "NorESM1-M"
  )
  scenarios <- c("rcp85", "rcp45")
}


make_spatraster <- function(f, var, year, reference) {
  ## Takes a filename f for a npz file generated from start_wb_v_1_5.py and
  ## creates a SpatRaster with crs and extent set from the reference
  ## and date properly set
  ## npz files are daily arrays of water balance outputs

  yday <- str_split_i(f, pattern = "/", i = -1) %>%  ## Pick off base file name after directory /s
    str_split_i(pattern = "_", i = 4) %>% ## Yday should be 4th element in filenames in output from collate script
    as.numeric() + 1 ## Ydays are indexed at 0
  date <- as_date(glue("{year}-{yday}"), format = "%Y-%j")

  npz <- np$load(f)

  new_rast <- rast(npz$f[["param"]],
                   crs = crs(reference),
                   extent = ext(reference)
                   )
  time(new_rast) <- date
  names(new_rast) <- var
  units(new_rast) <- var_units[[var]]

  return(new_rast)
}


make_collation <- function(options) {
  ## Collate together npz files (daily arrays of wb model outputs) in input dir
  ## Convert each npz into a spatraster, then append to output raster
  var <- options[1]
  year <- options[2]
  model <- options[3]
  scenario <- options[4]
  site <- options[5]

  script_data_dir <- file.path(glue("../data/input/{site}/"))
  input_data_dir <- file.path(glue("~/out/{site}/wb/"))
  output_data_dir <- file.path(glue("~/out/{site}/collated/"))
  reference <- rast(file.path(script_data_dir, "1980_dayl_resampled.nc4"))

  out_file <- file.path(output_data_dir, glue("{model}_{scenario}_{var}_{year}.nc"))
  print(out_file)

  ## Fix missing files
  ## if (file.exists(out_file)) {
  ##   print(glue("File exists {out_file}"))
  ##   return(2)
  ## }

  output_rast <- rast(
    nrows = nrow(reference),
    ncols = ncol(reference),
    crs = crs(reference),
    extent = ext(reference),
    resolution = res(reference)    
  )

  in_files <- list.files(path = input_data_dir, pattern = glue("{model}_{scenario}_{year}_.*_{var}.npz"), full.names = TRUE) |> str_sort(numeric = TRUE)

  for (f in in_files) {
    print(f)
    new_rast <- make_spatraster(f, var, year, reference)
    add(output_rast) <- new_rast
    ## file.remove(f)
  }

  ## Fix NAs now so I don't have to later
  output_rast <- subst(output_rast, -9999, NA)
  
  writeCDF(output_rast, filename = out_file, overwrite = TRUE, compression = 3)

  ## rm(output_rast)
  return(1)
}


options <- expand.grid(
  var = keys(var_units),
  year = years,
  model = models,
  scenario = scenarios,
  site = site
) %>%
  t() %>%
  data.frame()

mclapply(options,
         FUN = make_collation,
         mc.cores = 1
         )
