library(terra)
library(glue)
library(hash)
library(reticulate)
library(tidyverse)
library(parallel)

terraOptions(verbose = TRUE)

### Output from start_wb_v_1_5.py script is a bunch of npz files,
### one for each day/year/model/scenario.  Import python numpy library
### to read in these files as arrays then work with r spatial libraries
### to georeference and collate them into netCDFs.
use_condaenv(condaenv = "nps-wb", conda = "auto", required = NULL)
np <- import("numpy")

script_data_dir <- file.path("../data/StephenHuysman_GRTE_WBP_ModelingAreas/static_basin/")
input_data_dir <- file.path("~/out/static_east/wb")
##input_data_dir <- file.path("/media/smithers/shuysman/data/out/nps-wb/avalanche/historical_025/")
output_data_dir <- file.path("~/out/static_east/collated/")
##output_data_dir <- file.path("/media/smithers/shuysman/data/out/nps-wb/avalanche/historical_025/collated/")
reference <- rast(file.path(script_data_dir, "dem/static_east_dem_clipped_nad83.tif"))

historical_years <- 1979:2022 ## 2023 data is incomplete
## historical_years <- 2022
projection_years <- 2006:2099

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

historical_models <- c("historical")
historical_scenarios <- c("gridmet")

projection_models <- c(
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
projection_scenarios <- c("rcp85", "rcp45")


make_spatraster <- function(f, var, year) {
  ## Takes a filename f for a npz file generated from start_wb_v_1_5.py and
  ## creates a SpatRaster with crs and extent set from the reference
  ## and date properly set
  ## npz files are daily arrays of water balance outputs

  yday <- str_split_i(f, pattern = "_", i = 4) %>% as.numeric() + 1
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

  out_file <- file.path(output_data_dir, glue("{model}_{scenario}_{var}_{year}.nc"))
  print(out_file)

  ## Fix missing files
  if (file.exists(out_file)) {
    print(glue("File exists {out_file}"))
    return(2)
  }

  output_rast <- rast(
    nrows = nrow(reference),
    ncols = ncol(reference),
    ## xmin = xmin(reference),
    ## xmax = xmax(reference),
    ## ymin = ymin(reference),
    ## ymax = ymax(reference),
    crs = crs(reference),
    extent = ext(reference),
    resolution = res(reference),
    )

  in_files <- str_subset(wb_files, pattern = glue("{model}_{scenario}_{year}_.*_{var}.npz")) |> str_sort(numeric = TRUE)

  for (f in in_files) {
    print(f)
    new_rast <- make_spatraster(f, var, year)
    add(output_rast) <- new_rast
    ## file.remove(f)
  }

  writeCDF(output_rast, filename = out_file, overwrite = TRUE, compression = 3)

  ## rm(output_rast)
  return(1)
}


wb_files <- list.files(input_data_dir, full.names = TRUE)

historical_options <- expand.grid(
  var = keys(var_units),
  year = historical_years,
  model = historical_models,
  scenario = historical_scenarios
) %>%
  t() %>%
  data.frame()

mclapply(historical_options,
         FUN = make_collation,
         mc.cores = 1
         )

## projection_options <- expand.grid(var = keys(var_units),
##                                   year = projection_years,
##                                   model = projection_models,
##                                   scenario = projection_scenarios) %>%
##     t() %>%
##     data.frame()

## mclapply(projection_options,
##          FUN = make_collation,
##          mc.cores = 4)
