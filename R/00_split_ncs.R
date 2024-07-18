library(tidyverse)
library(glue)
library(terra)
library(parallel)

terraOptions(verbose = TRUE)

## r <- rast("/media/smithers/shuysman/data/MACA/gye/forecasts/daily/pr_BNU-ESM_rcp45_2006-2099_daily_gye.nc")

## reference <- rast("/home/steve/OneDrive/burroughs_wb/data/burroughs_creek_USGS1m_clipped_nad83.tif")
reference <- rast("../data/StephenHuysman_GRTE_WBP_ModelingAreas/surprise/dem/surprise_USGS1m_nad83.tif")

##in_dir <- file.path("/media/smithers/shuysman/data/MACA/gye/forecasts/daily/")
in_dir <- file.path("~/data/MACA/gye/forecasts/daily/")
##out_dir <- file.path("/media/smithers/shuysman/data/MACA/gye/forecasts/daily-split/")
out_dir <- file.path("~/out/surprise/daily-split/")
##out_dir <- file.path("/tmp/test/")

process_gcm <- function(options) {
    variable <- options[1]
    gcm <- options[2]
    scenario <- options[3]

    filename <- file.path(in_dir,
                          glue("{variable}_{gcm}_{scenario}_2006-2099_daily_gye.nc"))

    r <- rast(filename)
    
    for (n in 1:nlyr(r)) {
        lyr <- subset(r, n)

        first_year <- str_split_i(sources(lyr)[1], i = 4, pattern = '_') %>% str_split_i(i = 1, pattern = "-")
        second_year <- str_split_i(sources(lyr)[1], i = 4, pattern = '_') %>% str_split_i(i = 2, pattern = "-")
        
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
    }
}

models <- c(
    ## "bcc-csm1-1-m",
    ## "bcc-csm1-1",
    ## "BNU-ESM",
    "CanESM2", ## Hot Wet
    ## "CNRM-CM5",
    ## "CSIRO-Mk3-6-0",
    ## "GFDL-ESM2G",
    ## "GFDL-ESM2M",
    "HadGEM2-CC365", ## Hot Dry
    ## "HadGEM2-ES365",
    ## "inmcm4",
    ## "IPSL-CM5A-LR",
    ## "IPSL-CM5A-MR",
    ## "IPSL-CM5B-LR",
    ## "MIROC5",
    ## "MIROC-ESM-CHEM",
    ## "MIROC-ESM",
    "MRI-CGCM3" ## Warm Dry/ Warm Wet
    ## "NorESM1-M"
)
### CCSM4 tmmn data missing
scenarios <- c('rcp85', 'rcp45')
## variables <- c("tmmx", "tmmn", "pr")
variables <- c("tmmx", "tmmn", "pr")

options <- expand.grid(variables = variables, model = models, scenario = scenarios) %>%
    t() %>%
    data.frame()  ### How to use (mc)lapply with expand.grid https://gist.github.com/ihrke/c54e8d7c82cd91ae51d1d62f8d29b936

mclapply(options,
         FUN = process_gcm,
         mc.cores = parallel::detectCores())

    #MACA File Format = {YEAR_CHUNK}_macav2metdata_{PARAM}_{MODEL_PART1}_{MODEL_PART2}_{SCENARIO}_{FIRST_YEAR_OF_CHUNK}_{LAST_YEAR_OF_CHUNK}_CONUS_dail_reprojected_with_extent{DAYNUMBER}_resampled.tif


## models <- c("inmcm4")
## scenarios <- c("rcp85")
## variables <- c("tmmx")
