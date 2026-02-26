NPS Gridded Water Balance Model - VM Setup and Workflow
=======================================================

This document describes how to run the water balance model on the virtual
machine using the 'nps_wb' conda environment.


!! IMPORTANT - OUTPUT UNITS !!
All water balance model outputs (AET and Deficit, both daily and annual
sums) are stored as mm * 10 (16-bit integer compression). You MUST divide
raster values by 10 to obtain actual millimeter values. For example, a
stored pixel value of 3500 = 350.0 mm. This applies to all NetCDF files
in output/{site}/wb/ and output/{site}/sums/.


PREREQUISITES
-------------
1. Conda/Mamba installed
2. The 'nps_wb' environment (includes Python, R, CDO, GNU Parallel)
3. Internet connection (for downloading climate data and SSURGO soils)

Activate the environment before running any commands:

    conda activate nps_wb


SETTING UP A NEW SITE
---------------------
Before running the model, you must set up the input directory structure
and gather the required data for your site.

Step 1: Choose a site name

    Use a short, machine-readable name with no spaces or special characters.
    Examples: "avalanche_peak", "shoshone_1"

Step 2: Create the directory structure

    mkdir -p data/input/{site_name}/shapefile
    mkdir -p data/input/{site_name}/dem
    mkdir -p data/input/{site_name}/soil

    Example for a site called "my_site":

    mkdir -p data/input/my_site/shapefile
    mkdir -p data/input/my_site/dem
    mkdir -p data/input/my_site/soil

Step 3: Add your site boundary shapefile

    Place your site boundary shapefile (all components) in the shapefile folder:
        data/input/{site_name}/shapefile/boundary.shp
        data/input/{site_name}/shapefile/boundary.shx
        data/input/{site_name}/shapefile/boundary.dbf
        data/input/{site_name}/shapefile/boundary.prj

    The shapefile should contain only ONE polygon object defining the area
    to be studied. It can be in any CRS - it will be reprojected automatically.

    Due to memory constraints, site size should be below ~100-150 acres.
    At 100 acres, sites need approximately 10 GB of RAM.

Step 4: Download the USGS 1m LiDAR DEM

    The model requires a 1-meter resolution DEM from the USGS 3DEP program.

    Download instructions:

    1. Go to the National Map download service:
       https://apps.nationalmap.gov/downloader/

    2. Click "Define Area of Interest" and either:
       - Draw a box around your site, OR
       - Upload your shapefile to set the bounding box automatically

    3. In the left panel, select:
       Data > Elevation Products (3DEP) > 1 meter DEM

    4. Set file format to: GeoTIFF

    5. Click "Search Products" to see available tiles

    6. Add the required tiles to your cart and download
       - If your site spans multiple tiles, you will need to merge them
         using GIS software (QGIS: Raster > Miscellaneous > Merge)

    7. Save the final GeoTIFF to:
       data/input/{site_name}/dem/USGS_1m.tif

    IMPORTANT:
    - All other inputs will be reprojected to the CRS of the DEM. Ensure that it is a suitable CRS for your final outputs.
    - The DEM must cover the ENTIRE site boundary
    - Not all areas have 1m coverage; check availability before selecting sites
    - Alternative DEM resolutions should work. However, this repository has only been tested in the 1m LiDAR DEM.

Step 5: Check for metdata gridcell overlap

    Before registering the site, check if your site polygon overlaps multiple
    metdata grid cells. This step must be performed manually using GIS software
    (QGIS recommended).

    1. Load your site polygon and data/metdata_elevationdata.nc into QGIS
    2. Visually examine whether the polygon falls within a single gridcell
       or crosses cell boundaries

    If NO overlap (site is entirely within one gridcell):
        - Use the site centroid coordinates for sites.csv
        - Use the elevation value from that gridcell for metdata_elev

    If site OVERLAPS multiple gridcells:
        - Choose a representative point that falls inside a SINGLE gridcell
        - Select the gridcell whose elevation most closely matches your
          site's mean DEM elevation (this improves lapse rate corrections)
        - Use that point's coordinates and the gridcell's elevation value

Step 6: Register the site in sites.csv

    The file src/sites.csv stores site metadata needed for climate data
    retrieval and lapse rate corrections. Format:

        site,lon,lat,metdata_elev

    - site: Your site name (must match directory name exactly)
    - lon: Longitude of representative point (decimal degrees, negative for west)
    - lat: Latitude of representative point (decimal degrees)
    - metdata_elev: Elevation from the 4km metdata grid cell (meters)

    Edit src/sites.csv and add a row with your site information:

        Example row:
        my_site,-110.1234,44.5678,2450

    To get the metdata_elev value, use the "Identify Features" tool in QGIS
    to click on the metdata_elevationdata.nc raster at your chosen point.

    Example sites.csv:

        site,lon,lat,metdata_elev
        avalanche_peak,-110.0832,44.4921,2847
        test,-109.8745,44.2156,2650
        my_site,-110.1234,44.5678,2450

    TODO: Interactive registration mode (01_clim_data.R run interactively)
    is incomplete and not recommended. Use manual entry instead.


ALSO REQUIRED (included in repository)
--------------------------------------
These files must exist in the data/ directory:

    data/merged_jennings2.tif        - Jennings T50 snow coefficients
    data/metdata_elevationdata.nc    - 4km elevation grid for lapse rates


WORKFLOW OVERVIEW
-----------------
Phase 1: Data Preparation    (R)      - Derive terrain layers, download soils
Phase 2: Climate Retrieval   (R)      - Download GridMET/MACA climate data
Phase 3: Water Balance Model (Python) - Calculate daily AET and Deficit
Phase 3b: Annual Aggregation (Bash)   - Sum daily grids to annual totals


PHASE 1: DATA PREPARATION
-------------------------
Prepares DEM derivatives (slope, aspect, hillshade) and downloads SSURGO
soil water holding capacity data.

    Rscript src/00_prep_data.R \
        --name={site_name} \
        --shapefile=data/input/{site_name}/shapefile/{name}.shp \
        --dem=data/input/{site_name}/dem/USGS_1m.tif

Outputs created in data/input/{site_name}/:
    dem/dem_nad83.tif           - Cropped 1m DEM
    dem/slope_nad83.tif         - Slope (radians)
    dem/aspect_nad83.tif        - Aspect (degrees)
    dem/hillshade_nad83.tif     - Hillshade for visualization
    jennings_t50_coefficients.tif
    soil/soil_whc_025.tif       - Soil water holding capacity (mm)


PHASE 2: CLIMATE DATA RETRIEVAL
-------------------------------
Downloads daily climate data from GridMET (historical) and MACA (projections).

Option A - Interactive mode (recommended for new sites):

    Rscript src/01_clim_data.R

    Follow prompts to enter site name. This will:
    - Extract site coordinates from the shapefile
    - Look up metdata grid elevation
    - Optionally append site info to src/sites.csv (say 'y' to register)

Option B - Batch mode (after site is registered in sites.csv):

    Rscript src/01_clim_data.R {site_name}

    Or for all registered sites:

    bash src/01_get_climate_data_batch.sh

Outputs created in data/input/{site_name}/:
    gridmet_1979_2023.csv       - Historical daily climate (1979-2024)
    macav2metdata_2006_2099.csv - Projected daily climate (2006-2099)

NOTE: Climate downloads can be slow. Consider running in screen/tmux.


PHASE 3: WATER BALANCE MODEL
----------------------------
Runs the gridded water balance model to calculate daily AET and Deficit.

Single run:

    python src/02_start_wb_v_1_5.py {model} {scenario} {site_name}

Arguments:
    model    - "historical" | "CanESM2" | "HadGEM2-CC365" | "MRI-CGCM3"
    scenario - "gridmet" (for historical) | "rcp45" | "rcp85"

Examples:

    # Historical run (1979-2024)
    python src/02_start_wb_v_1_5.py historical gridmet my_site

    # Future projection (2006-2099)
    python src/02_start_wb_v_1_5.py CanESM2 rcp85 my_site

Run all scenarios for all sites (requires GNU Parallel):

    bash src/02_run_all_historical_projections.sh

Outputs created in output/{site_name}/wb/:
    {model}_{scenario}_{year}_AET.nc     - Daily actual evapotranspiration
    {model}_{scenario}_{year}_Deficit.nc - Daily climatic water deficit

!! IMPORTANT: Values are stored as mm * 10 (16-bit integer compression).
   Divide by 10 to get actual mm values. See OUTPUT UNITS warning above.


PHASE 3b: ANNUAL AGGREGATION
----------------------------
Sums daily grids to annual totals using CDO.

    bash src/03_annual_sum.sh

Outputs created in output/{site_name}/sums/:
    {model}_{scenario}_AET_annual_sum.nc
    {model}_{scenario}_Deficit_annual_sum.nc


QUICK START EXAMPLE
-------------------
Run the complete workflow using the included "test" site.

The following site setup steps are already done for "test", but are shown
here as an example for new sites:

    # Activate environment
    conda activate nps_wb

    # Create directory structure (already exists for test)
    mkdir -p data/input/test/shapefile
    mkdir -p data/input/test/dem
    mkdir -p data/input/test/soil

    # Copy input files into place (already exists for test)
    cp /path/to/boundary.* data/input/test/shapefile/
    cp /path/to/USGS_1m.tif data/input/test/dem/

    # Add site to sites.csv (already done for test)
    # See Step 5-6 above for how to determine these values
    echo "test,-109.8745,44.2156,2650" >> src/sites.csv

Run the workflow (run from project root directory /home/ubuntu/nps-microclimate-water-balance):

    # Phase 1: Prepare data (DEM derivatives and soil data)
    Rscript src/00_prep_data.R \
        --name=test \
        --shapefile=data/input/test/shapefile/sample.shp \
        --dem=data/input/test/dem/USGS_1m.tif

    # Phase 2: Download climate data
    bash src/01_get_climate_data_batch.sh

    # Phase 3: Run water balance (historical only for quick test)
    python src/02_start_wb_v_1_5.py historical gridmet test

    # Phase 3b: Annual sums
    bash src/03_annual_sum.sh

NOTE: Phase 2 climate downloads can take ~1 hour due to server rate limits.


VISUALIZING RESULTS
-------------------
After the model run completes, use the example R Markdown notebook to visualize
inputs and outputs:

    reports/example-washburn.Rmd

This notebook produces diagnostic maps of model inputs (DEM, slope, aspect,
soil), boxplots of microclimatic variability across scenarios, spatial maps
of AET and CWD, and algorithmically selected planting locations.

To render the notebook from the command line:

    Rscript -e 'rmarkdown::render("reports/example-washburn.Rmd")'

Or open it in RStudio and click "Knit".

If running the model on a remote VM, retrieve the output files first:

    # From your local machine, copy output from the VM
    scp -r user@vm-host:/path/to/nps-microclimate-water-balance/output/test/sums \
        output/test/sums

The notebook is preconfigured for the "test" site. See the Configuration
section at the top of the notebook for paths and analysis parameters. It can
be adapted for additional sites by changing site_name, display_name, and the
corresponding paths.


TROUBLESHOOTING
---------------
- "Latitude outside of range": Site must be between 30-60°N latitude
- Memory errors: Reduce site size or increase VM RAM (~10GB per 100 hectares)
- CDO not found: Ensure nps_wb environment is activated
- SSURGO download fails: Check internet connection, try again
- Site not found: Check that site name matches directory name and sites.csv entry


FILE STRUCTURE
--------------
src/
    00_prep_data.R                  - Phase 1: Data preparation
    01_clim_data.R                  - Phase 2: Climate retrieval
    01_get_climate_data_batch.sh    - Phase 2: Batch wrapper
    02_start_wb_v_1_5.py            - Phase 3: Water balance model
    02_run_all_historical_projections.sh - Phase 3: Run all scenarios
    03_annual_sum.sh                - Phase 3b: Annual aggregation
    sites.csv                       - Site registry (name, lon, lat, elev)

data/
    merged_jennings2.tif            - Jennings T50 coefficients (required)
    metdata_elevationdata.nc        - 4km elevation grid (required)
    input/{site}/                   - Site-specific inputs
        shapefile/                  - Site boundary polygon
        dem/                        - USGS 1m LiDAR DEM
        soil/                       - SSURGO data (created by Phase 1)

output/{site}/
    wb/                             - Daily water balance grids
    sums/                           - Annual aggregated grids
