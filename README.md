# NPS gridded water balance model

This repository contains a high-resolution, gridded water balance model developed for ecological applications within the National Park Service. The project is a modified version of the [NPS Gridded Water Balance Model](http://www.yellowstoneecology.com/research/Gridded_Water_Balance_Model_Version_2_User_Manual.pdf) (Tercek et al. 2021b; Tercek et al. 2023) designed to run at finer resolutions and with some additional features such as empirically derived lapse rate corrections. The core water balance logic is the same between the two projects.

The model is designed to downscale coarse climate datasets (gridMET for historical analysis and MACA for future projections) to a 1-meter resolution. It integrates fine-scale topographic data derived from LiDAR and soil water holding capacity from SSURGO to simulate water balance components—such as Actual Evapotranspiration (AET) and Climatic Water Deficit—at the microsite level.

This detailed, localized output supports critical conservation and management decisions, such as identifying optimal planting locations for climate-sensitive species like Whitebark Pine. The workflow is built for batch processing multiple sites on high-performance computing (HPC) clusters, but can also be run on a workstation or laptop with sufficient resources.

[2002-2022 Average Annual Climatic Water Deficit - Surprise and Amphitheater, GRTE](https://github.com/user-attachments/assets/c366ca37-40a1-4cf6-9676-012ead12c62b)

# Requirements

This project relies on a combination of Python, R, and shell scripting and has only been tested on Linux environments (Debian 12/13 and Rocky Linux 8). The provided sbatch files can be used to submit intensive parts of the pipeline to an HPC cluster running Slurm Workload Manager with dependencies managed by Conda/Mamba (deprecated).

Some manual data preparation is required. QGIS is recommended for these steps but any GIS providing basic raster and vector manipulations should also work.

**R Environment:**

This code was tested with `R` 4.5.0. Initial data prep, climate data retrieval steps, and included reports require R. An R environment with the required packages can be installed using `renv` using the included `renv.lock` file. If `renv` is installed, run `renv::install()` from the project root directory to install required R packages.

**Python Environment (`nps-wb`):**

The water balance model (`02_start_wb_v_1_5.py`) requires Python 3. Use `spec-file.txt` with Conda/Mamba (deprecated) or `requirements.txt` with `venv` (recommended) to set up a python environment with the required packages.

**Other Recommended Software:**

*   **GIS:** [QGIS](https://qgis.org/) for initial data preparation.
*   **Environment Management:** 
	* [renv](https://rstudio.github.io/renv/) for managing R packages.
	* [Mamba](https://mamba.readthedocs.io/en/latest/installation.html), [Conda](https://docs.conda.io/en/latest/miniconda.html), or python virtual environments for managing software packages and environments for python scripts.
*   **Command-Line Utilities:**
    *   [Climate Data Operators (CDO)](https://code.mpimet.mpg.de/projects/cdo/): For processing NetCDF files.
    *   [GNU Parallel](https://www.gnu.org/software/parallel/): For executing jobs in parallel.

# Model Run Instructions
An example is included in `data/input/test/` which can be used to run the model for testing. The following steps can be used to reproduce the test files in order to demonstrate how to set up new sites for the model. The test site demonstrates a site polygon that *does not* overlap a metdata gridcell, an issue described in more detail below. 

## Site Setup
Create (or receive) a shapefile (ESRI Shapefile format) for a single area of interest to run the water balance model. Depending on computer memory constraints, the site size should be to below around 100-150 hectares. At around 100 hectares, sites need approximately 10 GB of RAM to run the water balance model. 

Create a directory in data/input with the desired site name, i.e., `data/input/holly_lake_small/` and create subdirectories `shapefile`, `dem`, and `soil` in this directory. The site name used for the input directory needs be kept consistent with others steps (e.g., `sites.csv` if using batching). Place the input shapefile (extracted if compressed) in the `shapefile` directory. The shapefile should only have **1 polygon object** defining the area to be studied or undefined behavior may occur.

Retrieve 1 m USGS LiDAR data for the area of interest
   - use National Map download service ( https://apps.nationalmap.gov/downloader/ )
   - Upload shapefile to set bounding box for data retrieval. You can use the original NPS/USFS planting area polygon shapefile here instead of bbox.gpkg, because the 1m LiDAR is provided in big chunks and we'll trim it down to the bbox later.
   - Select `Data` > `Elevation Products (3D Elevation Program Products and Services)` > `Subcategories` > `1 meter DEM`
   - File format:  `GeoTIFF, IMG`
   - Click `Search Products`, click the little shopping cart to add one layer to cart or add all to cart and stitch together if AOI overlaps multiple files.
   - Save the GeoTIFF file to `dem/USGS_1m.tif`
   
Once the input shapefile and 1 m DEM are saved in the correct locations, run `00_prep_data.R` which automates the rest of the data preparation steps. The script requires arguments for site name (matching site data directory name in data/input/), shapefile, and source USGS 1 m DEM GeoTIFF file in the following format:

`Rscript src/00_prep_data.R --name=test --shapefile=data/input/test/shapefile/sample.shp --dem=data/input/test/dem/`

This script:
1. Crops the 1 m DEM to the extent of the site polygon input. 
2. Generates and saves slope, aspect, and hillshade layers from the cropped 1 m DEM.
3. Resamples and crops Jennings T<sub>50</sub> coefficients to the resolution and area covered by the cropped 1 m DEM.
 4. Retrieves soil water holding capacity at 25 cm soil depth (`aws025wta`) from the SSURGO database and creates rasterized version compatible with other GIS layers.
   
## Climate Data
### Check for gridcell overlap
**Before running any of the following steps** you need to check if the site polygon is completely within a metdata grid cell, or if it overlaps the edges/corners between cells. This step currently needs to be performed manually using a GIS of your choice (QGIS recommended). Load the site polygon and included `data/metdata_elevationdata.nc` file into your GIS. Visually example the polygon for potential overlap with the `metdata_elevationdata.nc` raster. The following image illustrates two polygons, one overlapping and one not overlapping with the metdata gridcells:

![Overlapping and non-overlapping site polygons](./docs/gridcell-overlap-example.png)

If no overlap, proceed to downloading the climate data. If there is overlap, select a point in the polygon that is within a gridcell that is representative of the site and enter this point for the coordinates in `site.csv`. I recommend selecting a point within the metdata gridcell that most closely matches the elevation of the site (average elevation of all 1 m DEM pixels within the site polygon). If using the `01_clim_data.R` script in interactive mode, you will be prompted with choices to help facilitate this selection (not yet implemented). `01_get_climate_data_batch.sh` can be used to non-interactively batch download climate data for all sites in `sites.csv`. 

### (Optional) Batching with sites.csv
Set up sites.csv for batching. Enter the (machine-readable) site name, longitude, latitude, and metdata elevation. If the site overlaps multiple metdata cells, choose a cell that is representative (the metdata cell has similar elevation to mean 1 m USGS DEM elevation for the site) and use coordinates for a point within that cell and the site polygon for `lon` and `lat` and enter the metdata cell elevation for that cell in `metdata_elev`. Metdata elevation is used to calculate lapse rate adjustments in the water balance model. Latitude/longitude coordinates are used to retrieve climate data in `01_clim_data.R` when run in batch mode and in some of the included site report files.

Example `sites.csv`:
```
site,               lon,            lat,        metdata_elev
holly_lake_small,   -110.8011392,   43.7922368, 2912.56
burroughs,          -109.674010,    43.705940,  2813.28
avalanche,          -110.134319,    44.482919,  2790
static_west,        -110.805497,    43.675967,  3056.44
static_east,        -110.805497,    43.675967,  3056.44
surprise,           -110.777570,    43.729726,  2872.52
```


### Download climate data

Run `01_clim_data.R` to download gridMET and MACA climate data for one site or `01_get_climate_data_batch.sh` to batch download for all sites in `sites.csv`. Climate data will be retrieved from 1979-01-01 to 2024-12-31 (all complete years on record) and saved to `gridmet_1979_2023.csv` and `macav2metdata_2006_2099.csv`. This step can take a long time (around 1 hour per site) because of rate limiting from the [MACA THREDDS server](http://thredds.northwestknowledge.net:8080/thredds/catalog.html).

You can either:
* Run the script in **batch mode** (recommended) — set up `sites.csv` and call `bash src/01_get_climate_data_batch.sh` to download all sites in sites.csv or `Rscript src/01_clim_data.R $site_name` to download climate data for one site.
* Run the script in **interactive mode** — run `source("src/01_clim_data.R")` in an R session to use interactive mode. This runs some prompts to help things like setting up finding required files and `sites.csv` . 

## Run water balance model
After running `00_prep_data.R` and `01_clim_data.R` (or `01_get_climate_data_batch.sh`), your `data/input/{site}/` directory should look like this:
```
test/
├── shapefile/
│   └── sample.shp
├── dem/
│   ├── aspect_nad83.tif
│   ├── dem_nad83.tif
│   ├── hillshade_nad83.tif
│   └── slope_nad83.tif
├── gridmet_1979_2023.csv
├── macav2metdata_2006_2099.csv
├── jennings_t50_coefficients.tif
└── soil/
    └── soil_whc_025.tif
```

(If running on a SLURM cluster) run `02_batch_wb_historical.sbatch` and `02_batch_wb_projections.sbatch`. This runs the water balance model for historical and projected scenarios, saves daily AET and CWD grids, and uses `cdo` to generate the annual sum files for AET and CWD.

(If running on a workstation or laptop) Run `python 02_start_wb_v_1_5.py $model $scenario $site` for each model and scenario combination you wish to run for each site. Historical gridMET runs use `historical` for model and `gridmet` for scenario. For projections, enter the GCM for model and either rcp45 or rcp85 for scenario.


**NOTE** `02_start_wb_v_1_5.py` contains a routine to adjust daily temperature values based on a lapse rate correction factor estimated for Mt. Washburn by Tercek et al. 2021a. This correction increases or decreases daily temperature based on the difference between the 1 m USGS LiDAR DEM pixel elevation and the 4 km metdata cell elevation. The different between original and adjusted temperature can be significant if the 1 m DEM pixel is significantly higher or lower in elevation than the metdata cell. The lapse rate corrections are likely to be inaccurate for locations distant from Mt. Washburn and should be altered or disabled for locations not in its vicinity.

Model output will be saved `output/$SITE/wb` (daily water balance grids) and `$output/$SITE/sums` (annual sum grids). The annual sum grids can be easier to work with for final reporting unless metrics are needed at a finer timestep, e.g., seasonal (summer) sum of AET. Currently, only AET and CWD are saved as model output. If other variables are desired, `02_start_wb_v_1_5.py` can be modified. 

The Rmd in `reports/` give examples of how to generate reports from the final output. Much code from these reports can be reused for new sites but in general new sites will require new Rmd reports to be written.

# References
Tercek, M.T., Rodman, A., Woolfolk, S., Wilson, Z., Thoma, D., and Gross, J., 2021a, Correctly applying lapse rates in ecological studies: comparing temperature observations and gridded data in Yellowstone: Ecosphere, v. 12, p. e03451, doi:10.1002/ecs2.3451.

Tercek, M.T., Thoma, D., Gross, J.E., Sherrill, K., Kagone, S., and Senay, G., 2021b, Historical changes in plant water use and need in the continental United States: PLOS ONE, v. 16, p. e0256586, doi:10.1371/journal.pone.0256586.

Tercek, M.T., Gross, J.E., and Thoma, D.P., 2023, Robust projections and consequences of an expanding bimodal growing season in the western United States: Ecosphere, v. 14, p. e4530, doi:10.1002/ecs2.4530.
