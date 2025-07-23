# NPS gridded water balance model

This repository contains a high-resolution, gridded water balance model developed for ecological applications within the National Park Service. 

The model is designed to downscale coarse climate datasets (GridMET for historical analysis and MACA for future projections) to a 1-meter resolution. It integrates fine-scale topographic data derived from LiDAR and soil water holding capacity from SSURGO to simulate key water balance components—such as Actual Evapotranspiration (AET) and Climatic Water Deficit—at the microsite level.

This detailed, localized output supports critical conservation and management decisions, such as identifying optimal planting locations for climate-sensitive species like Whitebark Pine. The workflow is built for batch processing multiple sites on high-performance computing (HPC) clusters.

[2002-2022 Average Annual Climatic Water Deficit - Surprise and Amphitheater, GRTE](https://github.com/user-attachments/assets/c366ca37-40a1-4cf6-9676-012ead12c62b)

# Requirements

This project relies on a combination of Python, R, and shell scripting and has only been tested on Linux environments (Debian 12/13 and Rocky Linux 8). The provided sbatch files can be used to submit intensive parts of the pipeline to an HPC cluster running Slurm Workload Manager with dependencies managed by Conda/Mamba. The provided spec-file.txt can be used to create a Conda environment with the required dependencies for the python components of pipeline.

Some manual data preparation is required. QGIS is recommended for these steps but any GIS providing basic raster and vector manipulations should also work.

**Software:**

*   **GIS:** [QGIS](https://qgis.org/) for initial data preparation.
*   **Environment Management:** [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) or [Conda](https://docs.conda.io/en/latest/miniconda.html) for managing software packages and environments.
*   **Command-Line Utilities:**
    *   [Climate Data Operators (CDO)](https://code.mpimet.mpg.de/projects/cdo/): For processing NetCDF files.
    *   [GNU Parallel](https://www.gnu.org/software/parallel/): For executing jobs in parallel.


**R Environment:**

Initial data prep and climate data retrieval steps require R and the following packages:
*   `climateR`
*   `glue`
*   `sf`
*   `terra`
*   `tidyverse`
*   `FedData`
*   `optparse`
*   `janitor`

**Python Environment (`nps-wb`):**

The water balance model (`02_start_wb_v_1_5.py`) requires Python 3 and the following packages:
*   `gdal`
*   `netcdf4`
*   `numpy`
*   `pandas`
*   `rioxarray`
*   `utm`
*   `xarray`


# Model Run Instructions
An example is included in `/data/input/test` which can be used to run the model for testing. The following steps can be used to reproduce the test files in order to demonstrate how to set up new sites for the model. The test site demonstrates a site polygon that *does not* overlap a metdata gridcell, an issue described in more detail below. This is a modified version of Mike Tercek's [NPS Gridded Water Balance Model](http://www.yellowstoneecology.com/research/Gridded_Water_Balance_Model_Version_2_User_Manual.pdf) designed to run at finer resolutions and with some additional features such as empirically derived lapse rate corrections. The core water balance logic is the same between the two scripts.

## Site Setup
Create (or receive) a shapefile (ESRI Shapefile format) for a single area of interest to run the water balance model. Depending on computer memory constraints, the site size should be to below around 100-150 hectares. At around 100 hectares, sites need approximately 10 GB of RAM to run the water balance model. 

Create a directory in data/input with the desired site name, i.e., `/data/input/holly_lake_small` and create subdirectories `shapefile`, `dem`, and `soil` in this directory. The site name used for the input directory needs be kept consistent with others steps (e.g., `sites.csv` if using batching). Place the input shapefile (extracted if compressed) in the `shapefile` directory. 

Retrieve 1 m USGS LiDAR data for the area of interest
   - use National Map download service ( https://apps.nationalmap.gov/downloader/ )
   - Upload shapefile to set bounding box for data retrieval. You can use the original NPS/USFS planting area polygon shapefile here instead of bbox.gpkg, because the 1m LiDAR is provided in big chunks and we'll trim it down to the bbox later.
   - Select `Data` > `Elevation Products (3D Elevation Program Products and Services)` > `Subcategories` > `1 meter DEM`
   - File format:  `GeoTIFF, IMG`
   - Click `Search Products`, click the little shopping cart to add one layer to cart or add all to cart and stitch together if AOI overlaps multiple files.
   - Save the GeoTIFF file to `dem/USGS_1m.tif`
   
Once the input shapefile and 1 m DEM are saved in the correct locations, run `00_prep_data.R` which automates the rest of the data preparation steps. This script:
1. Crops the 1 m DEM to the extent of the site polygon input. 
2. Generates and saves slope, aspect, and hillshade layers from the cropped 1 m DEM.
3. Resamples and crops Jennings T<sub>50</sub> coefficients to the resolution and area covered by the cropped 1 m DEM.
 4. Retrieves soil water holding capacity at 25 cm soil depth (`aws025wta`) from the SSURGO database and creates rasterized version compatible with other GIS layers.
   
## Climate Data
### Check for gridcell overlap
**Before running any of the following steps** you need to check if the site polygon is completely within a metdata grid cell, or if it overlaps the edges/corners between cells. This step currently needs to be performed manually using a GIS of your choice (QGIS recommended). Load the site polygon and included `data/metdata_elevationdata.nc` file into your GIS. Visually example the polygon for potential overlap with the `metdata_elevationdata.nc` raster. The following image illustrates two polygons, one overlapping and one not overlapping with the metdata gridcells:

![Overlapping and non-overlapping site polygons](./docs/gridcell-overlap-example.png)

If no overlap, proceed to downloading the climate data. If there is overlap, select a point in the polygon that is within a gridcell that is representative of the site. I recommend selecting the metdata gridcell that most closely matches the elevation of the site (average elevation of all 1 m DEM pixels within the site polygon). If using the `00_clim_data.R` script in interactive mode, you will be prompted with choices to help facilitate this selection (not yet implemented). `00_get_climate_data_batch.sh` can be used to non-interactively batch download climate data for all sites in `sites.csv`.

### (Optional) Batching with sites.csv
Set up sites.csv for batching. Enter the (machine-readable) site name, longitude, latitude, and metdata elevation. If the site overlaps multiple metdata cells, choose a representative one (similar elevation to mean 1 m USGS DEM elevation for the site) and use that elevation here. Metdata elevation is used to calculate lapse rate adjustments in the water balance model.
example:
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


## Run water balance model
After running `00_prep_data.R` and `00_clim_data.R` (or `00_get_climate_data_batch.sh`), your `../data/input/{site}` directory should look like this:
```
test/
├── shapefile
│   └── sample.shp
├── dem
│   ├── aspect_nad83.tif
│   ├── dem_nad83.tif
│   ├── hillshade_nad83.tif
│   └── slope_nad83.tif
├── gridmet_1979_2023.csv
├── macav2metdata_2006_2099.csv
├── jennings_t50_coefficients.tif
└── soil
    └── soil_whc_025.tif
```

(If running on a SLURM cluster) run `02_batch_wb_historical.sbatch` and `02_batch_wb_projections.sbatch`

Model output will be saved `$HOME/out/$SITE/wb` (daily water balance grids) and `$HOME/out/$SITE/sums` (annual sum grids). The annual sum grids can be easier to work with for final reporting unless metrics are needed at a finer timestep, e.g., seasonal (summer) sum of AET. Currently, only AET and CWD are saved as model output. If other variables are desired, `02_start_wb_v_1_5.py` can be modified. 

The Rmd in `reports/` give examples of how to generate reports from the final output. Much code from these reports can be reused for new sites but in general new sites will require new Rmd reports to be written.
