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
## Site Setup
Create a directory in data/input with the desired site name, i.e., `/data/input/holly_lake_small` and create subdirectories `shapefile`, `dem`, and `soil` in this directory. The site name used for the input directory needs be kept consistent with others steps (e.g., `sites.csv` if using batching). Place the input shapefile (extracted if compressed) in the `shapefile` directory. 

Retrieve 1 m USGS LiDAR data for the area of interest
   - use National Map download service ( https://apps.nationalmap.gov/downloader/ )
   - Upload shapefile to set bounding box for data retrieval. You can use the original NPS/USFS planting area polygon shapefile here instead of bbox.gpkg, because the 1m LiDAR is provided in big chunks and we'll trim it down to the bbox later.
   - Select `Data` > `Elevation Products (3D Elevation Program Products and Services)` > `Subcategories` > `1 meter DEM`
   - File format:  `GeoTIFF, IMG`
   - Click `Search Products`, click the little shopping cart to add one layer to cart or add all to cart and stitch together if AOI overlaps multiple files.
   - Save the tiff to `dem/USGS_1m.tif`
   
Once the input shapefile and 1 m DEM are saved in the correct locations, run `00_prep_data.R` which automates the rest of the data preparation steps. This script:
1. Crops the 1 m DEM to the extent of the site polygon input. 
2. Generates and saves slope, aspect, and hillshade layers from the cropped 1 m DEM.
3. Resamples and crops Jennings T~50~ coefficients to the resolution and area covered by the cropped 1 m DEM.
4. Retrieves soil water holding capacity at 25 cm soil depth (`aws025wta`) from the SSURGO database and creates rasterized version compatible with other GIS layers.
   
## Climate Data
### Check for gridcell overlap
**Before running any of the following steps** you need to check if the site polygon is completely within a metdata grid cell, or if it overlaps the edges/corners between cells. This step currently needs to be performed manually using a GIS of your choice (QGIS recommended). Load the site polygon and included `data/metdata_elevationdata.nc` file into your GIS. Visually example the polygon for potential overlap with the `metdata_elevationdata.nc` raster. The following image illustrates two polygons, one overlapping and one not overlapping with the metdata gridcells:

![Overlapping and non-overlapping site polygons](./docs/gridcell-overlap-example.png)

If no overlap, proceed to downloading the climate data. If there is overlap, select a point in the polygon that is within a gridcell that is representative of the site. I recommend selecting the metdata gridcell that most closely matches the elevation of the site (average elevation of all 1 m DEM pixels within the site polygon). If using the `00_clim_data.R` script in interactive mode, you will be prompted with choices to help facilitate this selection (not yet implemented). `00_get_climate_data_batch.sh` can be used to non-interactively batch download climate data for all sites in `sites.csv`.

### Download climate data



## Steps
The initial data prep is currently performed in QGIS. ArcGIS or any other GIS with basic geospatial operations should also work fine. You need to manually prepare a bounding box, 1 m DEM layer, 1 m slope layer, 1 m aspect layer, 1 m hillshade layer, and 1 m soil WHC raster following these instructions.
1. Set up a directory in `input/` for each AOI. Pick a machine-readable name for each site, and use this naming consistently throughout. I.e., "Avalanche Peak" -> `input/avalanche_peak`
2. Make bbox.gpkg for bounding box from input shapefile for Area of Interest (AOI). This will be the area of 1m pixels to run the water balance model. Note that the outer border of 1 m pixels will be cut off the final water balance calculations, because determination of slope "eats" up those pixels. You can buffer the shapefile before this step if you are concerned about information loss for those pixels, which is likely to be unnoticeable at the final analysis extent.
   - At this step check if the bounding box overlaps with a gridMET grid cell. If so, you want to determine which of the overlapping grid cells most closely resembles the average elevation of the AOI.
   1. Calculate average elevation using 1m DEM
   2. Compare with elevation in `metdata_elevationdata.nc`
3. 
4. Crop raster using AOI bboxes
   - "Clip raster by extent" tool in QGIS
	 - Input layer: USGS_1m.tif
	 - Clipping extent: bbox.gpkg
	 - save as `input/[site]/dem/dem_nad83.tif`
5. Make hillshade from dem_nad83.tif
   - save as `input/[site]/dem/hillshade_nad83.tif`
   - elevation layer: `input/[site]/dem/dem_nad83.tif`
   - zfactor: 1
   - azimuth: 300
   - vertical angle: 40
6. Make slope from dem_nad83.tif
   - save as `input/[site]/dem/slope_nad83.tif`
   - elevation layer: `input/[site]/dem/dem_nad83.tif`
   - zfactor: 1
7. Make aspect from dem_nad83.tif
   - save as `input/[site]/dem/aspect_nad83.tif`
   - elevation layer: `input/[site]/dem/dem_nad83.tif`
   - zfactor: 1
8. Make soil layer
   - I use the "Curve Number Generator" plugin in QGIS
   - Processing Toolbox> Curve Number Generator > Curve Number Generator (CONUS)
   - Set AOI to bbox
   - Soils [optional] > save to `input/[site]/soil/ssurgo.gpkg`
9. Copy files to `data/`
   - 1980_dayl_na.nc4: This might not be needed anymore... Need to review code.
   - merged_jennings2.tif

10. Set up sites.csv for batching
example:
site,lon,lat, metdata_elev
holly_lake_small,  -110.8011392,43.7922368,2912.56
burroughs,-109.674010,43.705940,2813.28
avalanche,-110.134319, 44.482919,2790
static_west, -110.805497, 43.675967,3056.44
static_east, -110.805497, 43.675967,3056.44
surprise,-110.777570,  43.729726,2872.52

  - If polygons extend past gridMET cell boundaries, set a point here that is within a grid cell that is representative of the site (i.e., similar elevation). The actual point location is not super important as long as it is located within the grid cell you are targetting.
  - metdata_elev comes from `metdata_elevationdata.nc`
  
11. Run `src/00_clim_data_batch.sh`
  - retrieves gridMET and MACA data for each point in sites.csv

12. Run `src/01_resample_layers.R`
 - Resamples 1980_dayl_na.nc4, merged_jennings2.tif, and soils data using 1m DEM as reference

13. (If running on a SLURM cluster) run `02_batch_wb_historical.sbatch` and `02_batch_wb_projections.sbatch`

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
