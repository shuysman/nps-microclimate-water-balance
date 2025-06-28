# 1m Water Balance model for WBP core area analysis


https://github.com/user-attachments/assets/c366ca37-40a1-4cf6-9676-012ead12c62b

2002-2022 Average Annual Climatic Water Deficit - Surprise and Amphitheater, GRTE

# Site Setup
Create a directory in data/input with the desired site name, i.e., `/data/input/holly_lake_small`
Create subdirectories dem and soil

dem, aspect, slope, hillshade layers are created manually, as well as soils gpkg with soil polygons (in QGIS)
run `01_resample_layers.R` to create resampled versions of `1980_dayl` (not needed? may be able to skip this file),
jennings, soil_whc layers.

Add site id and lat long coordinates to `sites.csv`.  The whole site will be simulated with the climate at that point
Run `00_clim_data.R` to pull climate data for that site.  

```
metdata_elevationdata.nc
burroughs/
├── 1980_dayl_resampled.nc4
├── dem
│   ├── aspect_nad83.tif
│   ├── dem_nad83.tif
│   ├── hillshade_nad83.tif
│   └── slope_nad83.tif
├── gridmet_1979_2023.csv
├── bbox.gpkg
├── jennings_t50_coefficients.tif
├── macav2metdata_2006_2099.csv
└── soil
    ├── soil_whc_025.tif
    ├── soil_whc_100.tif
    └── ssurgo_soils.gpkg
```

Example site data available at https://huysman.net/research/core_areas/data.zip


# Steps
The initial data prep is currently performed in QGIS. ArcGIS or any other GIS with basic geospatial operations should also work fine. You need to manually prepare a bounding box, 1 m DEM layer, 1 m slope layer, 1 m aspect layer, 1 m hillshade layer, and 1 m soil WHC raster following these instructions.
1. Set up a directory in `input/` for each AOI. Pick a machine-readable name for each site, and use this naming consistently throughout. I.e., "Avalanche Peak" -> `input/avalanche_peak`
2. Make bbox.gpkg for bounding box from input shapefile for Area of Interest (AOI). This will be the area of 1m pixels to run the water balance model. Note that the outer border of 1 m pixels will be cut off the final water balance calculations, because determination of slope "eats" up those pixels. You can buffer the shapefile before this step if you are concerned about information loss for those pixels, which is likely to be unnoticeable at the final analysis extent.
   - At this step check if the bounding box overlaps with a gridMET grid cell. If so, you want to determine which of the overlapping grid cells most closely resembles the average elevation of the AOI.
   1. Calculate average elevation using 1m DEM
   2. Compare with elevation in `metdata_elevationdata.nc`
3. Retrieve 1m USGS LiDAR data for AOI bbox
   - use National Map download service ( https://apps.nationalmap.gov/downloader/ )
   - Upload shapefile to set bounding box for data retrieval. You can use the original NPS/USFS planting area polygon shapefile here instead of bbox.gpkg, because the 1m LiDAR is provided in big chunks and we'll trim it down to the bbox later.
   - Select `Data` > `Elevation Products (3D Elevation Program Products and Services)` > `Subcategories` > `1 meter DEM`
   - File format:  `GeoTIFF, IMG`
   - Click `Search Products`, click the little shopping cart to add one layer to cart or add all to cart and stitch together if AOI overlaps multiple 
   - Save the tiff to [unit]/dem/USGS_1m.tif, add to QGIS project
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

