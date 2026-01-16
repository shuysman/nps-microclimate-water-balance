# NPS Gridded Water Balance Model

This repository contains a high-resolution, gridded water balance model developed for ecological applications within the National Park Service. The project is a modified version of the [NPS Gridded Water Balance Model](http://www.yellowstoneecology.com/research/Gridded_Water_Balance_Model_Version_2_User_Manual.pdf) (Tercek et al. 2021b; Tercek et al. 2023) designed to run at finer resolutions and with some additional features such as empirically derived lapse rate corrections. The core water balance logic is the same between the two projects.

The model is designed to downscale coarse climate datasets (gridMET for historical analysis and MACA for future projections) to a 1-meter resolution. It integrates fine-scale topographic data derived from LiDAR and soil water holding capacity from SSURGO to simulate water balance components—such as Actual Evapotranspiration (AET) and Climatic Water Deficit (CWD)—at the microsite level.

This detailed, localized output supports critical conservation and management decisions, such as identifying optimal planting locations for climate-sensitive species like Whitebark Pine. The workflow is built for batch processing multiple sites on high-performance computing (HPC) clusters, but can also be run on a workstation or laptop with sufficient resources.

![2002-2022 Average Annual Climatic Water Deficit - Surprise and Amphitheater, GRTE](https://github.com/user-attachments/assets/c366ca37-40a1-4cf6-9676-012ead12c62b)

> **Important: Geographic Limitations**
>
> This model contains a routine to adjust daily temperature values based on lapse rate correction factors estimated for Mt. Washburn by Tercek et al. 2021a. These corrections are likely to be inaccurate for locations distant from Mt. Washburn (Greater Yellowstone Ecosystem) and should be altered or disabled for other regions. See the [Lapse Rate Corrections](#lapse-rate-corrections) section for details.

## Table of Contents

- [Requirements](#requirements)
  - [R Environment](#r-environment)
  - [Python Environment](#python-environment)
  - [Other Recommended Software](#other-recommended-software)
- [Container Installation](#container-installation)
- [Quick Start](#quick-start)
- [Model Run Instructions](#model-run-instructions)
  - [Phase 1: Site Setup](#phase-1-site-setup)
  - [Phase 2: Climate Data](#phase-2-climate-data)
  - [Phase 3: Run Water Balance Model](#phase-3-run-water-balance-model)
  - [Phase 4: Post-Processing](#phase-4-post-processing)
- [Output](#output)
- [Lapse Rate Corrections](#lapse-rate-corrections)
- [Reporting](#reporting)
- [Troubleshooting](#troubleshooting)
- [Acknowledgments](#acknowledgments)
- [License](#license)
- [References](#references)

## Requirements

This project relies on a combination of Python, R, and shell scripting and has only been tested on Linux environments (Debian 12/13 and Rocky Linux 8). The provided sbatch files can be used to submit intensive parts of the pipeline to an HPC cluster running Slurm Workload Manager.

Some manual data preparation is required. QGIS is recommended for these steps but any GIS providing basic raster and vector manipulations should also work.

### R Environment

This code was tested with R 4.5.0. Initial data prep, climate data retrieval steps, and included reports require R. An R environment with the required packages can be installed using `renv` with the included `renv.lock` file:

```bash
R -e "renv::install()"
```

### Python Environment

The water balance model (`02_start_wb_v_1_5.py`) requires Python 3. Use `requirements.txt` with `venv` (recommended) to set up a Python environment:

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Alternatively, use `spec-file.txt` with Conda/Mamba (deprecated).

### Other Recommended Software

- **GIS:** [QGIS](https://qgis.org/) for initial data preparation
- **Environment Management:**
  - [renv](https://rstudio.github.io/renv/) for managing R packages
  - Python virtual environments (Recommended)
  - [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) or  [Conda](https://docs.conda.io/en/latest/miniconda.html) (HPC using SLURM)
- **Command-Line Utilities:**
  - [Climate Data Operators (CDO)](https://code.mpimet.mpg.de/projects/cdo/): For processing NetCDF files. Tested with CDO versions 2.1.1 and 2.5.1.
  - [GNU Parallel](https://www.gnu.org/software/parallel/): For executing jobs in parallel.

## Container Installation

A container image is available as an alternative to manual environment setup. The container includes all Python and R dependencies pre-installed. Requires [Podman](https://podman.io/) or [Docker](https://www.docker.com/).

### Building the Container

```bash
# Standard build (4 cores for package compilation)
podman build -t nps-wb:latest .

# Custom core count for faster builds
podman build --build-arg NCPUS=16 -t nps-wb:latest .
```

Replace `podman` with `docker` if using Docker.

### Running the Container

Mount the `data/` and `output/` directories to persist input data and results:

```bash
# Interactive shell
podman run -it --rm \
  -v ./data:/app/data \
  -v ./output:/app/output \
  nps-wb:latest bash

# Run a specific command
podman run --rm \
  -v ./data:/app/data \
  -v ./output:/app/output \
  nps-wb:latest \
  python src/02_start_wb_v_1_5.py historical gridmet test
```

### Container Quick Start

```bash
# Build the container
podman build -t nps-wb:latest .

# Run the complete workflow for the test site
podman run --rm -v ./data:/app/data -v ./output:/app/output nps-wb:latest \
  bash -c "Rscript src/00_prep_data.R --name=test --shapefile=data/input/test/shapefile/sample.shp --dem=data/input/test/dem/USGS_1m.tif && \
           bash src/01_get_climate_data_batch.sh && \
           python src/02_start_wb_v_1_5.py historical gridmet test && \
           bash src/03_annual_sum.sh"
```

## Quick Start

An example site is included in `data/input/test/`. To run the complete workflow:

```bash
# Set up environments (only needed once)
R -e "renv::install()"
pip install -r requirements.txt

# Run the workflow
Rscript src/00_prep_data.R --name=test --shapefile=data/input/test/shapefile/sample.shp --dem=data/input/test/dem/USGS_1m.tif
bash src/01_get_climate_data_batch.sh
python src/02_start_wb_v_1_5.py historical gridmet test
bash src/03_annual_sum.sh
```

## Model Run Instructions

Data prep and model code is located in `src/`. The model run is separated into phases:

1. Site setup and data prep (R): `00_prep_data.R`
2. Climate data retrieval (R): `01_clim_data.R`
3. Water balance model (Python): `02_start_wb_v_1_5.py`
4. Post-processing (Bash): `03_annual_sum.sh`

Steps after running `00_prep_data.R` require site information to be configured in `src/sites.csv`.

### Phase 1: Site Setup

#### Create Site Directory

Create (or receive) a shapefile (ESRI Shapefile format) for a single area of interest. Depending on computer memory constraints, the site size should be below around 100-150 hectares. At around 100 hectares, sites need approximately 10 GB of RAM to run the water balance model.

Create a directory structure in `data/input/` with the desired site name:

```
data/input/{site_name}/
├── shapefile/
├── dem/
└── soil/
```

The site name used for the input directory must be kept consistent with other steps (e.g., `sites.csv` if using batching). Place the input shapefile in the `shapefile/` directory. The shapefile should contain only **one polygon object** defining the area to be studied.

#### Retrieve 1m LiDAR DEM

1. Go to the [National Map download service](https://apps.nationalmap.gov/downloader/)
2. Upload your shapefile to set the bounding box for data retrieval
3. Select `Data` > `Elevation Products (3D Elevation Program Products and Services)` > `Subcategories` > `1 meter DEM`
4. File format: `GeoTIFF, IMG`
5. Click `Search Products`, add layers to cart (stitch together if AOI overlaps multiple files)
6. Save the GeoTIFF file to `dem/USGS_1m.tif`

#### Run Data Preparation Script

Run `00_prep_data.R` to automate the remaining data preparation steps:

```bash
Rscript src/00_prep_data.R --name=test --shapefile=data/input/test/shapefile/sample.shp --dem=data/input/test/dem/USGS_1m.tif
```

This script:
1. Crops the 1 m DEM to the extent of the site polygon
2. Generates slope, aspect, and hillshade layers
3. Resamples Jennings T₅₀ coefficients to the site resolution and extent
4. Retrieves soil water holding capacity at 25 cm depth (`aws025wta`) from SSURGO

### Phase 2: Climate Data

#### Check for Gridcell Overlap

**Before downloading climate data**, check if the site polygon is completely within a metdata grid cell or if it overlaps cell boundaries. This step must be performed manually using a GIS (QGIS recommended).

1. Load the site polygon and `data/metdata_elevationdata.nc` into your GIS
2. Visually examine the polygon for overlap with the metdata raster cells

![Overlapping and non-overlapping site polygons](./docs/gridcell-overlap-example.png)

- **No overlap**: Proceed to downloading climate data
- **Overlap**: Select a representative point within the polygon that falls inside a single gridcell. Choose the gridcell whose elevation most closely matches the site's mean DEM elevation. Enter this point's coordinates in `sites.csv`.

#### Configure sites.csv

Set up `src/sites.csv` with site metadata:

| Column | Description |
|--------|-------------|
| `site` | Machine-readable site name (matches directory name) |
| `lon` | Longitude of representative point |
| `lat` | Latitude of representative point |
| `metdata_elev` | Elevation of the metdata cell (for lapse rate corrections) |

Example `sites.csv`:

```csv
site,lon,lat,metdata_elev
holly_lake_small,-110.8011392,43.7922368,2912.56
burroughs,-109.674010,43.705940,2813.28
avalanche,-110.134319,44.482919,2790
static_west,-110.805497,43.675967,3056.44
static_east,-110.805497,43.675967,3056.44
surprise,-110.777570,43.729726,2872.52
```

#### Download Climate Data

Run `01_clim_data.R` to download gridMET and MACA climate data. Data is retrieved from 1979-01-01 to 2024-12-31 and saved to `gridmet_1979_2024.csv` and `macav2metdata_2006_2099.csv`.

> **Note**: This step can take around 1 hour per site due to rate limiting from the [MACA THREDDS server](http://thredds.northwestknowledge.net:8080/thredds/catalog.html).

**Batch mode (recommended):**

```bash
# All sites in sites.csv
bash src/01_get_climate_data_batch.sh

# Single site
Rscript src/01_clim_data.R $site_name
```

**Interactive mode:**

```r
source("src/01_clim_data.R")
```

### Phase 3: Run Water Balance Model

After completing data preparation and climate data retrieval, your site directory should look like this:

```
data/input/test/
├── shapefile/
│   └── sample.shp
├── dem/
│   ├── aspect_nad83.tif
│   ├── dem_nad83.tif
│   ├── hillshade_nad83.tif
│   └── slope_nad83.tif
├── gridmet_1979_2024.csv
├── macav2metdata_2006_2099.csv
├── jennings_t50_coefficients.tif
└── soil/
    └── soil_whc_025.tif
```

#### Workstation or Laptop

Run the water balance model with:

```bash
python src/02_start_wb_v_1_5.py {model} {scenario} {site}
```

| Run Type | Model | Scenario |
|----------|-------|----------|
| Historical | `historical` | `gridmet` |
| Future RCP 4.5 | GCM name (e.g., `CanESM2`) | `rcp45` |
| Future RCP 8.5 | GCM name (e.g., `HadGEM2-CC365`) | `rcp85` |

Available GCMs: `CanESM2`, `HadGEM2-CC365`, `MRI-CGCM3`

Example:

```bash
# Historical run
python src/02_start_wb_v_1_5.py historical gridmet test

# Future projection
python src/02_start_wb_v_1_5.py CanESM2 rcp85 test
```

#### HPC Cluster (SLURM)

Submit batch jobs for all scenarios:

```bash
sbatch src/02_batch_wb_historical.sbatch
sbatch src/02_batch_wb_projections.sbatch
```

These scripts run the water balance model for historical and projected scenarios, save daily AET and CWD grids, and use `cdo` to generate annual sum files.

### Phase 4: Post-Processing

Generate gridded annual sums of AET and CWD:

```bash
bash src/03_annual_sum.sh
```

## Output

Model output is saved to:

- `output/{site}/wb/` — Daily water balance grids (NetCDF)
- `output/{site}/sums/` — Annual aggregated grids (NetCDF)

The annual sum grids are easier to work with for final reporting due to lower memory requirements. Use daily grids only if finer temporal resolution is needed (e.g., June-July-August sums).

Currently, only AET and CWD are saved as output. To save additional variables, modify `02_start_wb_v_1_5.py`.

The [User guide for the NPS Gridded Water Balance Model Dataset](http://www.yellowstoneecology.com/research/Gridded_Water_Balance_Model_Version_2_User_Manual.pdf) contains detailed information on output file format and metadata. This repository is equivalent to version 1.5 in that document.

## Lapse Rate Corrections

`02_start_wb_v_1_5.py` adjusts daily temperature values using lapse rate correction factors estimated for Mt. Washburn (Tercek et al. 2021a). This correction modifies temperature based on the elevation difference between each 1 m DEM pixel and the 4 km metdata cell elevation.

**Important considerations:**

- Temperature adjustments can be significant when the 1 m DEM pixel elevation differs substantially from the metdata cell elevation
- These corrections were empirically derived for Mt. Washburn and may be inaccurate for distant locations
- For sites outside the Greater Yellowstone Ecosystem, consider modifying or disabling the lapse rate corrections in the Python script

## Reporting

The `reports/` directory contains R Markdown examples for generating site reports from model output. While much code can be reused across sites, new sites generally require custom reports.

## Troubleshooting

### Memory Issues

- **Symptom**: Process killed or out-of-memory errors
- **Solution**: Reduce site size to under 100 hectares, or run on a machine with more RAM (10+ GB recommended for 100 ha sites)

### Climate Data Download Failures

- **Symptom**: Incomplete or missing climate CSV files
- **Solution**: The MACA THREDDS server has rate limits. Wait and retry, or run downloads overnight

### Missing Dependencies

- **Symptom**: Import errors in Python or package not found in R
- **Solution**: Ensure environments are activated and dependencies installed:
  ```bash
  # Python
  source venv/bin/activate
  pip install -r requirements.txt

  # R
  R -e "renv::install()"
  ```

### CDO Errors During Post-Processing

- **Symptom**: `cdo: command not found` or NetCDF errors
- **Solution**: Install CDO via your package manager (e.g., `apt install cdo` or `conda install -c conda-forge cdo`)

## Acknowledgments

This work was supported by funding provided by the National Park Service through an agreement with the [Northern Rockies Conservation Cooperative](https://nrccooperative.org/).

## License

This project is provided for research and conservation purposes. Contact the authors for licensing inquiries.

## References

Tercek, M.T., Rodman, A., Woolfolk, S., Wilson, Z., Thoma, D., and Gross, J., 2021a, Correctly applying lapse rates in ecological studies: comparing temperature observations and gridded data in Yellowstone: Ecosphere, v. 12, p. e03451, doi:10.1002/ecs2.3451.

Tercek, M.T., Thoma, D., Gross, J.E., Sherrill, K., Kagone, S., and Senay, G., 2021b, Historical changes in plant water use and need in the continental United States: PLOS ONE, v. 16, p. e0256586, doi:10.1371/journal.pone.0256586.

Tercek, M.T., Gross, J.E., and Thoma, D.P., 2023, Robust projections and consequences of an expanding bimodal growing season in the western United States: Ecosphere, v. 14, p. e4530, doi:10.1002/ecs2.4530.
