# Containerfile for NPS Gridded Water Balance Model
# Build: podman build -t nps-wb:latest .
# Run:   podman run -it --rm -v ./data:/app/data -v ./output:/app/output nps-wb:latest

FROM docker.io/rocker/geospatial:4.5.0

LABEL maintainer="NPS Water Balance Project"
LABEL description="High-resolution (1m) water balance model for the National Park Service"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Python and pip
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    # CDO for NetCDF operations
    cdo \
    # GNU Parallel for job parallelization
    parallel \
    # GDAL Python bindings dependencies
    libgdal-dev \
    gdal-bin \
    # NetCDF libraries
    libnetcdf-dev \
    libhdf5-dev \
    # Additional geospatial libs
    libproj-dev \
    libgeos-dev \
    libudunits2-dev \
    # Build tools
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set GDAL environment variables for Python bindings
ENV GDAL_CONFIG=/usr/bin/gdal-config
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

# Create app directory
WORKDIR /app

# Copy requirements first for better layer caching
COPY requirements.txt .

# Create Python virtual environment and install packages
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python dependencies
# Note: GDAL Python bindings must match system GDAL version, so we install it separately
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    grep -v "^GDAL" requirements.txt > requirements-nogdal.txt && \
    pip install --no-cache-dir -r requirements-nogdal.txt && \
    pip install --no-cache-dir "GDAL==$(gdal-config --version)"

# Copy renv files for R package restoration
COPY renv.lock .
COPY .Rprofile .
COPY renv/ renv/

# Restore R packages using renv
RUN R -e "renv::restore(prompt = FALSE)"

# Copy the rest of the application
COPY src/ src/

# Create data and output directories
RUN mkdir -p data/input data/metdata output

# Set the default working directory
WORKDIR /app

# Default command shows help
CMD ["bash", "-c", "echo 'NPS Water Balance Model Container' && echo '' && echo 'Usage examples:' && echo '  # Data preparation:' && echo '  Rscript src/00_prep_data.R --name=test --shapefile=data/input/test/shapefile/sample.shp --dem=data/input/test/dem/USGS_1m.tif' && echo '' && echo '  # Climate data retrieval:' && echo '  bash src/01_get_climate_data_batch.sh' && echo '' && echo '  # Water balance modeling:' && echo '  python src/02_start_wb_v_1_5.py historical gridmet test' && echo '' && echo '  # Annual aggregation:' && echo '  bash src/03_annual_sum.sh' && echo '' && echo 'Run with bash to get an interactive shell.'"]
