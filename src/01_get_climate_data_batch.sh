#!/bin/bash
#
# 01_get_climate_data_batch.sh - Batch Climate Data Retrieval
#
# Wrapper script that runs 01_clim_data.R for all sites listed in src/sites.csv.
# Use this after adding site entries to sites.csv (via interactive mode of
# 01_clim_data.R or manually).
#
# INPUTS:
#   src/sites.csv - CSV with columns: site, lon, lat, metdata_elev
#
# OUTPUTS:
#   For each site, creates climate CSVs in data/input/{site}/:
#     - gridmet_1979_2023.csv
#     - macav2metdata_2006_2099.csv
#
# USAGE:
#   bash src/01_get_climate_data_batch.sh
#
# NOTE: Downloads can be slow. Consider running in screen/tmux for many sites.
#

set -euxo pipefail

sites=$(awk -F, 'NR > 1 { print $1 }' src/sites.csv)

for site in $sites; do
    Rscript src/01_clim_data.R $site
done
