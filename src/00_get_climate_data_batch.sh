#!/bin/bash

set -euxo pipefail

sites="avalanche_peak cub_creek chittenden"

for site in $sites; do
    Rscript 00_clim_data.R $site
done
