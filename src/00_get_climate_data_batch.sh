#!/bin/bash

set -euxo pipefail

sites="burroughs avalanche static_west static_east surprise holly_lake_small"

for site in $sites; do
    Rscript 00_clim_data.R $site
done
