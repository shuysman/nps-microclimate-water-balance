#!/bin/bash

set -euxo pipefail

sites=$(awk -F, 'NR > 1 { print $1 }' src/sites.csv)

for site in $sites; do
    Rscript src/01_clim_data.R $site
done
