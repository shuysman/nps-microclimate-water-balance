#!/bin/bash

set -euxo pipefail

sites=$(awk -F, 'NR > 1 { print $1 }' sites.csv)

for site in $sites; do
    Rscript 00_clim_data.R $site
done
