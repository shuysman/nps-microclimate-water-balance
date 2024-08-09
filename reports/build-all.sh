#!/bin/bash

set -e
set -u
set -x
set -o pipefail


parallel -j 3 Rscript -e "rmarkdown::render\(\'{}\'\)" ::: *.Rmd

scp -B burroughs.html web:/www/data/research/burroughs/index.html
scp -B avalanche.html web:/www/data/research/core_areas/avalanche/index.html
scp -B static_east.html web:/www/data/research/core_areas/static_east/index.html
scp -B static_west.html web:/www/data/research/core_areas/static_west/index.html
scp -B surprise.html web:/www/data/research/core_areas/surprise/index.html
scp -B holly_lake_small.html web:/www/data/research/core_areas/holly_lake_small/index.html

cd ../for_libby/
rm grte_yell_wb.zip
zip -9 grte_yell_wb.zip */*.tif
cd ../reports/
