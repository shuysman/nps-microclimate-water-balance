#!/bin/bash

set -e
set -u
set -x
set -o pipefail

for file in *.Rmd; do
    Rscript -e "rmarkdown::render('$file')"
done

cd ../for_libby/
rm grte_yell_wb.zip
zip -9 grte_yell_wb.zip */*.tif
cd ../reports/
