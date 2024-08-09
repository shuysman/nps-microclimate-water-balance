#!/bin/bash

set -e
set -u
set -x
set -o pipefail

for file in *.Rmd; do
    Rscript -e "rmarkdown::render('$file')"
done
