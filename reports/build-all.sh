#!/bin/bash

for file in *.Rmd; do
    Rscript -e "rmarkdown::render('$file')"
done
