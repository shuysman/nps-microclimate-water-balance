#!/usr/bin/env bash
set -euo pipefail

sites=$(awk -F, 'NR > 1 { print $1 }' sites.csv)

## GCM/RCP combinations chosen to represent extremes in area of interest, following
## Lawrence, David J., Amber N. Runyon, John E. Gross, Gregor W. Schuurman, and Brian W. Miller. 2021.
## “Divergent, Plausible, and Relevant Climate Futures for near- and Long-Term Resource Planning.”
## Climatic Change 167 (3): 38. https://doi.org/10.1007/s10584-021-03169-y.

parallel -j32 python src/02_start_wb_v_1_5.py {1} {2} {3} \
  ::: historical CanESM2 HadGEM2-CC365 MRI-CGCM3 MRI-CGCM3 \
  :::+ gridmet rcp85 rcp85 rcp45 rcp85 \
  ::: $sites
