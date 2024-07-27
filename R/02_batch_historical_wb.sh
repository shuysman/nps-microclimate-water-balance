#!/bin/bash

set -euxo pipefail

export sites="avalanche surprise static_west static_east"
export models="historical"
export scenarios="gridmet"

conda run -n nps-wb parallel -j 4 python 02_start_wb_v_1_5.py {} ::: $models ::: $scenarios ::: $sites
