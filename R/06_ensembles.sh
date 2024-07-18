#!/bin/bash
set -euxo pipefail

#data_dir=/media/smithers/shuysman/data/burroughs/sums
#out_dir=/media/smithers/shuysman/data/burroughs/sums/ensembles

data_dir=~/out/sums
out_dir=~/out/sums/ensembles

mkdir $out_dir

cdo -P 64 ensmean ${data_dir}/*_rcp45_Deficit_*.nc ${out_dir}/ensemble_rcp45_Deficit_annual_sum_burroughs.nc
cdo -P 64 ensmean ${data_dir}/*_rcp85_Deficit_*.nc ${out_dir}/ensemble_rcp85_Deficit_annual_sum_burroughs.nc
cdo -P 64 ensmean ${data_dir}/*_rcp45_AET_*.nc ${out_dir}/ensemble_rcp45_AET_annual_sum_burroughs.nc
cdo -P 64 ensmean ${data_dir}/*_rcp85_AET_*.nc ${out_dir}/ensemble_rcp85_AET_annual_sum_burroughs.nc
