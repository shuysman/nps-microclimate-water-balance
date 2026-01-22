#!/bin/bash
#
# 03_annual_sum.sh - Aggregate Daily Water Balance to Annual Sums
#
# Converts daily AET and Deficit NetCDF grids from the water balance model
# into annual sum grids using CDO (Climate Data Operators). This is the final
# processing step before analysis and visualization.
#
# PROCESS:
#   1. For each year's daily NetCDF: compute annual sum (cdo yearsum)
#   2. Merge all years into a single time-series NetCDF (cdo mergetime)
#   3. Clean up intermediate files
#
# SCENARIOS PROCESSED:
#   Historical:  historical/gridmet (1979-2024)
#   Projections: CanESM2, HadGEM2-CC365, MRI-CGCM3 × rcp45, rcp85 (2006-2099)
#
# USAGE:
#   bash src/03_annual_sum.sh
#
# INPUTS:
#   src/sites.csv - List of sites to process
#   output/{site}/wb/{model}_{scenario}_{year}_{var}.nc - Daily water balance grids
#
# OUTPUTS:
#   output/{site}/sums/{model}_{scenario}_{var}_annual_sum.nc
#     - Variables: AET, Deficit (annual totals in mm*10)
#     - One file per model/scenario combination
#
# DEPENDENCIES:
#   - CDO (Climate Data Operators) must be installed and in PATH
#   - GNU Parallel for concurrent processing
#
# NOTE: Values are stored as mm*10 (multiplied by 10 for integer storage).
#       Divide by 10 when analyzing to get mm/year.
#
set -euxo pipefail

sites=$(awk -F, 'NR > 1 { print $1 }' src/sites.csv)

vars="AET Deficit"

#models="bcc-csm1-1-m bcc-csm1-1 BNU-ESM CanESM2 CNRM-CM5 CSIRO-Mk3-6-0 GFDL-ESM2G GFDL-ESM2M HadGEM2-CC365 HadGEM2-ES365 inmcm4 IPSL-CM5A-LR IPSL-CM5A-MR IPSL-CM5B-LR MIROC5 MIROC-ESM-CHEM MIROC-ESM MRI-CGCM3 NorESM1-M"
proj_models="CanESM2 HadGEM2-CC365 MRI-CGCM3"
proj_scenarios="rcp45 rcp85"

hist_models="historical"
hist_scenarios="gridmet"

calc_annual_sum () {
    model=$1
    scenario=$2
    var=$3
    site=$4
    echo $site $model $scenario $var

    in_dir="output/${site}/wb"
    out_dir="output/${site}/sums"
    
    if [ ! -d $out_dir ]; then
	mkdir $out_dir    
    fi

    for file in ${in_dir}/${model}_${scenario}_*_${var}.nc; do
	cdo yearsum $file "${file}_sum.nc"
    done

    cdo mergetime ${in_dir}/${model}_${scenario}_*_${var}.nc_sum.nc "${out_dir}/${model}_${scenario}_${var}_annual_sum.nc"
    rm ${in_dir}/${model}_${scenario}_*_${var}.nc_sum.nc
}

export -f calc_annual_sum

parallel -j 4 calc_annual_sum {1} {2} {3} {4} ::: $hist_models ::: $hist_scenarios ::: $vars ::: $sites
parallel -j 4 calc_annual_sum {1} {2} {3} {4} ::: $proj_models ::: $proj_scenarios ::: $vars ::: $sites
