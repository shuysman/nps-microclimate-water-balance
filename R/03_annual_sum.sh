#!/bin/bash
set -euxo pipefail

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

    in_dir="${HOME}/out/${site}/wb"
    out_dir="${HOME}/out/${site}/sums"
    
    if [ ! -d $out_dir ]; then
	mkdir $out_dir    
    fi

    for file in ${in_dir}/${model}_${scenario}_*_${var}.nc; do
	cdo yearsum $file "${file}_sum.nc"
    done

    cdo mergetime ${in_dir}/${model}_${scenario}_*_${var}_sum.nc "${out_dir}/${model}_${scenario}_${var}_annual_sum.nc"
    rm ${in_dir}/${model}_${scenario}_*_${var}_sum.nc
}

export -f calc_annual_sum

parallel -j 16 calc_annual_sum {} ::: $hist_models ::: $hist_scenarios ::: $vars ::: avalanche burroughs static_west static_east surprise holly_lake_small
parallel -j 64 calc_annual_sum {} ::: $proj_models ::: $proj_scenarios ::: $vars ::: avalanche burroughs static_west static_east surprise holly_lake_small
