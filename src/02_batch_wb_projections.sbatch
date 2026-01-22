#!/bin/bash
##
## 02_batch_wb_projections.sbatch - SLURM Batch Script for Projected Water Balance
##
## Runs the water balance model for all sites using MACA downscaled climate
## projections (2006-2099). Processes 3 GCMs × 2 RCPs = 6 scenarios per site.
## Automatically runs annual aggregation (03_annual_sum.sh) after completion.
##
## GCM/SCENARIO MATRIX:
##   CanESM2       - rcp45, rcp85 (warm/hot-dry futures)
##   HadGEM2-CC365 - rcp45, rcp85 (hot futures)
##   MRI-CGCM3     - rcp45, rcp85 (warm-wet futures)
##
## USAGE:
##   sbatch src/02_batch_wb_projections.sbatch
##
## INPUTS:
##   src/sites.csv - List of sites to process
##   For each site: data/input/{site}/macav2metdata_2006_2099.csv
##
## OUTPUTS:
##   output/{site}/wb/{model}_{scenario}_{year}_{var}.nc - Daily NetCDFs
##   output/{site}/sums/ - Annual aggregated grids (via 03_annual_sum.sh)
##   wb-projections-{jobid}.out/.err - SLURM logs
##
## RESOURCE REQUIREMENTS:
##   64 CPUs, 256GB RAM, 1.5 day runtime
##   Runs 36 jobs in parallel (-j36): 6 scenarios × 6 concurrent sites
##
## PREREQUISITES:
##   - Mamba environments: "nps-wb" (Python), "cdo" (annual aggregation)
##   - Completed Phase 1 and Phase 2 for all sites
##
#SBATCH --account=priority-briansmithers        #specify the account to use
#SBATCH --job-name=wb-projections            # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
##SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64              # number of cores to allocate
#SBATCH --mem=256G                     # 2000 MB of Memory allocated; set --mem with care
##SBATCH --mem-per-cpu=100G
#SBATCH --time=1-12:00:00                 # Maximum job run time
##SBATCH --array=0-1847%10                  # Number of jobs in array
#SBATCH --mail-user=shuysman@gmail.com
#SBATCH --mail-type=ALL
##SBATCH --output=fire-%A-%a.out
##SBATCH --error=fire-%A-%a.err
#SBATCH --output=wb-projections-%j.out
#SBATCH --error=wb-projections-%j.err

date
hostname -s
module load Mamba
source $HOME/.bashrc
mamba activate nps-wb
sites=$(awk -F, 'NR > 1 { print $1 }' sites.csv)
parallel -j36 python 02_start_wb_v_1_5.py {1} {2} {3} :::  CanESM2 HadGEM2-CC365 MRI-CGCM3 ::: rcp45 rcp85 ::: $sites
mamba deactivate
mamba activate cdo
./03_annual_sum.sh
date
