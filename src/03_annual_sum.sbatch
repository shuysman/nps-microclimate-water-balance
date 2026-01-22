#!/bin/bash
##
## 03_annual_sum.sbatch - SLURM Batch Script for Annual Aggregation
##
## Wrapper script to run 03_annual_sum.sh on HPC cluster. Aggregates daily
## water balance outputs into annual sum grids using CDO.
##
## USAGE:
##   sbatch src/03_annual_sum.sbatch
##
## INPUTS:
##   output/{site}/wb/*.nc - Daily water balance NetCDFs from Phase 2
##
## OUTPUTS:
##   output/{site}/sums/{model}_{scenario}_{var}_annual_sum.nc
##   sum-{jobid}.out/.err - SLURM logs
##
## RESOURCE REQUIREMENTS:
##   128 CPUs, 512GB RAM, 1 day runtime
##   High memory needed for merging large time-series
##
## PREREQUISITES:
##   - Mamba environment "cdo" with CDO installed
##   - Completed water balance runs (02_start_wb_v_1_5.py)
##
## NOTE: This is automatically called at the end of 02_batch_wb_projections.sbatch,
##       but can be run separately if needed.
##
#SBATCH --account=priority-briansmithers        #specify the account to use
#SBATCH --job-name=annual-sum            # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
##SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128             # number of cores to allocate
#SBATCH --mem=512G                     # 2000 MB of Memory allocated; set --mem with care
##SBATCH --mem-per-cpu=100G
#SBATCH --time=1-00:00:00                 # Maximum job run time
##SBATCH --array=0-1847%10                  # Number of jobs in array
#SBATCH --mail-user=shuysman@gmail.com
#SBATCH --mail-type=ALL
##SBATCH --output=fire-%A-%a.out
##SBATCH --error=fire-%A-%a.err
#SBATCH --output=sum-%j.out
#SBATCH --error=sum-%j.err

date
hostname -s
module load Mamba
source $HOME/.bashrc
mamba activate cdo
./03_annual_sum.sh
date
