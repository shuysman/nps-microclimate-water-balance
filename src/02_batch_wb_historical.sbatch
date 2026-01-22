#!/bin/bash
##
## 02_batch_wb_historical.sbatch - SLURM Batch Script for Historical Water Balance
##
## Runs the water balance model for all sites using GridMET historical climate
## data (1979-2024). Uses GNU Parallel to process multiple sites concurrently.
##
## USAGE:
##   sbatch src/02_batch_wb_historical.sbatch
##
## INPUTS:
##   src/sites.csv - List of sites to process (reads site names from column 1)
##   For each site: data/input/{site}/ must contain prepared inputs
##
## OUTPUTS:
##   output/{site}/wb/*_historical_gridmet_*.nc - Daily AET/Deficit NetCDFs
##   wb-historical-{jobid}.out/.err - SLURM logs
##
## RESOURCE REQUIREMENTS:
##   16 CPUs, 128GB RAM, 1 day runtime
##   Runs 6 sites in parallel (-j6)
##
## PREREQUISITES:
##   - Mamba environment "nps-wb" with Python dependencies
##   - Completed Phase 1 (00_prep_data.R) and Phase 2 (01_clim_data.R) for all sites
##
#SBATCH --account=priority-briansmithers        #specify the account to use
#SBATCH --job-name=wb-historical             # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
##SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16              # number of cores to allocate
#SBATCH --mem=128G                     # 2000 MB of Memory allocated; set --mem with care
##SBATCH --mem-per-cpu=100G
#SBATCH --time=1-00:00:00                 # Maximum job run time
##SBATCH --array=0-1847%10                  # Number of jobs in array
#SBATCH --mail-user=shuysman@gmail.com
#SBATCH --mail-type=ALL
##SBATCH --output=fire-%A-%a.out
##SBATCH --error=fire-%A-%a.err
#SBATCH --output=wb-historical-%j.out
#SBATCH --error=wb-historical-%j.err

date
hostname -s
module load Mamba
source $HOME/.bashrc
mamba activate nps-wb
sites=$(awk -F, 'NR > 1 { print $1 }' sites.csv)
parallel -j6 python 02_start_wb_v_1_5.py {1} {2} {3} :::  historical ::: gridmet ::: $sites
date
