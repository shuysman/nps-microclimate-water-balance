#!/usr/bin/env bash
#
# 02_run_all_historical_projections.sh - Run All Water Balance Scenarios (Local)
#
# Runs the water balance model for all sites across historical and projected
# climate scenarios. Use this for local/workstation runs; for HPC, use the
# SLURM batch scripts instead.
#
# SCENARIOS RUN:
#   historical/gridmet  - 1979-2024 observed climate
#   CanESM2/rcp85       - Hot-dry future (worst case)
#   HadGEM2-CC365/rcp85 - Hot future
#   MRI-CGCM3/rcp45     - Warm-wet future (best case)
#   MRI-CGCM3/rcp85     - Moderate warming
#
# USAGE:
#   bash src/02_run_all_historical_projections.sh
#
# INPUTS:
#   src/sites.csv - List of sites to process
#
# OUTPUTS:
#   output/{site}/wb/ - Daily water balance NetCDFs for each scenario
#
# NOTE: Requires GNU Parallel and significant compute resources.
#       Runs 32 jobs concurrently (-j32).
#
# GCM/RCP combinations chosen to represent extremes in area of interest, following
# Lawrence, David J., Amber N. Runyon, John E. Gross, Gregor W. Schuurman, and Brian W. Miller. 2021.
# "Divergent, Plausible, and Relevant Climate Futures for near- and Long-Term Resource Planning."
# Climatic Change 167 (3): 38. https://doi.org/10.1007/s10584-021-03169-y.
#
set -euo pipefail

sites=$(awk -F, 'NR > 1 { print $1 }' src/sites.csv)

parallel -j32 python src/02_start_wb_v_1_5.py {1} {2} {3} \
  ::: historical CanESM2 HadGEM2-CC365 MRI-CGCM3 MRI-CGCM3 \
  :::+ gridmet rcp85 rcp85 rcp45 rcp85 \
  ::: $sites
g
