#!/bin/bash
# Job name:
#SBATCH --job-name=test
#
# Partition:
#SBATCH --partition=vulcan_gpu
# Wall clock limit:
#SBATCH --time=0:1:00
#

## Run command
./mu_STEM
