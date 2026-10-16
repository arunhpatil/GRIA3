#!/bin/bash

#SBATCH --job-name=maher_harmonize
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/4_harmonize_%j.out
#SBATCH --error=logs/4_harmonize_%j.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166

module purge
module load conda_R/4.5.x
module list
conda activate hnoca_gpu
#conda activate aanri

set -euo pipefail

python 4_QC_harmony.py 
