#!/bin/bash

#SBATCH --job-name=combineH5
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --output=logs/3_cellranger_combineH5_%j.out
#SBATCH --error=logs/3_cellranger_combineH5_%j.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166

module purge
module load conda_R/4.5.x
module list
conda activate aanri
ml cellranger/10.1.0 

set -euo pipefail

python /dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/source/3_creatingH5AD.py
