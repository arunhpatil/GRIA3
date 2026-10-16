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

## Note: This requires us to have gene_qc metics computed!.
python QC_threshold_optimization.py --adata /dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_maher/3_analysis/01_qc_normalization/1.3_organoid_harmonized_BM.h5ad --sample-key sample --output-dir /dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_maher/3_analysis/01_qc_normalization/QC_plots/OPTIMIZATION
