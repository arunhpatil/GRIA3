#!/bin/bash

#SBATCH --job-name=fai_scanvi
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --time=4:00:00
#SBATCH --output=logs/5_annotation_%j.out
#SBATCH --error=logs/5_annotation_%j.err
#SBATCH --cpus-per-task=5
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166


module purge
module load conda_R/4.5.x
module list
conda activate hnoca_gpu

set -euo pipefail

#quarto render 5_cellAnnotation.ipynb --execute -M jupyter:kernel=python3
#quarto render 5_cellAnnotation.ipynb --execute

jupyter nbconvert --to html --execute 5_cellAnnotation.ipynb

#jupyter nbconvert --to html 5_cellAnnotation.ipynb && \
#jupyter nbconvert --to pdf 5_cellAnnotation.ipynb

###SBATCH --gres=gpu:l40s:1
