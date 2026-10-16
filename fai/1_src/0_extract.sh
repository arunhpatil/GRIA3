#!/bin/bash
#SBATCH --job-name=extract
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --output=logs/extract.out
#SBATCH --error=logs/extract.err
#SBATCH --mem=10G

cd /dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/0_raw_data
tar -xf AN00030578_10X_RawData_Outs.tar

