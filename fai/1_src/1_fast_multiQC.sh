#!/bin/bash

#SBATCH --job-name=fai_fastq_qc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/1_fastqc_%j.out
#SBATCH --error=logs/1_fastqc_%j.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166

module purge
module load conda_R/4.5.x
module list
conda activate aanri
ml cellranger/10.1.0 
ml fastqc/0.12.1

set -euo pipefail

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

# Create a for loop with directory names and sample names and run 
FASTQ_ROOT="/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/0_raw_data/AN00030578_10X_RawData_Outs" 
REFERENCE="/dcs04/lieber/hwanglab/Arun/bmaher/Annotations/refdata-gex-GRCh38-2024-A"
QC_DIR="/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/1_qualityCtrl"
MULTIQC_DIR="${QC_DIR}/multiqc"

mkdir -p "${QC_DIR}"
mkdir -p "${MULTIQC_DIR}"

# ---------------------------------------------------------
# Samples
# ---------------------------------------------------------

samples=(
    "1017G3"
    "1017NC"
    "4001G3"
    "4001NC"
    "4002G3"
    "4002NC"
)

for sample in "${samples[@]}"; do

    FASTQ_DIR="${FASTQ_ROOT}/${sample}/23VTCWLT4"

    echo "======================================"
    echo "FastQC: ${sample}"
    echo "FASTQ:  ${FASTQ_DIR}"
    echo "======================================"

    if [[ ! -d "${FASTQ_DIR}" ]]; then
        echo "ERROR: Directory does not exist:"
        echo "${FASTQ_DIR}"
        exit 1
    fi

    fastqc \
        "${FASTQ_DIR}"/*.fastq.gz \
        --outdir="${QC_DIR}" \
        --threads="${SLURM_CPUS_PER_TASK}"

done

echo "Running MultiQC..."

multiqc \
    "${QC_DIR}" \
    --outdir="${MULTIQC_DIR}" \
    --filename="multiqc_report.html"

echo "=========================================="
echo "QC completed."
echo "Finished ${sample}"
echo "=========================================="
