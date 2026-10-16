#!/bin/bash

#SBATCH --job-name=cellranger
#SBATCH --array=0-5
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --output=logs/2_cellranger_%A_%a.out
#SBATCH --error=logs/2_cellranger_%A_%a.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166

module purge
module load conda_R/4.5.x
module list
conda activate aanri
ml cellranger/10.1.0 

set -euo pipefail

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

# Create a for loop with directory names and sample names and run 
FASTQ_ROOT="/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/0_raw_data/AN00030578_10X_RawData_Outs" 
REFERENCE="/dcs04/lieber/hwanglab/Arun/bmaher/Annotations/refdata-gex-GRCh38-2024-A"
OUTDIR="/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/2_cellranger"

mkdir -p "${OUTDIR}"

cd ${OUTDIR}

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

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

FASTQ_DIR="${FASTQ_ROOT}/${sample}/23VTCWLT4"


# ---------------------------------------------------------
# Print information
# ---------------------------------------------------------

echo "=========================================="
echo "Sample:     ${sample}"
echo "FASTQ dir:  ${FASTQ_DIR}"
echo "Reference:  ${REFERENCE}"
echo "Output:     ${OUTDIR}"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Array ID:   ${SLURM_ARRAY_TASK_ID}"
echo "=========================================="


# ---------------------------------------------------------
# Check FASTQ directory
# ---------------------------------------------------------

if [[ ! -d "${FASTQ_DIR}" ]]; then
    echo "ERROR: FASTQ directory does not exist:"
    echo "${FASTQ_DIR}"
    exit 1
fi

echo "FASTQ files:"
find "${FASTQ_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" -print


# ---------------------------------------------------------
# Run Cell Ranger
# ---------------------------------------------------------

cellranger count \
    --id="${sample}" \
    --transcriptome="${REFERENCE}" \
    --fastqs="${FASTQ_DIR}" \
    --sample="${sample}" \
    --localcores="${SLURM_CPUS_PER_TASK}" \
    --localmem=60 \
    --create-bam=true


echo "=========================================="
echo "Finished ${sample}"
echo "=========================================="
