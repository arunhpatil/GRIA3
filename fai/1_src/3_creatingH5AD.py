#!/usr/bin/env python

from pathlib import Path
import scanpy as sc
import anndata as ad

# ============================================================
# Paths
# ============================================================

OUTDIR = Path("/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/2_cellranger")

# Directory containing the Cell Ranger outputs
# Example:
# OUTDIR/
# ├── Sample1/outs/filtered_feature_bc_matrix.h5
# ├── Sample2/outs/filtered_feature_bc_matrix.h5
# └── ...

CELLRANGER_OUTDIR = OUTDIR

# Final location for combined AnnData object
FINAL_H5AD_OUTDIR = Path("/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/3_analysis/00_AnnData")
#
#├── 3_analysis/
#│   ├── 00_AnnData/
#│   ├── 01_qc/
#│   ├── 02_normalization/
#│   ├── 03_clustering/
#│   ├── 04_annotation/
#│   ├── 05_pseudotime/
#│   ├── 06_DE/
#│   ├── 07_TF_activity/
#│   └── 08_GO_pathways/

FINAL_H5AD_OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Samples
# ============================================================

samples = [
"1017G3",
"1017NC",
"4001G3",
"4001NC",
"4002G3",
"4002NC",
]


# ============================================================
# Sample metadata
# ============================================================

# Fill these in according to your experiment

sample_metadata = {
    "1017G3": {
        "cell_line": "1017G3",
        "genome": "sow", # A female otter is a sow
        "genotype": "KO",
	"pair": "stem_1017DF", # Organoids 1017 Danielle Fai (org_1017DF)
    },

    "1017NC": {
        "cell_line": "1017NC", 
        "genome": "doe", # A female koala is a doe
        "genotype": "WT",
        "pair": "stem_1017DF",
    },

    "4001G3": {
        "cell_line": "4001G3", 
        "genome": "troop", # Social group of female Lemurs (there is no name for a female lemur)
        "genotype": "KO",
        "pair": "stem_4001DF",
    },

    "4001NC": {
        "cell_line": "4001NC", 
        "genome": "cow", # A female whale is a cow
        "genotype": "WT",
        "pair": "stem_4001DF",
    },

    "4002G3": {
        "cell_line": "4002G3", 
        "genome": "mare", # A female zebra is a mare
        "genotype": "KO",
        "pair": "stem_4002DF",
    },

    "4002NC": {
        "cell_line": "4002NC", 
        "genome": "slowbie", # A baby sloth is a slowbie (there is no name for female sloth)
        "genotype": "WT",
        "pair": "stem_4002DF",
    },
}


# ============================================================
# Read Cell Ranger H5 files
# ============================================================

adatas = []

for sample in samples:

    print(f"\nReading sample: {sample}")

    # Cell Ranger output
    h5_file = (
        CELLRANGER_OUTDIR
        / sample
        / "outs"
        / "filtered_feature_bc_matrix.h5"
    )

    print(f"Input: {h5_file}")

    if not h5_file.exists():
        raise FileNotFoundError(
            f"Cell Ranger H5 file not found: {h5_file}"
        )

    # Read 10X H5
    adata = sc.read_10x_h5(h5_file)

    # Make gene names unique
    adata.var_names_make_unique()

    print(adata)


    # ========================================================
    # Add sample metadata
    # ========================================================

    adata.obs["sample"] = sample

    adata.obs["cell_line"] = sample_metadata[sample]["cell_line"]

    adata.obs["genome"] = sample_metadata[sample]["genome"]
    
    adata.obs["genotype"] = sample_metadata[sample]["genotype"]

    adata.obs["pair"] = sample_metadata[sample]["pair"]


    # ========================================================
    # Make cell barcodes unique across samples
    # ========================================================

    adata.obs_names = [
        f"{sample}_{barcode}"
        for barcode in adata.obs_names
    ]

    adata.obs_names_make_unique()


    # Store object
    adatas.append(adata)

    print(f"Finished: {sample}")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")


# ============================================================
# Concatenate all samples
# ============================================================

print("\nConcatenating samples...")

adata = ad.concat(
    adatas,
    axis=0,
    join="outer",
    label="sample_batch",
    keys=samples,
    index_unique=None,
)

print(adata)


# ============================================================
# Check metadata
# ============================================================

print("\nSample counts:")
print(adata.obs["sample"].value_counts())

print("\nGenotype counts:")
print(adata.obs["genotype"].value_counts())

print("\nCell-line counts:")
print(adata.obs["cell_line"].value_counts())

print("\nPair counts:")
print(adata.obs["pair"].value_counts())


# ============================================================
# Save combined AnnData object
# ============================================================

FINAL_H5AD = (FINAL_H5AD_OUTDIR / "organoid_scRNA_combined_DF.h5ad")

print(f"\nSaving AnnData to:")
print(FINAL_H5AD)

adata.write_h5ad(FINAL_H5AD)

print("\nDone!")
