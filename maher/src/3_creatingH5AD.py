#!/usr/bin/env python

from pathlib import Path
import scanpy as sc
import anndata as ad

# ============================================================
# Paths
# ============================================================

OUTDIR = Path("/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_maher/2_cellranger/")

CELLRANGER_OUTDIR = OUTDIR

# Final location for combined AnnData object
FINAL_H5AD_OUTDIR = Path("/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_maher/3_analysis/00_AnnData")

FINAL_H5AD_OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Samples
# ============================================================

samples = [
"1017_G3_D2",
"1017_NC_D1",
"4001_G3_D4",
"4001_NC_D3",
"4002_G3-1_D6",
"4002_G3-2_D7",
"4002_G3-3_D8",
"4002_NC_D5"
]


# ============================================================
# Sample metadata
# ============================================================

# Fill these in according to your experiment

sample_metadata = {
    "1017_G3_D2": {
        "cell_line": "1017_G3_D2",
        "genome": "otter", 
        "genotype": "KO",
	    "pair": "stem_1017BM", 
    },
    "1017_NC_D1": {
        "cell_line": "1017_NC_D1", 
        "genome": "Koala2", 
        "genotype": "WT",
        "pair": "stem_1017BM",
    },
    "4001_G3_D4": {
        "cell_line": "4001_G3_D4", 
        "genome": "lemur", 
        "genotype": "KO",
        "pair": "stem_4001BM",
    },
    "4001_NC_D3": {
        "cell_line": "4001_NC_D3", 
        "genome": "whale", 
        "genotype": "WT",
        "pair": "stem_4001BM",
    },
    "4002_G3-1_D6": {
        "cell_line": "4002_G3-1_D6", 
        "genome": "zebra1", 
        "genotype": "KO",
        "pair": "stem_4002BM",
    },
    "4002_G3-2_D7": {
        "cell_line": "4002_G3-2_D7", 
        "genome": "zebra2", 
        "genotype": "KO",
        "pair": "stem_4002BM",
    },
    "4002_G3-3_D8": {
        "cell_line": "4002_G3-3_D8", 
        "genome": "zebra3", 
        "genotype": "KO",
        "pair": "stem_4002BM",
    },
    "4002_NC_D5": {
        "cell_line": "4002_NC_D5", 
        "genome": "sloth", 
        "genotype": "WT",
        "pair": "stem_4002BM",
    }
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
    label="batch",
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

print("\nbatch counts:")
print(adata.obs["batch"].value_counts())
# ============================================================
# Save combined AnnData object
# ============================================================

FINAL_H5AD = (FINAL_H5AD_OUTDIR / "organoid_scRNA_combined_BM.h5ad")

print(f"\nSaving AnnData to:")
print(FINAL_H5AD)

adata.write_h5ad(FINAL_H5AD)

print("\nDone!")
