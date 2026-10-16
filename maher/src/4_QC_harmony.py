#!/usr/bin/env python

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import harmonypy as hm

# ============================================================
# OPTIONAL DEPENDENCIES
# ============================================================
import soupx
import scrublet as scr
import scvi

# ============================================================
# SETTINGS
# ============================================================

INPUT_H5AD = Path(
    "/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_maher/3_analysis/"
    "00_AnnData/organoid_scRNA_combined_BM.h5ad"
)

OUTPUT_DIR = Path(
    "/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_maher/3_analysis/"
    "01_qc_normalization"
)

PLOT_DIR = OUTPUT_DIR / "QC_plots"
# SOUPX_DIR = OUTPUT_DIR / "soupx"
# DOUBLET_DIR = OUTPUT_DIR / "doublets"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)
# SOUPX_DIR.mkdir(parents=True, exist_ok=True)
# DOUBLET_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# RAW 10X DIRECTORIES FOR SOUPX
#
# IMPORTANT:
# Each path must point to the Cell Ranger raw_feature_bc_matrix
# directory for that sample, containing matrix.mtx(.gz),
# barcodes.tsv(.gz), and features.tsv(.gz).
#
# Replace the paths below with your actual sample directories.
# ============================================================

# samples = ["1017G3", "1017NC", "4001G3", "4001NC", "4002G3", "4002NC"]
samples = ["1017_G3_D2", "1017_NC_D1", "4001_G3_D4", "4001_NC_D3", "4002_G3-1_D6", "4002_G3-2_D7", "4002_G3-3_D8", "4002_NC_D5"]

# str(Path("/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_maher/2_cellranger/") / f"{samples[0]}/outs/raw_feature_bc_matrix")
bc_path_list = [Path("/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_maher/2_cellranger/") / f"{i}/outs/raw_feature_bc_matrix"  for i in samples]

RAW_10X_DIRS = dict(zip(samples, bc_path_list))

# ============================================================
# PARAMETERS
# ============================================================

# SoupX
SOUPX_METHOD = "multinomial"
SOUPX_ROUND_TO_INT = True

# Minimum number of cells in which a gene must be detected
MIN_CELLS_PER_GENE = 5

# Minimum number of detected genes per cell
MIN_GENES = 500 

# IQR filtering
IQR_K = 1.5

# Scrublet
SCRUBLET_EXPECTED_DOUBLET_RATE = 0.06

# SOLO fallback
SOLO_MAX_EPOCHS = 400
SOLO_BATCH_SIZE = 128
SOLO_DOUBLET_RATIO = 2

# HVGs
N_HVG = 3000

# PCA / neighbors / UMAP
N_PCS = 50
N_UMAP_PCS = 20
N_NEIGHBORS = 15

# Harmony
HARMONY_KEY = "cell_line"

# Leiden
LEIDEN_RESOLUTION = 0.5


# ============================================================
# HELPERS
# ============================================================

def save_plot(filename):
    plt.savefig(PLOT_DIR / f"{filename}.pdf", bbox_inches="tight")
    plt.savefig(PLOT_DIR / f"{filename}.png", dpi=300, bbox_inches="tight")
    plt.close()


def add_qc_annotations(adata):
    adata.var["mt"] = (
        adata.var_names.str.upper().str.startswith("MT-")
    )
    adata.var["ribo"] = (
        adata.var_names.str.upper().str.startswith("RPS")
        | adata.var_names.str.upper().str.startswith("RPL")
    )
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        inplace=True,
    )
    return adata


def iqr_mask(values, k=1.5):
    values = np.asarray(values, dtype=float)
    q1 = np.nanpercentile(values, 25)
    q3 = np.nanpercentile(values, 75)
    iqr = q3 - q1
    lower = q1 - k * iqr
    upper = q3 + k * iqr
    return (values >= lower) & (values <= upper), lower, upper


def per_sample_iqr_filter(adata, sample_key="sample", k=1.5):
    """
    Apply IQR filtering separately within each sample to:
      1. log1p(total_counts)
      2. pct_counts_mt

    A cell is retained only if it passes both filters.
    """
    if sample_key not in adata.obs.columns:
        raise KeyError(f"{sample_key!r} is not present in adata.obs")
    keep = np.zeros(adata.n_obs, dtype=bool)
    print("\nPer-sample IQR filtering:")
    print("----------------------------------------------")
    for sample in adata.obs[sample_key].astype(str).unique():
        idx = np.asarray(adata.obs[sample_key].astype(str) == sample)
        log_counts = np.log1p(
            adata.obs.loc[idx, "total_counts"].to_numpy()
        )
        mt = adata.obs.loc[idx, "pct_counts_mt"].to_numpy()
        keep_counts, lo_counts, hi_counts = iqr_mask(log_counts, k)
        keep_mt, lo_mt, hi_mt = iqr_mask(mt, k)
        sample_keep = keep_counts & keep_mt
        keep[idx] = sample_keep
        print(
            f"{sample}: "
            f"{idx.sum():,} cells -> {sample_keep.sum():,} cells | "
            f"log1p(counts): [{lo_counts:.2f}, {hi_counts:.2f}] | "
            f"mt: [{lo_mt:.2f}, {hi_mt:.2f}]"
        )
    return keep


def run_soupx_sample(sample_adata, sample_name, raw_dir):
    """
    Run SoupX on one sample.
    The filtered matrix comes from sample_adata.
    The raw droplets come from the Cell Ranger raw_feature_bc_matrix.
    """
    if soupx is None:
        raise ImportError(
            "soupx is not installed. Install a SoupX Python implementation "
            "such as soupx-python or pysoupx."
        )
    raw_dir = Path(raw_dir)
    if not raw_dir.exists():
        raise FileNotFoundError(
            f"Raw 10X directory does not exist for {sample_name}: {raw_dir}"
        )
    print(f"\nRunning SoupX: {sample_name}")
    print(f"Raw 10X: {raw_dir}")
    # soupx-python API.
    # The raw matrix is genes x droplets and the filtered matrix is
    # genes x cells.
    raw_adata = sc.read_10x_mtx(
        raw_dir,
        var_names="gene_symbols",
        make_unique=True,
    )
    raw_adata.var_names_make_unique()
    sample_adata.var_names_make_unique()
    common_genes = sample_adata.var_names[
        sample_adata.var_names.isin(raw_adata.var_names)
    ]
    if len(common_genes) == 0:
        raise ValueError(
            f"No overlapping genes between raw and filtered matrices "
            f"for {sample_name}."
        )
    filtered = sample_adata[:, common_genes].copy()
    raw = raw_adata[:, common_genes].copy()
    soup_channel = soupx.SoupChannel(
        tod=raw.X.T.tocsr(),
        toc=filtered.X.T.tocsr(),
        metaData=pd.DataFrame(index=filtered.obs_names),
    )
    # autoEstCont() requires clustering information for the preferred
    # cluster-level contamination estimation.
    #
    # These preliminary clusters are ONLY for SoupX contamination
    # estimation. They are not the final integrated clusters.
    tmp = filtered.copy()
    if "counts" in tmp.layers:
        tmp.X = tmp.layers["counts"].copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(
        tmp,
        n_top_genes=min(2000, tmp.n_vars),
        flavor="seurat",
    )
    if tmp.n_obs >= 50 and tmp.n_vars >= 50:
        sc.pp.scale(tmp, max_value=10)
        n_pcs = min(20, tmp.n_vars - 1, tmp.n_obs - 1)
        if n_pcs >= 2:
            sc.tl.pca(tmp, n_comps=n_pcs, use_highly_variable=True, svd_solver="arpack",)
            sc.pp.neighbors(tmp, n_neighbors=min(15, tmp.n_obs - 1))
            sc.tl.leiden(tmp, resolution=0.5, key_added="soupx_cluster")
            clusters = tmp.obs["soupx_cluster"].astype(str).to_numpy()
            soup_channel.setClusters(clusters)
    # Automatic contamination estimation.
    soup_channel = soupx.autoEstCont(
        soup_channel,
        verbose=True,
    )
    # Multinomial correction as requested.
    corrected = soupx.adjustCounts(
        soup_channel,
        method=SOUPX_METHOD,
        roundToInt=SOUPX_ROUND_TO_INT,
    )
    corrected = corrected.T.tocsr()
    corrected_adata = filtered.copy()
    corrected_adata.X = corrected
    if hasattr(soup_channel, "metaData"):
        if "rho" in soup_channel.metaData.columns:
            rho = soup_channel.metaData["rho"]
            corrected_adata.obs["soupx_contamination"] = rho.reindex(
                corrected_adata.obs_names
            ).to_numpy()
    corrected_adata.obs["soupx_sample"] = sample_name
    return corrected_adata


def run_scrublet(sample_adata, sample_name):
    """
    Run Scrublet on one sample.
    Returns:
        mask, score, threshold
    """
    if scr is None:
        raise ImportError("scrublet is not installed.")
    counts = sample_adata.X
    scrub = scr.Scrublet(counts, expected_doublet_rate=SCRUBLET_EXPECTED_DOUBLET_RATE,)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    threshold = scrub.threshold_
    return (np.asarray(predicted_doublets, dtype=bool), np.asarray(doublet_scores, dtype=float), float(threshold))


def run_solo(sample_adata, sample_name):
    """
    SOLO fallback for one sample.

    SOLO is trained from a sample-specific SCVI model because SOLO
    should be run one batch/sample at a time.
    """
    if scvi is None:
        raise ImportError("scvi-tools is not installed.")
    solo_adata = sample_adata.copy()
    # Ensure integer count data.
    if "counts" in solo_adata.layers:
        solo_adata.X = solo_adata.layers["counts"].copy()
    scvi.model.SCVI.setup_anndata(
        solo_adata,
        batch_key=None,
    )
    vae = scvi.model.SCVI(
        solo_adata,
        n_latent=20,
    )
    vae.train(
        max_epochs=100,
        early_stopping=True,
    )
    solo = scvi.external.SOLO.from_scvi_model(
        vae,
        doublet_ratio=SOLO_DOUBLET_RATIO,
    )
    solo.train(
        max_epochs=SOLO_MAX_EPOCHS,
        batch_size=SOLO_BATCH_SIZE,
        early_stopping=True,
    )
    pred = solo.predict(
        soft=False,
        include_simulated_doublets=False,
    )
    # scvi-tools SOLO prediction column is typically "prediction".
    if "prediction" in pred.columns:
        labels = pred["prediction"].astype(str)
    else:
        labels = pred.iloc[:, 0].astype(str)
    doublet_mask = labels.str.lower().eq("doublet").to_numpy()
    if "doublet" in pred.columns:
        score = pred["doublet"].to_numpy()
    else:
        score = np.full(sample_adata.n_obs, np.nan)
    return doublet_mask, score


# ============================================================
# LOAD DATA
# ============================================================

print("==============================================")
print("Loading AnnData")
print("==============================================")

adata = sc.read_h5ad(INPUT_H5AD)

print(adata)
print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")
print("\nMetadata:")
print(adata.obs.columns.tolist())

if "sample" not in adata.obs.columns:
    raise KeyError(
        "The AnnData object must contain a 'sample' column for "
        "per-sample SoupX/QC/doublet processing."
    )


# ============================================================
# PRESERVE ORIGINAL COUNTS
# ============================================================

if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

print("\nOriginal count layer preserved as adata.layers['counts'].")


# ============================================================
# QC GENE ANNOTATION
# ============================================================
add_qc_annotations(adata)

# ============================================================
# AMBIENT RNA CORRECTION WITH SOUPX
# ============================================================

print("\n==============================================")
print("Ambient RNA correction with SoupX")
print("==============================================")

if not RAW_10X_DIRS:
    raise RuntimeError(
        "RAW_10X_DIRS is empty. True SoupX autoEstCont() requires "
        "the raw droplet matrix for each sample. Add the Cell Ranger "
        "raw_feature_bc_matrix path for every sample before running."
    )

samples_in_adata = adata.obs["sample"].astype(str).unique()

missing_raw = [
    sample for sample in samples_in_adata
    if sample not in RAW_10X_DIRS
]

if missing_raw:
    raise RuntimeError(
        "Missing RAW_10X_DIRS entries for samples: "
        + ", ".join(missing_raw)
    )

soupx_objects = []

for sample in samples_in_adata:
    sample_mask = adata.obs["sample"].astype(str) == sample
    sample_adata = adata[sample_mask].copy()
    # Use the preserved raw counts as the filtered Cell Ranger matrix.
    sample_adata.X = sample_adata.layers["counts"].copy()
    corrected = run_soupx_sample(
        sample_adata,
        sample,
        RAW_10X_DIRS[sample],
    )
    soupx_objects.append(corrected)

adata = sc.concat(
    soupx_objects,
    axis=0,
    join="inner",
    merge="same",
    uns_merge="unique",
    index_unique=None,
)

adata.var_names_make_unique()

print(
    f"\nAfter SoupX: {adata.n_obs:,} cells x "
    f"{adata.n_vars:,} genes"
)

SOUPX_H5AD = OUTPUT_DIR / "1.1_organoid_SoupX_corrected_BM.h5ad"
adata.write_h5ad(SOUPX_H5AD)


# ============================================================
# QC AFTER SOUPX
# ============================================================

print("\n==============================================")
print("QC after SoupX")
print("==============================================")

add_qc_annotations(adata)

print("\nQC summary:")
print(
    adata.obs[
        [
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
        ]
    ].describe()
)

if "sample" in adata.obs.columns:
    sc.pl.violin(
        adata,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"],
        groupby="sample",
        rotation=90,
        multi_panel=True,
        stripplot=False,
        show=False,
    )
    save_plot("1_QC_after_SoupX_by_sample")


# ============================================================
# BASIC CELL/GENE FILTERING
# ============================================================

print("\n==============================================")
print("Minimum gene/cell filtering")
print("==============================================")

cells_before = adata.n_obs
genes_before = adata.n_vars

sc.pp.filter_genes(
    adata,
    min_cells=MIN_CELLS_PER_GENE,
)

# Recalculate cell QC after gene filtering.
add_qc_annotations(adata)

adata = adata[
    adata.obs["n_genes_by_counts"] >= MIN_GENES
].copy()

print(
    f"Cells after minimum gene filter: "
    f"{cells_before:,} -> {adata.n_obs:,}"
)

print(
    f"Genes after min_cells={MIN_CELLS_PER_GENE}: "
    f"{genes_before:,} -> {adata.n_vars:,}"
)


# ============================================================
# PER-SAMPLE IQR QC
# ============================================================

print("\n==============================================")
print("Per-sample IQR QC")
print("==============================================")

iqr_keep = per_sample_iqr_filter(
    adata,
    sample_key="sample",
    k=IQR_K,
)

cells_before_iqr = adata.n_obs
adata = adata[iqr_keep].copy()

print(
    f"\nCells after per-sample IQR filtering: "
    f"{cells_before_iqr:,} -> {adata.n_obs:,}"
)

add_qc_annotations(adata)

# QC_H5AD = OUTPUT_DIR / "2.1_organoid_QC_filtered_BM.h5ad"
# adata.write_h5ad(QC_H5AD)


# ============================================================
# DOUBLETS: SCRUBLET FIRST, SOLO FALLBACK
# ============================================================

print("\n==============================================")
print("Doublet detection")
print("Scrublet first; SOLO fallback if Scrublet fails")
print("==============================================")

adata.obs["doublet_method"] = "not_run"
adata.obs["doublet_score"] = np.nan
adata.obs["predicted_doublet"] = False

doublet_objects = []

for sample in adata.obs["sample"].astype(str).unique():
    sample_mask = adata.obs["sample"].astype(str) == sample
    sample_adata = adata[sample_mask].copy()
    
    # Doublet detection must use counts.
    if "counts" in sample_adata.layers:
        sample_adata.X = sample_adata.layers["counts"].copy()
    
    print(f"\nDoublet detection: {sample}")
    
    try:
        doublet_mask, scores, threshold = run_scrublet(
            sample_adata,
            sample,
        )
        
        print(
            f"Scrublet successful: "
            f"{doublet_mask.sum():,} doublets / "
            f"{sample_adata.n_obs:,} cells "
            f"(threshold={threshold:.4f})"
        )
        
        sample_adata.obs["doublet_method"] = "Scrublet"
        sample_adata.obs["doublet_score"] = scores
        sample_adata.obs["predicted_doublet"] = doublet_mask
    except Exception as e:
        warnings.warn(
            f"Scrublet failed for {sample}: {e}. "
            "Falling back to SOLO."
        )
        doublet_mask, scores = run_solo(
            sample_adata,
            sample,
        )
        
        print(
            f"SOLO successful: "
            f"{doublet_mask.sum():,} doublets / "
            f"{sample_adata.n_obs:,} cells"
        )
        
        sample_adata.obs["doublet_method"] = "SOLO"
        sample_adata.obs["doublet_score"] = scores
        sample_adata.obs["predicted_doublet"] = doublet_mask
    
    doublet_objects.append(sample_adata)

adata = sc.concat(
    doublet_objects,
    axis=0,
    join="inner",
    merge="same",
    uns_merge="unique",
    index_unique=None,
)

print(
    f"\nTotal predicted doublets: "
    f"{adata.obs['predicted_doublet'].sum():,}"
)

sc.pl.violin(
    adata,
    "doublet_score",
    groupby="sample",
    rotation=90,
    stripplot=False,
    show=False,
)
save_plot("2_Doublet_scores_by_sample")


# ------------------------------------------------------------
# Remove predicted doublets
# ------------------------------------------------------------

cells_before_doublet = adata.n_obs

adata = adata[
    ~adata.obs["predicted_doublet"]
].copy()

print(
    f"Cells after doublet removal: "
    f"{cells_before_doublet:,} -> {adata.n_obs:,}"
)

DOUBLETS_H5AD = OUTPUT_DIR / "1.2_organoid_doublet_filtered_BM.h5ad"
adata.write_h5ad(DOUBLETS_H5AD)


# ============================================================
# PRESERVE CORRECTED COUNTS
# ============================================================

adata.layers["counts"] = adata.X.copy()


# ============================================================
# NORMALIZATION
# ============================================================

print("\n==============================================")
print("Normalization")
print("==============================================")

sc.pp.normalize_total(
    adata,
    target_sum=10000,
)

sc.pp.log1p(adata)

adata.raw = adata


# ============================================================
# CELL CYCLE SCORING
# ============================================================

s_genes = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1",
    "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1",
    "E2F8",
]

g2m_genes = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "CENPF", "SMC4", "CCNB2", "CKS2", "NUF2", "CKS1B", "MKI67",
    "TMPO", "CENPA", "TACC3", "FAM64A", "SMC2", "CCNB1", "CDC20",
    "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5",
    "CDCA3", "HMMR", "AURKB", "KIF20B", "CTCF", "NEK2", "G2E3",
    "GAS2L3", "CBX5", "CENPE", "AURKA", "GTSE1", "KIF11", "ANLN",
    "UBE2S", "BUB1", "CENPM", "CDC25A", "KIF23", "PIM1", "KIF4A",
    "ECT2", "KIF18B", "MCM10", "TK1",
]

s_genes_present = [g for g in s_genes if g in adata.var_names]
g2m_genes_present = [g for g in g2m_genes if g in adata.var_names]

print(f"S-phase genes present: {len(s_genes_present)}")
print(f"G2M genes present: {len(g2m_genes_present)}")

if len(s_genes_present) > 0 and len(g2m_genes_present) > 0:
    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes_present,
        g2m_genes=g2m_genes_present,
    )
    print("Cell-cycle scoring complete.")


# ============================================================
# HIGHLY VARIABLE GENES
# ============================================================

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_HVG,
    flavor="seurat_v3",
    layer="counts",
)

print(
    "Number of HVGs:",
    int(adata.var["highly_variable"].sum()),
)

sc.pl.highly_variable_genes(
    adata,
    show=False,
)
save_plot("3_highly_variable_genes")


# ============================================================
# SCALE DATA
# ============================================================

print("\nScaling data...")
sc.pp.scale(
    adata,
    max_value=10,
)


# ============================================================
# PCA
# ============================================================

print("\n==============================================")
print("PCA")
print("==============================================")

sc.tl.pca(
    adata,
    n_comps=N_PCS,
    use_highly_variable=True,
    svd_solver="arpack",
)

if "sample" in adata.obs.columns:
    sc.pl.pca(
        adata,
        color="sample",
        show=False,
    )
    save_plot("4_PCA_by_sample")

if "genotype" in adata.obs.columns:
    sc.pl.pca(
        adata,
        color="genotype",
        show=False,
    )
    save_plot("5_PCA_by_genotype")

if "pair" in adata.obs.columns:
    sc.pl.pca(
        adata,
        color="pair",
        show=False,
    )
    save_plot("5_PCA_by_pair")

# ============================================================
# HARMONY INTEGRATION
# ============================================================

print("\n==============================================")
print("Harmony integration")
print("==============================================")

if HARMONY_KEY not in adata.obs.columns:
    raise KeyError(
        f"Harmony key {HARMONY_KEY!r} not found in adata.obs."
    )

pca = adata.obsm["X_pca"]

ho = hm.run_harmony(
    adata.obsm["X_pca"],
    adata.obs,
    HARMONY_KEY,
)

harmony = ho.Z_corr

print("PCA shape:", pca.shape)
print("Raw Harmony shape:", harmony.shape)

if harmony.shape == pca.shape:
    adata.obsm["X_pca_harmony"] = harmony
elif harmony.T.shape == pca.shape:
    adata.obsm["X_pca_harmony"] = harmony.T
else:
    raise ValueError(
        f"Unexpected Harmony shape: {harmony.shape}; "
        f"expected {pca.shape}"
    )

print(
    "X_pca_harmony shape:",
    adata.obsm["X_pca_harmony"].shape,
)


# ============================================================
# NEIGHBORS ON HARMONY LATENT SPACE
# ============================================================

sc.pp.neighbors(
    adata,
    use_rep="X_pca_harmony",
    n_neighbors=N_NEIGHBORS,
    n_pcs=N_UMAP_PCS,
)


# ============================================================
# UMAP
# ============================================================

sc.tl.umap(adata)


# ============================================================
# LEIDEN CLUSTERING ON HARMONY LATENT SPACE
# ============================================================

print("\n==============================================")
print("Leiden clustering")
print("==============================================")

sc.tl.leiden(
    adata,
    resolution=LEIDEN_RESOLUTION,
    key_added="leiden",
)

print(
    "Number of clusters:",
    adata.obs["leiden"].nunique(),
)


# ============================================================
# INTEGRATED UMAP PLOTS
# ============================================================

sc.pl.umap(
    adata,
    color=[
        "leiden",
        "cell_line",
        "genotype",
        "pair",
        "genome",
    ],
    ncols=2,
    show=False,
)
save_plot("6_UMAP_integrated_combined")

if "sample" in adata.obs.columns:
    sc.pl.umap(
        adata,
        color="sample",
        show=False,
    )
    save_plot("7_UMAP_by_sample_integrated")

if "genotype" in adata.obs.columns:
    sc.pl.umap(
        adata,
        color="genotype",
        show=False,
    )
    save_plot("8_UMAP_by_genotype_integrated")

if "genome" in adata.obs.columns:
    sc.pl.umap(
        adata,
        color="genome",
        show=False,
    )
    save_plot("9_UMAP_by_genome_integrated")

if "pair" in adata.obs.columns:
    sc.pl.umap(
        adata,
        color="pair",
        show=False,
    )
    save_plot("10_UMAP_by_pair_integrated")

sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    show=False,
)
save_plot("11_UMAP_by_LeidenCluster_integrated")


# ============================================================
# GRIA3 EXPRESSION
# ============================================================

if "GRIA3" in adata.var_names:
    sc.pl.umap(
        adata,
        color="GRIA3",
        use_raw=True,
        show=False,
    )
    save_plot("12_UMAP_GRIA3")
else:
    print("WARNING: GRIA3 not found.")


if "GRIA3" in adata.var_names and "genotype" in adata.obs:
  sc.pl.umap(
      adata,
      color=["GRIA3", "genotype"],
      use_raw=True,
      show=False,
  )
  save_plot(
      "12_UMAP_GRIA3_vs_genotype.jpeg")
else:
  print("WARNING: 'GRIA3' not found in var_names or 'genotype' not found in obs.")


# ============================================================
# SAVE FINAL OBJECT
# ============================================================

FINAL_H5AD = OUTPUT_DIR / "1.3_organoid_harmonized_BM.h5ad"

adata.write_h5ad(FINAL_H5AD)

# ============================================================
# SUMMARY
# ============================================================

print("\n==============================================")
print("FINAL SUMMARY")
print("==============================================")

print(f"Final cells : {adata.n_obs:,}")
print(f"Final genes : {adata.n_vars:,}")
print(
    f"HVGs        : "
    f"{int(adata.var['highly_variable'].sum()):,}"
)

print("\nOutputs:")
print(OUTPUT_DIR)

print("\nPlots:")
print(PLOT_DIR)

print("\nFinal H5AD:")
print(FINAL_H5AD)

print("\nDone.")
