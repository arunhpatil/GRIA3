# %%
import os, sys
from pathlib import Path
# This is a workaround for a bug in ipykernel that causes typing_extensions to be imported from the wrong location. 
sys.path = [p for p in sys.path if "ipykernel/py3" not in p] + [p for p in sys.path if "ipykernel/py3" in p]
for m in [m for m in sys.modules if m == "typing_extensions"]:
    del sys.modules[m]

from typing import Literal

import anndata as ad
import lamindb as ln
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
from scipy.sparse import csc_matrix
from scipy.special import logit
from scipy.stats import median_abs_deviation

# %%
QUERY_PATH = Path(
    "/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_maher/3_analysis/"
    "01_qc_normalization/1.3_organoid_harmonized_BM.h5ad"
)
query_adata = sc.read_h5ad(QUERY_PATH)
print(query_adata)

# %%
qc_features = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']

# 4. Generate the diagnostic violin plots grouped by your "sample" column
# multi_panel=True forces them into separate, side-by-side subplots
sc.pl.violin(
    query_adata, 
    keys=qc_features, 
    groupby='sample', 
    multi_panel=True,
    rotation=45,  # Rotates sample names on X-axis for better readability
    stripplot=True,  # Overlays individual points
    jitter=True,  # Adds jitter to prevent overplotting
    size=3,  # Point size
    show=False    # Prevents premature rendering so we can adjust layout
)

# 5. Fine-tune layout adjustments
# plt.tight_layout()
plt.show()


# %%
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

# 1. Run your function to get the boolean mask array
cells_to_keep = per_sample_iqr_filter(query_adata, sample_key="sample", k=1.5)

# 2. SUBSET THE ANNDATA OBJECT (This step applies the changes)
query_adata = query_adata[cells_to_keep].copy()

sc.pl.violin(
    query_adata, 
    keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
    groupby='sample', 
    multi_panel=True,
    rotation=45
)
plt.show()


# %%
# import pandas as pd
# df = pd.read_json("/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/source/Chats/QC_summaize.json")
# print(df.head())
# print(df.to_string())
# print(df.to_markdown(index=False))
# !pip install tabulate


# %% 
adata = sc.read_h5ad(Path(
    "/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/3_analysis/"
    "01_qc_normalization/1.3_organoid_harmonized_DF.h5ad"
))
# Assess organoid maturity first
print("Sequencing quality summary:")
print(adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']].describe())

# # Then choose thresholds:
# if adata.obs['n_genes_by_counts'].median() < 800:
#     print("Early-stage organoids detected. Using conservative MIN_GENES=300")
#     MIN_GENES = 300
# elif adata.obs['n_genes_by_counts'].median() < 1200:
#     print("Mid-stage organoids. Using MIN_GENES=500")
#     MIN_GENES = 500
# else:
#     print("Mature/well-differentiated organoids. Using MIN_GENES=800")
#     MIN_GENES = 800

# print(f"Selected MIN_GENES = {MIN_GENES}")

# Apply filter
cells_before = adata.n_obs
adata = adata[adata.obs['n_genes_by_counts'] >= MIN_GENES].copy()
print(f"Cells after filtering: {cells_before:,} → {adata.n_obs:,}")


# First, check your actual data
print("Current data overview:")
print(adata.obs['n_genes_by_counts'].describe())
print(f"\nMedian genes per sample:")
print(adata.obs.groupby('sample')['n_genes_by_counts'].median())

# Choose MIN_GENES based on median
median_global = adata.obs['n_genes_by_counts'].median()
if median_global > 1500:
    MIN_GENES = 800
elif median_global > 1000:
    MIN_GENES = 600
else:
    MIN_GENES = 400

print(f"\n✓ Recommended MIN_GENES = {MIN_GENES}")

# Apply with tracking
print(f"\nBefore filtering: {adata.n_obs:,} cells")
adata = adata[adata.obs['n_genes_by_counts'] >= MIN_GENES].copy()
print(f"After filtering: {adata.n_obs:,} cells")


# %%
import scipy.stats as stats

print("\nSequencing depth by sample:")
for sample in sorted(adata.obs['sample'].unique()):
    subset = adata.obs[adata.obs['sample'] == sample]
    print(f"\n{sample}:")
    print(f"  Median genes: {subset['n_genes_by_counts'].median():.0f}")
    print(f"  Median counts: {subset['total_counts'].median():.0f}")
    print(f"  % MT: {subset['pct_counts_mt'].mean():.1f}%")
    print(f"  % cells: {len(subset) / len(adata) * 100:.1f}%")

# 2. Check if Harmony is actually working
sc.pl.pca(adata, color='sample', show=False)  # Before Harmony
sc.pl.umap(adata, color='sample', show=False)  # After Harmony
# Are samples still separated? If yes → Harmony failed (might be MIN_GENES issue)

# 3. Check if samples are biologically different (genotype/treatment)
if 'genotype' in adata.obs.columns:
    sc.pl.umap(adata, color=['sample', 'genotype'], show=False)


# %%
# Is it ACTUALLY different, or just lower quality?
print("\n=== Sequencing Platform Check ===")
print("Sample\t\tMedian Genes\tGenes per UMI\tInterpretation")
print("-" * 60)

for sample in sorted(adata.obs['sample'].unique()):
    subset = adata.obs[adata.obs['sample'] == sample]
    median_genes = subset['n_genes_by_counts'].median()
    umi_per_gene = (subset['total_counts'] / subset['n_genes_by_counts']).median()
    
    print(f"{sample}\t{median_genes:.0f}\t\t{umi_per_gene:.1f}\t\t", end="")
    
    if median_genes < 3500:
        print("⚠️ LOW")
    elif median_genes < 4500:
        print("OK")
    else:
        print("✓ HIGH")

# Check cell cycle / mitochondrial stress
# sc.tl.score_genes_cell_cycle(
#     adata,
#     s_genes=[...],  # your S-phase genes
#     g2m_genes=[...]  # your G2M genes
# )

sc.pl.umap(adata, color=['sample', 'S_score', 'G2M_score', 'pct_counts_mt'], show=False)
# If 4001NC has high S_score or G2M_score → cells are dividing (less mature)
# If 4001NC has high stress markers → real QC issue

# %%

# Before any filtering, check what happens at different thresholds per sample
print("CURRENT FILTERING IMPACT BY SAMPLE:")
print(f"Total cells: {adata.n_obs:,}")
print(f"Total genes: {adata.n_vars:,}\n")

for sample in sorted(adata.obs['sample'].unique()):
    sample_adata = adata[adata.obs['sample'] == sample]
    print(f"\n{sample}: {sample_adata.n_obs:,} cells")
    print("-" * 50)
    
    # Test MIN_GENES thresholds
    for min_genes in [50, 100, 200, 500, 600, 700, 800, 1000, 1500]:
        cells_kept = (sample_adata.obs['n_genes_by_counts'] >= min_genes).sum()
        pct = 100 * cells_kept / sample_adata.n_obs
        print(f"  MIN_GENES={min_genes:4d}: {cells_kept:,} cells ({pct:5.1f}%)")
    
    print()
    
    # # Test MIN_CELLS_PER_GENE thresholds
    # for min_cells in [1, 3, 5, 10, 20, 50]:
    #     # Simulate filtering genes
    #     gene_counts = (sample_adata.X > 0).sum(axis=0).A1 if hasattr(sample_adata.X, 'A1') else (sample_adata.X > 0).sum(axis=0)
    #     genes_kept = (gene_counts >= min_cells).sum()
    #     pct = 100 * genes_kept / sample_adata.n_vars
    #     print(f"  MIN_CELLS_PER_GENE={min_cells:2d}: {genes_kept:,} genes ({pct:5.1f}%)")

# %%
import matplotlib.pyplot as plt
import numpy as np

# Prepare data
samples = sorted(adata.obs['sample'].unique())
min_genes_range = [50, 100, 200, 500, 600, 700, 800, 1000, 1500]

# Calculate retention for each sample
retention_data = {}
for sample in samples:
    sample_adata = adata[adata.obs['sample'] == sample]
    retention = []
    for min_g in min_genes_range:
        kept = (sample_adata.obs['n_genes_by_counts'] >= min_g).sum()
        pct = 100 * kept / sample_adata.n_obs
        retention.append(pct)
    retention_data[sample] = retention

# Plot 1: Retention curves by sample
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Cell Retention Analysis by MIN_GENES Threshold', fontsize=14, fontweight='bold')

# Full view
ax = axes[0, 0]
colors = plt.cm.Set2(np.linspace(0, 1, len(samples)))
for sample, color in zip(samples, colors):
    ax.plot(min_genes_range, retention_data[sample], marker='o', label=sample, linewidth=2, color=color)
ax.axvline(600, color='orange', linestyle='--', linewidth=2, alpha=0.7, label='MIN_GENES=600')
ax.axvline(800, color='red', linestyle='--', linewidth=2, alpha=0.7, label='MIN_GENES=800')
ax.set_xlabel('MIN_GENES Threshold')
ax.set_ylabel('Cells Retained (%)')
ax.set_title('Full View: All Samples')
ax.legend(loc='best', fontsize=9)
ax.grid(True, alpha=0.3)

# Zoomed: 500-1000 range
ax = axes[0, 1]
min_genes_zoom = [g for g in min_genes_range if 500 <= g <= 1000]
for sample, color in zip(samples, colors):
    idx = [i for i, g in enumerate(min_genes_range) if 500 <= g <= 1000]
    retention_zoom = [retention_data[sample][i] for i in idx]
    ax.plot(min_genes_zoom, retention_zoom, marker='o', label=sample, linewidth=2, color=color)
ax.axvline(600, color='orange', linestyle='--', linewidth=2, alpha=0.7, label='MIN_GENES=600 (recommended)')
ax.axvline(800, color='red', linestyle='--', linewidth=2, alpha=0.7, label='MIN_GENES=800 (strict)')
ax.axhline(90, color='gray', linestyle=':', alpha=0.5, label='90% retention threshold')
ax.set_xlabel('MIN_GENES Threshold')
ax.set_ylabel('Cells Retained (%)')
ax.set_title('Zoomed: 500-1000 Range (Decision Zone)')
ax.legend(loc='best', fontsize=9)
ax.grid(True, alpha=0.3)

# Heatmap: retention at key thresholds
ax = axes[1, 0]
key_thresholds = [500, 600, 800, 1000]
heatmap_data = []
for sample in samples:
    row = []
    for t in key_thresholds:
        idx = min_genes_range.index(t)
        row.append(retention_data[sample][idx])
    heatmap_data.append(row)

im = ax.imshow(heatmap_data, cmap='RdYlGn', aspect='auto', vmin=70, vmax=100)
ax.set_xticks(range(len(key_thresholds)))
ax.set_xticklabels([f'MIN={t}' for t in key_thresholds])
ax.set_yticks(range(len(samples)))
ax.set_yticklabels(samples)
ax.set_title('Cell Retention Heatmap (%)')
plt.colorbar(im, ax=ax, label='% Cells Retained')

# Add text annotations
for i in range(len(samples)):
    for j in range(len(key_thresholds)):
        text = ax.text(j, i, f'{heatmap_data[i][j]:.1f}%',
                      ha="center", va="center", color="black", fontsize=9, fontweight='bold')

# Summary statistics table
ax = axes[1, 1]
ax.axis('off')

# Create summary text
summary_text = "RECOMMENDATION:\n\n"
summary_text += "Option A (Adaptive - BEST):\n"
summary_text += "  4002G3: MIN_GENES=600\n"
summary_text += "  Others: MIN_GENES=800\n"
summary_text += "  Rationale: Balanced retention across all samples\n\n"

summary_text += "Option B (Moderate - Simpler):\n"
summary_text += "  All samples: MIN_GENES=600\n"
summary_text += "  Rationale: Uniform, minimal cell loss\n\n"

summary_text += "Avoid:\n"
summary_text += "  MIN_GENES=800 uniformly\n"
summary_text += f"  → 4002G3 loses {100-81:.0f}% of cells ⚠️\n\n"

summary_text += f"MIN_CELLS_PER_GENE = 5\n"
summary_text += "(All samples: ~76-77% genes retained)"

ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
        fontsize=11, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/source/Chats/MIN_GENES_threshold_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Plot saved!")

# %%
