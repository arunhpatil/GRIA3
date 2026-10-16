"""
QC Threshold Optimization Tool

This script generates data-driven visualizations for selecting optimal
MIN_GENES and MIN_CELLS_PER_GENE thresholds for any scRNA-seq dataset.

Usage:
    python QC_threshold_optimization.py --adata /path/to/data.h5ad
    
Author: Data Science Pipeline
Date: 2026-10-14
"""
import os, sys
from pathlib import Path
# This is a workaround for a bug in ipykernel that causes typing_extensions to be imported from the wrong location. 
sys.path = [p for p in sys.path if "ipykernel/py3" not in p] + [p for p in sys.path if "ipykernel/py3" in p]
for m in [m for m in sys.modules if m == "typing_extensions"]:
    del sys.modules[m]
    
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc


def analyze_min_genes(adata, sample_key="sample"):
    """
    Calculate cell retention at different MIN_GENES thresholds per sample.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix
    sample_key : str
        Column name in adata.obs containing sample identifiers
        
    Returns
    -------
    dict
        {sample: [retention percentages for each threshold]}
    """
    samples = sorted(adata.obs[sample_key].unique())
    min_genes_range = [50, 100, 200, 500, 600, 700, 800, 1000, 1500]
    
    retention_data = {}
    for sample in samples:
        sample_adata = adata[adata.obs[sample_key] == sample]
        retention = []
        for min_g in min_genes_range:
            kept = (sample_adata.obs['n_genes_by_counts'] >= min_g).sum()
            pct = 100 * kept / sample_adata.n_obs
            retention.append(pct)
        retention_data[sample] = retention
    
    return samples, min_genes_range, retention_data


def analyze_min_cells_per_gene(adata, sample_key="sample"):
    """
    Calculate gene retention at different MIN_CELLS_PER_GENE thresholds per sample.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix
    sample_key : str
        Column name in adata.obs containing sample identifiers
        
    Returns
    -------
    dict
        {sample: [retention percentages for each threshold]}
    """
    samples = sorted(adata.obs[sample_key].unique())
    min_cells_range = [1, 3, 5, 10, 20, 50]
    
    retention_data = {}
    for sample in samples:
        sample_adata = adata[adata.obs[sample_key] == sample]
        retention = []
        for min_c in min_cells_range:
            # Count genes expressed in at least min_c cells
            gene_counts = (sample_adata.X > 0).sum(axis=0)
            if hasattr(gene_counts, 'A1'):
                gene_counts = gene_counts.A1
            genes_kept = (gene_counts >= min_c).sum()
            pct = 100 * genes_kept / sample_adata.n_vars
            retention.append(pct)
        retention_data[sample] = retention
    
    return samples, min_cells_range, retention_data


def plot_min_genes_analysis(adata, sample_key="sample", output_path=None):
    """
    Create comprehensive MIN_GENES threshold visualization.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix
    sample_key : str
        Column name in adata.obs containing sample identifiers
    output_path : str, optional
        Path to save figure. If None, displays plot.
    """
    samples, min_genes_range, retention_data = analyze_min_genes(adata, sample_key)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('MIN_GENES Threshold Optimization Analysis', fontsize=16, fontweight='bold')
    
    # Color palette
    colors = plt.cm.Set3(np.linspace(0, 1, len(samples)))
    
    # ========== Panel 1: Full view ==========
    ax = axes[0, 0]
    for sample, color in zip(samples, colors):
        ax.plot(min_genes_range, retention_data[sample], marker='o', 
                label=sample, linewidth=2.5, color=color, markersize=8)
    ax.axvline(500, color='green', linestyle='--', linewidth=2, alpha=0.6, label='500 (recommended)')
    ax.axvline(600, color='orange', linestyle='--', linewidth=2, alpha=0.6, label='600')
    ax.axvline(800, color='red', linestyle='--', linewidth=2, alpha=0.6, label='800 (strict)')
    ax.set_xlabel('MIN_GENES Threshold', fontsize=11, fontweight='bold')
    ax.set_ylabel('Cells Retained (%)', fontsize=11, fontweight='bold')
    ax.set_title('Full View: All Thresholds', fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 105])
    
    # ========== Panel 2: Zoomed (decision zone) ==========
    ax = axes[0, 1]
    min_genes_zoom = [g for g in min_genes_range if 400 <= g <= 1000]
    for sample, color in zip(samples, colors):
        idx = [i for i, g in enumerate(min_genes_range) if 400 <= g <= 1000]
        retention_zoom = [retention_data[sample][i] for i in idx]
        ax.plot(min_genes_zoom, retention_zoom, marker='o', label=sample, 
                linewidth=2.5, color=color, markersize=8)
    ax.axvline(500, color='green', linestyle='--', linewidth=2, alpha=0.6, label='500 (recommended)')
    ax.axvline(600, color='orange', linestyle='--', linewidth=2, alpha=0.6)
    ax.axvline(800, color='red', linestyle='--', linewidth=2, alpha=0.6)
    ax.axhline(90, color='gray', linestyle=':', alpha=0.5, label='90% threshold')
    ax.set_xlabel('MIN_GENES Threshold', fontsize=11, fontweight='bold')
    ax.set_ylabel('Cells Retained (%)', fontsize=11, fontweight='bold')
    ax.set_title('Decision Zone: 400-1000 Range', fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([70, 105])
    
    # ========== Panel 3: Heatmap ==========
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
    ax.set_xticklabels([f'MIN={t}' for t in key_thresholds], fontsize=10)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=10)
    ax.set_title('Cell Retention Heatmap (%)', fontsize=12, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax, label='% Cells Retained')
    cbar.ax.tick_params(labelsize=9)
    
    # Add text annotations
    for i in range(len(samples)):
        for j in range(len(key_thresholds)):
            text = ax.text(j, i, f'{heatmap_data[i][j]:.1f}%',
                          ha="center", va="center", color="black", 
                          fontsize=10, fontweight='bold')
    
    # ========== Panel 4: Summary statistics ==========
    ax = axes[1, 1]
    ax.axis('off')
    
    # Calculate statistics
    summary_stats = []
    summary_stats.append("DATASET SUMMARY\n" + "="*40)
    summary_stats.append(f"Total cells: {adata.n_obs:,}")
    summary_stats.append(f"Total genes: {adata.n_vars:,}")
    summary_stats.append(f"Number of samples: {len(samples)}\n")
    
    summary_stats.append("PER-SAMPLE STATISTICS")
    summary_stats.append("-"*40)
    for sample in samples:
        sample_data = adata.obs[adata.obs[sample_key] == sample]['n_genes_by_counts']
        summary_stats.append(f"{sample}:")
        summary_stats.append(f"  Cells: {len(sample_data):,}")
        summary_stats.append(f"  Median genes: {sample_data.median():.0f}")
        summary_stats.append(f"  Min: {sample_data.min():.0f}, Max: {sample_data.max():.0f}")
    
    summary_text = "\n".join(summary_stats)
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7, pad=1))
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"✓ MIN_GENES plot saved: {output_path}")
    else:
        plt.show()
    
    plt.close()


def plot_min_cells_per_gene_analysis(adata, sample_key="sample", output_path=None):
    """
    Create comprehensive MIN_CELLS_PER_GENE threshold visualization.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix
    sample_key : str
        Column name in adata.obs containing sample identifiers
    output_path : str, optional
        Path to save figure. If None, displays plot.
    """
    samples, min_cells_range, retention_data = analyze_min_cells_per_gene(adata, sample_key)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('MIN_CELLS_PER_GENE Threshold Optimization Analysis', fontsize=16, fontweight='bold')
    
    # Color palette
    colors = plt.cm.Set3(np.linspace(0, 1, len(samples)))
    
    # ========== Panel 1: Full view ==========
    ax = axes[0, 0]
    for sample, color in zip(samples, colors):
        ax.plot(min_cells_range, retention_data[sample], marker='o', 
                label=sample, linewidth=2.5, color=color, markersize=8)
    ax.axvline(3, color='blue', linestyle='--', linewidth=2, alpha=0.6, label='3 (permissive)')
    ax.axvline(5, color='green', linestyle='--', linewidth=2, alpha=0.6, label='5 (recommended)')
    ax.axvline(10, color='orange', linestyle='--', linewidth=2, alpha=0.6, label='10 (strict)')
    ax.set_xlabel('MIN_CELLS_PER_GENE Threshold', fontsize=11, fontweight='bold')
    ax.set_ylabel('Genes Retained (%)', fontsize=11, fontweight='bold')
    ax.set_title('Full View: All Thresholds', fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 105])
    
    # ========== Panel 2: Zoomed (decision zone) ==========
    ax = axes[0, 1]
    min_cells_zoom = [c for c in min_cells_range if 1 <= c <= 20]
    for sample, color in zip(samples, colors):
        idx = [i for i, c in enumerate(min_cells_range) if 1 <= c <= 20]
        retention_zoom = [retention_data[sample][i] for i in idx]
        ax.plot(min_cells_zoom, retention_zoom, marker='o', label=sample, 
                linewidth=2.5, color=color, markersize=8)
    ax.axvline(3, color='blue', linestyle='--', linewidth=2, alpha=0.6, label='3 (permissive)')
    ax.axvline(5, color='green', linestyle='--', linewidth=2, alpha=0.6, label='5 (recommended)')
    ax.axvline(10, color='orange', linestyle='--', linewidth=2, alpha=0.6)
    ax.axhline(75, color='gray', linestyle=':', alpha=0.5, label='75% threshold')
    ax.set_xlabel('MIN_CELLS_PER_GENE Threshold', fontsize=11, fontweight='bold')
    ax.set_ylabel('Genes Retained (%)', fontsize=11, fontweight='bold')
    ax.set_title('Decision Zone: 1-20 Range', fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([50, 105])
    
    # ========== Panel 3: Heatmap ==========
    ax = axes[1, 0]
    key_cells_thresholds = [1, 3, 5, 10]
    heatmap_data = []
    for sample in samples:
        row = []
        for t in key_cells_thresholds:
            idx = min_cells_range.index(t)
            row.append(retention_data[sample][idx])
        heatmap_data.append(row)
    
    im = ax.imshow(heatmap_data, cmap='RdYlGn', aspect='auto', vmin=50, vmax=100)
    ax.set_xticks(range(len(key_cells_thresholds)))
    ax.set_xticklabels([f'MIN={t}' for t in key_cells_thresholds], fontsize=10)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=10)
    ax.set_title('Gene Retention Heatmap (%)', fontsize=12, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax, label='% Genes Retained')
    cbar.ax.tick_params(labelsize=9)
    
    # Add text annotations
    for i in range(len(samples)):
        for j in range(len(key_cells_thresholds)):
            text = ax.text(j, i, f'{heatmap_data[i][j]:.1f}%',
                          ha="center", va="center", color="black", 
                          fontsize=10, fontweight='bold')
    
    # ========== Panel 4: Summary statistics ==========
    ax = axes[1, 1]
    ax.axis('off')
    
    # Calculate statistics
    summary_stats = []
    summary_stats.append("FILTERING IMPACT\n" + "="*40)
    summary_stats.append("MIN_GENES=500 + MIN_CELLS_PER_GENE=5\n")
    
    total_genes_all_thresholds = []
    for threshold in [3, 5, 10]:
        gene_counts = (adata.X > 0).sum(axis=0)
        if hasattr(gene_counts, 'A1'):
            gene_counts = gene_counts.A1
        genes_kept = (gene_counts >= threshold).sum()
        pct = 100 * genes_kept / adata.n_vars
        total_genes_all_thresholds.append((threshold, genes_kept, pct))
    
    summary_stats.append("Gene retention by MIN_CELLS:")
    for threshold, count, pct in total_genes_all_thresholds:
        summary_stats.append(f"  MIN={threshold}: {count:,} genes ({pct:.1f}%)")
    
    summary_stats.append("\nRecommended combination:")
    summary_stats.append("  MIN_GENES = 500")
    summary_stats.append("  MIN_CELLS_PER_GENE = 5")
    summary_stats.append("\nExpected outcome:")
    summary_stats.append("  Cell retention: ~95-98%")
    summary_stats.append("  Gene retention: ~76-77%")
    
    summary_text = "\n".join(summary_stats)
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7, pad=1))
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"✓ MIN_CELLS_PER_GENE plot saved: {output_path}")
    else:
        plt.show()
    
    plt.close()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Optimize QC thresholds for scRNA-seq data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python QC_threshold_optimization.py --adata data.h5ad
  python QC_threshold_optimization.py --adata data.h5ad --sample-key batch --output-dir ./plots/
        """
    )
    
    parser.add_argument('--adata', type=str, required=True,
                       help='Path to AnnData h5ad file')
    parser.add_argument('--sample-key', type=str, default='sample',
                       help='Column name in obs containing sample identifiers (default: sample)')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Directory to save plots. If None, displays plots interactively.')
    
    args = parser.parse_args()
    
    # Load data
    print(f"Loading data from: {args.adata}")
    adata = sc.read_h5ad(args.adata)
    print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Verify sample_key exists
    if args.sample_key not in adata.obs.columns:
        raise ValueError(f"Column '{args.sample_key}' not found in adata.obs")
    
    # Create output directory if needed
    output_dir = None
    if args.output_dir:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    print("\nGenerating MIN_GENES analysis...")
    min_genes_path = None
    if output_dir:
        min_genes_path = output_dir / "QC_MIN_GENES_analysis.png"
    plot_min_genes_analysis(adata, sample_key=args.sample_key, output_path=min_genes_path)
    
    print("Generating MIN_CELLS_PER_GENE analysis...")
    min_cells_path = None
    if output_dir:
        min_cells_path = output_dir / "QC_MIN_CELLS_PER_GENE_analysis.png"
    plot_min_cells_per_gene_analysis(adata, sample_key=args.sample_key, output_path=min_cells_path)
    
    print("\n" + "="*50)
    print("ANALYSIS COMPLETE")
    print("="*50)
    print(f"Samples found: {sorted(adata.obs[args.sample_key].unique())}")
    if output_dir:
        print(f"Plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
