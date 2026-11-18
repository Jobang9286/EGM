#!/usr/bin/env python
"""
EGM scRNA-seq Preprocessing Pipeline
=====================================
- Load 10x mtx per sample (Y1, A1, R0-R9)
- Doublet detection (Scrublet)
- QC filtering (UMI, genes, mito%)
- Save individual h5ad per sample + merged raw h5ad
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ── Configuration ──
BASE_DIR = "."
OUT_DIR  = os.path.join(BASE_DIR, "preprocessed")
FIG_DIR  = os.path.join(OUT_DIR, "qc_figures")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLES = ["Y1", "A1", "R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9"]

# Condition mapping
CONDITION_MAP = {
    "Y1": "Young",
    "A1": "Aged",
    "R0": "EGM_D0", "R1": "EGM_D1", "R2": "EGM_D2", "R3": "EGM_D3",
    "R4": "EGM_D4", "R5": "EGM_D5", "R6": "EGM_D6", "R7": "EGM_D7",
    "R8": "EGM_D8", "R9": "EGM_D9"
}

# Timepoint ordering (for UNAGI staging later)
TIMEPOINT_MAP = {
    "Y1": -1,   # Young control (reference)
    "A1": 0,    # Aged baseline
    "R0": 1, "R1": 2, "R2": 3, "R3": 4, "R4": 5,
    "R5": 6, "R6": 7, "R7": 8, "R8": 9, "R9": 10
}

# QC thresholds
MIN_GENES = 200
MAX_GENES = 8000
MIN_COUNTS = 500
MAX_MITO_PCT = 25.0      # mouse fibroblast - slightly relaxed
SCRUBLET_THRESHOLD = 0.25 # default; auto-threshold preferred

sc.settings.verbosity = 2
sc.settings.n_jobs = 48

# ── Helper Functions ──
def load_sample(sample_name):
    """Load 10x mtx and add metadata."""
    path = os.path.join(BASE_DIR, sample_name)
    print(f"\n{'='*60}")
    print(f"Loading {sample_name} from {path}")
    adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=False)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    
    # Add metadata
    adata.obs['sample'] = sample_name
    adata.obs['condition'] = CONDITION_MAP[sample_name]
    adata.obs['timepoint'] = TIMEPOINT_MAP[sample_name]
    
    # Categorize condition type
    if sample_name == "Y1":
        adata.obs['group'] = "Young"
    elif sample_name == "A1":
        adata.obs['group'] = "Aged"
    else:
        adata.obs['group'] = "EGM_treated"
    
    print(f"  Raw: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def run_doublet_detection(adata, sample_name):
    """Run Scrublet for doublet detection."""
    print(f"  Running Scrublet on {sample_name}...")
    
    # Use raw counts for scrublet
    counts_matrix = adata.X.copy()
    
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2, min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=30
    )
    
    adata.obs['doublet_score'] = doublet_scores
    
    # Use auto-threshold if available, else fallback
    if scrub.threshold_ is not None and not np.isnan(scrub.threshold_):
        threshold = scrub.threshold_
        adata.obs['predicted_doublet'] = predicted_doublets
    else:
        threshold = SCRUBLET_THRESHOLD
        adata.obs['predicted_doublet'] = doublet_scores > threshold
    
    n_doublets = adata.obs['predicted_doublet'].sum()
    pct_doublets = 100 * n_doublets / adata.n_obs
    print(f"  Scrublet threshold: {threshold:.3f}")
    print(f"  Doublets detected: {n_doublets} ({pct_doublets:.1f}%)")
    
    # Save scrublet histogram
    try:
        fig, axes = scrub.plot_histogram()
        fig.suptitle(f'{sample_name} - Scrublet Scores', fontsize=14)
        fig.savefig(os.path.join(FIG_DIR, f'{sample_name}_scrublet_hist.png'), 
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
    except:
        pass
    
    return adata, threshold


def qc_filter(adata, sample_name):
    """Calculate QC metrics and filter cells."""
    print(f"  QC filtering {sample_name}...")
    
    # Mouse mitochondrial genes start with "mt-"
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    # Also check for "MT-" just in case
    adata.var['mt'] = adata.var['mt'] | adata.var_names.str.startswith('MT-')
    
    # Ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl', 'RPS', 'RPL'))
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], 
                                percent_top=None, log1p=False, inplace=True)
    
    n_before = adata.n_obs
    
    # Print QC stats before filtering
    print(f"    Pre-filter stats:")
    print(f"      Median genes/cell: {np.median(adata.obs['n_genes_by_counts']):.0f}")
    print(f"      Median UMI/cell:   {np.median(adata.obs['total_counts']):.0f}")
    print(f"      Median mito%:      {np.median(adata.obs['pct_counts_mt']):.1f}%")
    
    # Remove doublets first
    adata = adata[~adata.obs['predicted_doublet']].copy()
    n_after_doublet = adata.n_obs
    
    # Apply QC filters
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
    adata = adata[adata.obs['n_genes_by_counts'] < MAX_GENES].copy()
    adata = adata[adata.obs['pct_counts_mt'] < MAX_MITO_PCT].copy()
    
    # Filter genes (at least in 3 cells)
    sc.pp.filter_genes(adata, min_cells=3)
    
    n_after = adata.n_obs
    print(f"    Removed: {n_before - n_after_doublet} doublets, "
          f"{n_after_doublet - n_after} low-quality cells")
    print(f"    Retained: {n_after} cells ({100*n_after/n_before:.1f}%)")
    
    return adata


def plot_qc(adata, sample_name):
    """Generate QC violin plots."""
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    
    for ax, metric, title in zip(axes, 
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'doublet_score'],
        ['Genes per Cell', 'UMI Counts', 'Mito %', 'Doublet Score']):
        
        if metric in adata.obs.columns:
            ax.hist(adata.obs[metric], bins=50, color='steelblue', edgecolor='none', alpha=0.7)
            ax.set_title(f'{sample_name}: {title}')
            ax.set_xlabel(title)
            ax.set_ylabel('Cells')
    
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, f'{sample_name}_qc_violin.png'), 
                dpi=150, bbox_inches='tight')
    plt.close(fig)


# ── Main Pipeline ──
def main():
    print("=" * 60)
    print("EGM scRNA-seq Preprocessing Pipeline")
    print("=" * 60)
    
    summary_records = []
    all_adatas = []
    
    for sample in SAMPLES:
        # 1. Load
        adata = load_sample(sample)
        n_raw = adata.n_obs
        
        # 2. Doublet detection
        adata, scrub_thresh = run_doublet_detection(adata, sample)
        
        # 3. QC plots (pre-filter, but with scores)
        # Calculate QC metrics for plotting before filtering
        adata.var['mt'] = adata.var_names.str.startswith('mt-') | adata.var_names.str.startswith('MT-')
        adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl', 'RPS', 'RPL'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'],
                                    percent_top=None, log1p=False, inplace=True)
        plot_qc(adata, sample)
        
        # 4. QC filtering
        adata = qc_filter(adata, sample)
        n_filtered = adata.n_obs
        
        # 5. Save individual h5ad (raw counts preserved)
        out_path = os.path.join(OUT_DIR, f"{sample}_filtered.h5ad")
        adata.write_h5ad(out_path)
        print(f"  Saved: {out_path}")
        
        # Record summary
        summary_records.append({
            'sample': sample,
            'condition': CONDITION_MAP[sample],
            'group': adata.obs['group'].iloc[0],
            'timepoint': TIMEPOINT_MAP[sample],
            'n_raw': n_raw,
            'n_doublets': n_raw - adata.n_obs + (n_raw - n_filtered - (n_raw - adata.n_obs)),
            'n_filtered': n_filtered,
            'pct_retained': round(100 * n_filtered / n_raw, 1),
            'scrublet_threshold': round(scrub_thresh, 3),
            'median_genes': int(np.median(adata.obs['n_genes_by_counts'])),
            'median_umi': int(np.median(adata.obs['total_counts'])),
            'median_mito_pct': round(np.median(adata.obs['pct_counts_mt']), 2)
        })
        
        all_adatas.append(adata)
    
    # ── Merge all samples ──
    print(f"\n{'='*60}")
    print("Merging all samples...")
    
    merged = sc.concat(all_adatas, join='outer', label='sample', 
                       keys=SAMPLES, index_unique='-')
    
    # Ensure metadata is preserved
    print(f"Merged: {merged.n_obs} cells x {merged.n_vars} genes")
    print(f"Samples: {merged.obs['sample'].value_counts().to_dict()}")
    
    # Save merged
    merged_path = os.path.join(OUT_DIR, "EGM_all_samples_filtered.h5ad")
    merged.write_h5ad(merged_path)
    print(f"Saved merged: {merged_path}")
    
    # ── Summary table ──
    summary_df = pd.DataFrame(summary_records)
    summary_path = os.path.join(OUT_DIR, "preprocessing_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSaved summary: {summary_path}")
    print("\n" + summary_df.to_string(index=False))
    
    # ── Summary figure ──
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    x = range(len(summary_df))
    labels = summary_df['sample'].tolist()
    colors = ['#2ecc71' if s == 'Y1' else '#e74c3c' if s == 'A1' else '#3498db' 
              for s in labels]
    
    axes[0].bar(x, summary_df['n_raw'], color=colors, alpha=0.4, label='Raw')
    axes[0].bar(x, summary_df['n_filtered'], color=colors, alpha=0.9, label='Filtered')
    axes[0].set_xticks(x); axes[0].set_xticklabels(labels, rotation=45)
    axes[0].set_ylabel('Number of Cells'); axes[0].set_title('Cell Counts (Raw vs Filtered)')
    axes[0].legend()
    
    axes[1].bar(x, summary_df['median_genes'], color=colors)
    axes[1].set_xticks(x); axes[1].set_xticklabels(labels, rotation=45)
    axes[1].set_ylabel('Median Genes'); axes[1].set_title('Median Genes per Cell')
    
    axes[2].bar(x, summary_df['median_mito_pct'], color=colors)
    axes[2].set_xticks(x); axes[2].set_xticklabels(labels, rotation=45)
    axes[2].set_ylabel('Median Mito %'); axes[2].set_title('Median Mitochondrial %')
    axes[2].axhline(y=MAX_MITO_PCT, color='red', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'preprocessing_summary.png'), 
                dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print(f"Total cells retained: {merged.n_obs}")
    print(f"Output directory: {OUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
