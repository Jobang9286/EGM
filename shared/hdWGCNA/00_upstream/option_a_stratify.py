import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, mannwhitneyu
import os, gc, warnings
warnings.filterwarnings('ignore')

PIPE = "EGM_data/senepy_pipeline"
FIG_DIR = os.path.join(PIPE, "figures")
OUT_DIR = os.path.join(PIPE, "responder_analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load UMAP-integrated data + SCAP scores ──
print("Loading UMAP-integrated data...", flush=True)
adata = sc.read_h5ad(os.path.join(PIPE, "EGM_umap_integrated.h5ad"))
print(f"  {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

print("Loading SCAP scores...", flush=True)
scored = sc.read_h5ad(os.path.join(PIPE, "step4_scores/EGM_all_senepy_scored.h5ad"), backed='r')
for col in scored.obs.columns:
    if col.startswith('senepy_hub'):
        adata.obs[col] = scored.obs[col].values
del scored
gc.collect()

# ── Step 1: R3 SCAP-0 distribution — bimodal check ──
print("\n" + "="*60, flush=True)
print("Step 1: R3 SCAP-0 score distribution", flush=True)
print("="*60, flush=True)

r3_mask = adata.obs['sample'].values == 'R3'
r3_scores = adata.obs.loc[r3_mask, 'senepy_hub0'].values
y1_scores = adata.obs.loc[adata.obs['sample']=='Y1', 'senepy_hub0'].values
a1_scores = adata.obs.loc[adata.obs['sample']=='A1', 'senepy_hub0'].values

print(f"  R3 cells: {r3_mask.sum()}", flush=True)
print(f"  R3 SCAP-0: mean={r3_scores.mean():.4f}, median={np.median(r3_scores):.4f}", flush=True)
print(f"  R3 SCAP-0: Q25={np.percentile(r3_scores,25):.4f}, Q75={np.percentile(r3_scores,75):.4f}", flush=True)

# Young threshold for "rejuvenated": Q75 of Young
y1_q75 = np.percentile(y1_scores, 75)
y1_q90 = np.percentile(y1_scores, 90)
# Aged threshold: Q25 of Aged
a1_q25 = np.percentile(a1_scores, 25)

print(f"\n  Young Q75: {y1_q75:.4f}, Young Q90: {y1_q90:.4f}", flush=True)
print(f"  Aged Q25: {a1_q25:.4f}", flush=True)

# Define responder/non-responder using Young Q75 as cutoff
# Responder: SCAP-0 <= Young Q75 (fully rejuvenated)
# Intermediate: Young Q75 < SCAP-0 <= Aged Q25
# Non-responder: SCAP-0 > Aged Q25 (still senescent-like)
resp_mask_r3 = r3_scores <= y1_q75
inter_mask_r3 = (r3_scores > y1_q75) & (r3_scores <= a1_q25)
nonresp_mask_r3 = r3_scores > a1_q25

n_resp = resp_mask_r3.sum()
n_inter = inter_mask_r3.sum()
n_nonresp = nonresp_mask_r3.sum()
total = len(r3_scores)

print(f"\n  R3 Stratification (SCAP-0 threshold):", flush=True)
print(f"    Responder (≤ Young Q75={y1_q75:.4f}):      {n_resp:>6} ({100*n_resp/total:.1f}%)", flush=True)
print(f"    Intermediate:                               {n_inter:>6} ({100*n_inter/total:.1f}%)", flush=True)
print(f"    Non-responder (> Aged Q25={a1_q25:.4f}):    {n_nonresp:>6} ({100*n_nonresp/total:.1f}%)", flush=True)

# ── Step 2: Apply same classification to all timepoints ──
print("\n" + "="*60, flush=True)
print("Step 2: Responder fraction across all timepoints", flush=True)
print("="*60, flush=True)

sample_order = ['Y1','A1','R0','R1','R2','R3','R4','R5','R6','R7','R8','R9']
frac_rows = []

for sample in sample_order:
    s_mask = adata.obs['sample'].values == sample
    s_scores = adata.obs.loc[s_mask, 'senepy_hub0'].values
    n = len(s_scores)
    n_r = (s_scores <= y1_q75).sum()
    n_i = ((s_scores > y1_q75) & (s_scores <= a1_q25)).sum()
    n_nr = (s_scores > a1_q25).sum()
    
    frac_rows.append({
        'sample': sample, 'n_cells': n,
        'responder_n': n_r, 'responder_pct': 100*n_r/n,
        'intermediate_n': n_i, 'intermediate_pct': 100*n_i/n,
        'nonresponder_n': n_nr, 'nonresponder_pct': 100*n_nr/n,
    })
    print(f"  {sample:>3}: Resp={100*n_r/n:5.1f}%  Inter={100*n_i/n:5.1f}%  NonResp={100*n_nr/n:5.1f}%", flush=True)

frac_df = pd.DataFrame(frac_rows)
frac_df.to_csv(os.path.join(OUT_DIR, 'responder_fractions.csv'), index=False)

# ── Step 3: Annotate cells ──
print("\n" + "="*60, flush=True)
print("Step 3: Annotating all cells with responder status", flush=True)
print("="*60, flush=True)

scap0 = adata.obs['senepy_hub0'].values
status = np.full(len(scap0), 'intermediate', dtype=object)
status[scap0 <= y1_q75] = 'responder'
status[scap0 > a1_q25] = 'non-responder'
adata.obs['scap0_status'] = status

# ── Step 4: DEG — R3 responder vs non-responder ──
print("\n" + "="*60, flush=True)
print("Step 4: DEG analysis — R3 responders vs non-responders", flush=True)
print("="*60, flush=True)

r3_adata = adata[r3_mask].copy()

# Need normalized data for DEG
if 'counts' in r3_adata.layers:
    r3_adata.X = r3_adata.layers['counts'].copy()
sc.pp.normalize_total(r3_adata, target_sum=1e4)
sc.pp.log1p(r3_adata)

# Only compare responder vs non-responder (drop intermediate for cleaner signal)
r3_rn = r3_adata[r3_adata.obs['scap0_status'].isin(['responder', 'non-responder'])].copy()
print(f"  R3 responder: {(r3_rn.obs['scap0_status']=='responder').sum()}", flush=True)
print(f"  R3 non-responder: {(r3_rn.obs['scap0_status']=='non-responder').sum()}", flush=True)

sc.tl.rank_genes_groups(r3_rn, groupby='scap0_status', groups=['non-responder'],
                         reference='responder', method='wilcoxon', pts=True)

# Extract results
deg_df = sc.get.rank_genes_groups_df(r3_rn, group='non-responder')
deg_df = deg_df.sort_values('scores', ascending=False)
deg_df.to_csv(os.path.join(OUT_DIR, 'R3_nonresponder_vs_responder_DEG.csv'), index=False)

# Top upregulated in non-responder (gate-keeper candidates — genes that keep cells senescent)
up_in_nonresp = deg_df[(deg_df['logfoldchanges'] > 0.5) & (deg_df['pvals_adj'] < 0.05)]
down_in_nonresp = deg_df[(deg_df['logfoldchanges'] < -0.5) & (deg_df['pvals_adj'] < 0.05)]

print(f"\n  DEG results (adj.p < 0.05, |logFC| > 0.5):", flush=True)
print(f"    Up in non-responder:   {len(up_in_nonresp)} genes", flush=True)
print(f"    Down in non-responder: {len(down_in_nonresp)} genes", flush=True)

print(f"\n  Top 20 UP in non-responder (gate-keeper candidates):", flush=True)
for _, row in up_in_nonresp.head(20).iterrows():
    print(f"    {row['names']:15s}  logFC={row['logfoldchanges']:>6.2f}  padj={row['pvals_adj']:.2e}", flush=True)

print(f"\n  Top 20 DOWN in non-responder (rejuvenation markers):", flush=True)
for _, row in down_in_nonresp.head(20).iterrows():
    print(f"    {row['names']:15s}  logFC={row['logfoldchanges']:>6.2f}  padj={row['pvals_adj']:.2e}", flush=True)

# ── Step 5: Save annotated data ──
adata.write_h5ad(os.path.join(OUT_DIR, 'EGM_annotated_responder.h5ad'))
print(f"\nSaved annotated data: {OUT_DIR}/EGM_annotated_responder.h5ad", flush=True)

# ── Step 6: Plots ──
print("\nGenerating plots...", flush=True)

# Plot A: R3 SCAP-0 density with thresholds
fig, ax = plt.subplots(figsize=(5, 3), dpi=300)
x_grid = np.linspace(-0.02, 0.3, 500)

for scores, color, label in [(y1_scores, '#27AE60', 'Young'),
                               (r3_scores, '#F5B041', 'R3 (EGM D3)'),
                               (a1_scores, '#E8524A', 'Aged')]:
    kde = gaussian_kde(scores, bw_method=0.15)
    density = kde(x_grid)
    ax.fill_between(x_grid, density, alpha=0.3, color=color)
    ax.plot(x_grid, density, color=color, linewidth=1.5, label=label)

ax.axvline(y1_q75, color='#27AE60', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axvline(a1_q25, color='#E8524A', linestyle='--', linewidth=0.8, alpha=0.6)
ax.text(y1_q75, ax.get_ylim()[1]*0.9, 'Young\nQ75', fontsize=7, color='#27AE60', ha='center')
ax.text(a1_q25, ax.get_ylim()[1]*0.9, 'Aged\nQ25', fontsize=7, color='#E8524A', ha='center')

ax.set_xlabel('SCAP-0 Score', fontsize=12, fontweight='bold')
ax.set_ylabel('Density', fontsize=12, fontweight='bold')
ax.tick_params(labelsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=8, loc='upper right')
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'r3_scap0_density_thresholds.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

# Plot B: Stacked bar — responder fraction across timepoints
fig, ax = plt.subplots(figsize=(7, 3.5), dpi=300)
x_pos = np.arange(len(sample_order))
resp_pcts = frac_df['responder_pct'].values
inter_pcts = frac_df['intermediate_pct'].values
nonresp_pcts = frac_df['nonresponder_pct'].values

ax.bar(x_pos, resp_pcts, color='#27AE60', label='Responder', edgecolor='black', linewidth=0.3)
ax.bar(x_pos, inter_pcts, bottom=resp_pcts, color='#F5B041', label='Intermediate', edgecolor='black', linewidth=0.3)
ax.bar(x_pos, nonresp_pcts, bottom=resp_pcts+inter_pcts, color='#E8524A', label='Non-responder', edgecolor='black', linewidth=0.3)

ax.axvline(x=5, color='#D62828', linestyle='--', linewidth=1.0, alpha=0.5)  # R3 position in sample_order

ax.set_xticks(x_pos)
ax.set_xticklabels(sample_order, fontsize=10, fontweight='bold')
ax.set_ylabel('Cell Fraction (%)', fontsize=12, fontweight='bold')
ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
ax.tick_params(axis='y', labelsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=8, loc='upper right')
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'responder_fraction_stacked_bar.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

# Plot C: UMAP colored by responder status (R3 only)
fig, ax = plt.subplots(figsize=(5.5, 4.5), dpi=300)
umap = adata.obsm['X_umap']
# Background: all cells gray
ax.scatter(umap[:, 0], umap[:, 1], c='#f0f0f0', s=0.1, alpha=0.3, rasterized=True)
# R3 cells colored by status
r3_idx = np.where(r3_mask)[0]
status_colors = {'responder': '#27AE60', 'intermediate': '#F5B041', 'non-responder': '#E8524A'}
for st, col in status_colors.items():
    st_mask = adata.obs.iloc[r3_idx]['scap0_status'].values == st
    ax.scatter(umap[r3_idx[st_mask], 0], umap[r3_idx[st_mask], 1],
               c=col, s=0.5, alpha=0.5, rasterized=True, label=st)

from matplotlib.patches import Patch
handles = [Patch(facecolor=c, label=l) for l, c in status_colors.items()]
ax.legend(handles=handles, fontsize=8, loc='upper right', markerscale=2)
ax.set_xlabel('UMAP 1', fontsize=12, fontweight='bold')
ax.set_ylabel('UMAP 2', fontsize=12, fontweight='bold')
ax.tick_params(labelsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'umap_r3_responder_status.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

print("\nAll plots saved.", flush=True)
print("\n" + "="*60, flush=True)
print("OPTION A COMPLETE", flush=True)
print("="*60, flush=True)
