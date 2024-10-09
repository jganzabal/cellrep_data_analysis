# Filter out cells with high mitochondrial gene reads, and cells with sgRNAs
# with low representation

# %% 
# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# %% Import AnnData Object ---------------------------------------------------

tcells = sc.read_h5ad("data/tcells_all_unfilt.h5ad")
tcells
# %%
# Calculate RNA feature metrics and filter --------------------------------

# Calculate mitochondrial percentage
tcells.var['mt'] = tcells.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(tcells, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
tcells.var
# %%
# %%
# Calculate ribosomal percentage
tcells.var['ribo'] = tcells.var_names.str.match('^RP[SL]')
sc.pp.calculate_qc_metrics(tcells, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
# %%
# Plot statistics
sc.pl.violin(
    tcells,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
    jitter=0.4, groupby='condition', multi_panel=True, 
)
# %%
# Filter out dead cells (high MT) and multiplets (high N_features)
tcells = tcells[
    (tcells.obs['n_genes_by_counts'] > 400) &
    (tcells.obs['n_genes_by_counts'] < 6000) &
    (tcells.obs['pct_counts_mt'] < 25)
]
# %%
# Plot statistics after filtering
sc.pl.violin(tcells, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
             jitter=0.4, groupby='condition', multi_panel=True, 
)
# %%
# Plot CRISPR cell calls and filter ---------------------------------------

# Plot n cells per guide target 
guide_target_summary = tcells.obs.groupby(['condition', 'gene']).size().reset_index(name='n_cells')
print(guide_target_summary)

plt.figure(figsize=(12, 6))
sns.barplot(data=guide_target_summary.sort_values('n_cells'), 
            x='gene', y='n_cells', hue='condition')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
# %%
# Remove cells with guides targeting genes w < 100 cells in either group:
## IRX4, PRDM1, TCF7, and HELZ2

genes2filter = guide_target_summary[guide_target_summary['n_cells'] < 100]['gene'].unique()
print(genes2filter)

genes2keep = guide_target_summary[~guide_target_summary['gene'].isin(genes2filter)]['gene'].unique()
print(genes2keep)

tcells = tcells[tcells.obs['gene'].isin(genes2keep)]
print(tcells)
# %%
# Save AnnData Object -----------------------------------------------------

tcells.write_h5ad("data/tcells_all_filt.h5ad")
