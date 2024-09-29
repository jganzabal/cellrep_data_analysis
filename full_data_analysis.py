# %%
import scanpy as sc
import numpy as np
import re
import pandas as pd
# %%
adata = sc.read_10x_mtx(
    'data/GSE190604/',
    prefix='GSE190604_',
    var_names='gene_symbols',
    make_unique=False,
    cache=True,
    gex_only=False
)
# %%
adata.X.shape
# %%
# Extract the row indices, column pointers, and non-zero data
rows =  adata.X.indices     # Row indices of non-zero elements
cols =  adata.X.indptr      # Column pointers, indicating where each column starts in the row and data arrays
data =  adata.X.data        # Non-zero values
# %%
len(data)/(103805 * 36755), len(cols), len(rows)
# %%
(data == 0).sum()
# %%

# %%
# Get the column for each non-zero value by iterating through `indptr`
column_indices = []
for col in range(len(cols) - 1):  # Iterate over each column
    start_idx = cols[col]
    end_idx = cols[col + 1]
    column_indices.extend([col] * (end_idx - start_idx))  # Repeat column index for each non-zero value in the column

non_zero_entries = list(zip(rows, column_indices, data))
# %%
len(non_zero_entries)
# %%
# Print the non-zero entries in the format (row, column, value)
# %%
df_cell_lines = pd.DataFrame([{
    'cell_line': adata.obs.index[row],
    'gene': var_guide_calls[col],
    'rna_count':rna_count
} for row, col, rna_count in non_zero_entries
])
# %%
nonzero_indices = adata.X.nonzero()
# %%
len(adata)
# %%
adata.var['feature_types'].value_counts()
# %%
adata.var['guide_calls'] = adata.var['feature_types'] == 'CRISPR Guide Capture'
# %%
adata.var
# %%
# adata.obs['condition'] = adata.obs.index.map(lambda x: 'Nostim' if re.search(r'-[1-4]', x) else 'Stim')
# %%
# guide_calls = pd.read_csv("data/Perturb-seq/data/cellranger-guidecalls-aggregated-unfiltered.txt", sep="\t")
# guide_calls
# %%
# sc.pp.calculate_qc_metrics(
#     adata, qc_vars=["guide_calls"], inplace=True, log1p=True
# )
# %%
# adata.obs[['total_counts_guide_calls', 'total_counts']].loc['GGGAGATAGACCGTTT-1']
# %%
X_guide_calls = adata.X[:, adata.var['guide_calls']].copy()
var_guide_calls = list(adata.var[adata.var['guide_calls']].index.copy())

# %%
X_guide_calls = adata.X[:, :].copy()
var_guide_calls = list(adata.var[:].index.copy())
# %%
X_guide_calls.shape, len(var_guide_calls)
# %%
# %%
# Extract the row indices, column pointers, and non-zero data
rows = X_guide_calls.indices     # Row indices of non-zero elements
cols = X_guide_calls.indptr      # Column pointers, indicating where each column starts in the row and data arrays
data = X_guide_calls.data        # Non-zero values
# Get the column for each non-zero value by iterating through `indptr`
column_indices = []
for col in range(len(cols) - 1):  # Iterate over each column
    start_idx = cols[col]
    end_idx = cols[col + 1]
    column_indices.extend([col] * (end_idx - start_idx))  # Repeat column index for each non-zero value in the column

non_zero_entries = list(zip(rows, column_indices, data))
# Print the non-zero entries in the format (row, column, value)
# %%
df_cell_lines = pd.DataFrame([{
    'cell_line': adata.obs.index[row],
    'gene': var_guide_calls[col],
    'rna_count':rna_count
} for row, col, rna_count in non_zero_entries
])
# %%
df_cell_lines.set_index(['cell_line', 'gene']).loc[['CATCAAGGTTGTTGTG-8']]
# %%
