# Get top markers for each cluster

# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

# Import AnnData Objects ---------------------------------------------------

stim = sc.read_h5ad("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad")
resting = sc.read_h5ad("data/HuTcellsCRISPRaPerturbSeq_Resting.h5ad")

# Generate DataFrame with list column for testing ----------------------------

dat = pd.DataFrame({
    'condition': ['Resting', 'Stim'],
    'object': [resting, stim]
})
print(dat)

del resting, stim  # Free up memory

# Get Cluster Markers -----------------------------------------------------

def find_all_markers(adata):
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).melt()
    markers.columns = ['cluster', 'rank', 'gene']
    markers['scores'] = pd.DataFrame(adata.uns['rank_genes_groups']['scores']).melt()['value']
    markers['pvals'] = pd.DataFrame(adata.uns['rank_genes_groups']['pvals']).melt()['value']
    markers['pvals_adj'] = pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj']).melt()['value']
    markers['logfoldchanges'] = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges']).melt()['value']
    return markers

dat['markers'] = dat['object'].apply(find_all_markers)

print(dat)

markers_df = pd.concat([
    dat['markers'].iloc[0].assign(condition='Resting'),
    dat['markers'].iloc[1].assign(condition='Stim')
])

# Save markers
markers_df['condition'] = markers_df['condition'].replace('Stim', 'Re-stimulated')
markers_df.to_csv("data/clusterMarkers.txt", sep="\t", index=False)

