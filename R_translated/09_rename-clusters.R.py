# Reorder/name clusters based on dendogram (tree)

# %% 
# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage

# %%
# Import Data -------------------------------------------------------------

dat = pd.DataFrame({
    'condition': ['NoStim', 'Stim'],
    'object': [
        sc.read_h5ad("data/Nostim.h5ad"),
        sc.read_h5ad("data/Stim.h5ad")
    ]
})
print(dat)

# %% Get cluster order for plotting ------------------------------------------

# Use the order they appear in the cluster dendogram
# First build dendograms
def build_cluster_tree(adata):
    sc.tl.dendrogram(adata, groupby='leiden')
    return adata

dat['object'] = dat['object'].apply(build_cluster_tree)

# Now get tree order

def get_tree_order(adata):
    dendro = adata.uns['dendrogram_leiden']['dendrogram_info']
    leaves = dendro['leaves']
    return [str(leaf) for leaf in leaves]

dat['tree_order'] = dat['object'].apply(get_tree_order)
print(dat)
# %%
dat['object'].iloc[1]
# %% NoStim
sc.pl.dendrogram(dat['object'].iloc[0], groupby='leiden')
print(dat['tree_order'].iloc[0])

# %% Stim
sc.pl.dendrogram(dat['object'].iloc[1], groupby='leiden')
print(dat['tree_order'].iloc[1])

# %% Relevel cluster metadata

dat['object'].iloc[0].obs['leiden'] = pd.Categorical(
    dat['object'].iloc[0].obs['leiden'],
    categories=dat['tree_order'].iloc[0],
    ordered=True
)
print(dat['object'].iloc[0].obs['leiden'].head())

# %%
dat['object'].iloc[1].obs['leiden'] = pd.Categorical(
    dat['object'].iloc[1].obs['leiden'],
    categories=dat['tree_order'].iloc[1],
    ordered=True
)
display(dat['object'].iloc[1].obs['leiden'].head())

# %% Cluster conversion tables
nostim_cluster_conversion = pd.DataFrame({
    'original_id': dat['object'].iloc[0].obs['leiden'].cat.categories,
    'new_id': [str(i+1) for i in range(len(dat['object'].iloc[0].obs['leiden'].cat.categories))]
})
display(nostim_cluster_conversion)

# %%
stim_cluster_conversion = pd.DataFrame({
    'original_id': dat['object'].iloc[1].obs['leiden'].cat.categories,
    'new_id': [str(i+1) for i in range(len(dat['object'].iloc[1].obs['leiden'].cat.categories))]
})
display(stim_cluster_conversion)
# %%
cluster_conversion = pd.concat([
    nostim_cluster_conversion.assign(condition='NoStim'),
    stim_cluster_conversion.assign(condition='Stim')
])
print(cluster_conversion)
# %%
# %% Convert to new cluster ids based on tree order
nostim_new_cluster_ids = dict(zip(nostim_cluster_conversion['original_id'], nostim_cluster_conversion['new_id']))
dat['object'].iloc[0].obs['leiden'] = dat['object'].iloc[0].obs['leiden'].map(nostim_new_cluster_ids)
sc.pl.umap(dat['object'].iloc[0], color='leiden', legend_loc='on data', title='NoStim')
# %%
stim_new_cluster_ids = dict(zip(stim_cluster_conversion['original_id'], stim_cluster_conversion['new_id']))
dat['object'].iloc[1].obs['leiden'] = dat['object'].iloc[1].obs['leiden'].map(stim_new_cluster_ids)
sc.pl.umap(dat['object'].iloc[1], color='leiden', legend_loc='on data', title='Stim')
