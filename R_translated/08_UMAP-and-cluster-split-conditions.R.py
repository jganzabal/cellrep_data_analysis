# UMAP and cluster for downstream analysis
# Picked dimensions and clustering algorithms by examining outputs of 
# various numbers of dimensions and clustering algorithms

# NoStim and Stim will be split for separate analyses

# %% 
# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %% Import AnnData Object ---------------------------------------------------

tcells = sc.read_h5ad("data/tcells_all_filt.h5ad")
tcells
# %% Split into two objects based on conditions ------------------------------

nostim = tcells[tcells.obs['condition'] == 'Nostim']
stim = tcells[tcells.obs['condition'] == 'Stim']
del tcells

# %%
#  UMAP and Cluster --------------------------------------------------------

sc.pp.pca(nostim, n_comps=50)
sc.pp.neighbors(nostim, n_neighbors=15, n_pcs=20)
sc.tl.umap(nostim)
sc.tl.leiden(nostim, resolution=0.5)

# %%
sc.pp.pca(stim, n_comps=50)
sc.pp.neighbors(stim, n_neighbors=15, n_pcs=20)
sc.tl.umap(stim)
sc.tl.leiden(stim, resolution=0.4)

# %% Plots -------------------------------------------------------------------

sc.pl.umap(nostim, color='leiden', legend_loc='on data', title='NoStim')
sc.pl.umap(stim, color='leiden', legend_loc='on data', title='Stim')

# %%
# Combine similar clusters ------------------------------------------------

# Combine Stim clusters 2 and 15. In downstream analysis it was shown that
# these two clusters were very similar. Cluster 2 represents the bulk of cells
# with negative regulator guides, and cluster 15 is exclusively MUC1 guides

# Rename cluster 15 to 2
# stim.obs['leiden'] = stim.obs['leiden'].cat.rename_categories({'15': '2'}) # Esto no tiene sentido

sc.pl.umap(stim, color='leiden', legend_loc='on data', title='Stim (after merging)')

# %% Save data ---------------------------------------------------------------

nostim.write_h5ad("data/Nostim.h5ad")
stim.write_h5ad("data/Stim.h5ad")

