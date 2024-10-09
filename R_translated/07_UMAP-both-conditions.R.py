# %% Generate UMAP plot for all cells (both restim/resting conditions)

# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %% Import AnnData object ----------------------------------------------------

tcells = sc.read_h5ad("data/tcells_all_filt.h5ad")
print(tcells)

# %% UMAP and Cluster --------------------------------------------------------

sc.tl.pca(tcells, n_comps=100)
sc.pp.neighbors(tcells, n_neighbors=15, n_pcs=35)
sc.tl.umap(tcells)
sc.tl.leiden(tcells, resolution=1.2)

# %% Rasterized condition UMAP plot for manuscript ---------------------------

umap_df = pd.DataFrame(tcells.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'], index=tcells.obs_names)
umap_df['cell'] = umap_df.index
umap_df = umap_df.merge(tcells.obs, left_index=True, right_index=True)
umap_df
# %%

np.random.seed(1)
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(tcells, color='condition', ax=ax, show=False)
ax.set_xlabel('UMAP 1')
ax.set_ylabel('UMAP 2')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='both', length=0)
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.legend(title='Condition', loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
plt.tight_layout()
plt.show()
# %%
