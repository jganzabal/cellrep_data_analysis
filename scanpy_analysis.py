# %%
import scanpy as sc
import pandas as pd
import numpy as np
# %%
# adata = sc.read_h5ad("data/Perturb-seq/data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad")
adata = sc.read_h5ad("data/Perturb-seq/data/HuTcellsCRISPRaPerturbSeq_Resting.h5ad")
# %%
adata.obs
# %%
(adata.X >= 1).sum()/np.prod(adata.X.shape)
# %%
np.log10(216)
# %%
gene_indexes = (adata.var == 'ABCB10').values.reshape(-1)
cell_indexes = (adata.obs.index == 'GGGAGATAGACCGTTT-1')
# %%
adata.X[cell_indexes, gene_indexes]
# %%
adata.var['features'].apply(lambda x: x)
# %%
nonzero_elements = adata.X.data
nonzero_indices = adata.X.nonzero()
# %%
adata.X
# %%
# Create a DataFrame with nonzero elements and their indices
nonzero_df = pd.DataFrame({
    'row': nonzero_indices[0],
    'col': nonzero_indices[1],
    'RNA counts': nonzero_elements
})
# %%
