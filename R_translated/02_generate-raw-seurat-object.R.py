# Generate AnnData object and add metadata
# %%
# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np

# Import Data -------------------------------------------------------------
# %%
# Full aggregated data from 8 chromium Wells
aggr_dat = sc.read_10x_mtx(
    '../data/GSE190604/',
    prefix='GSE190604_',
    cache=True,
    gex_only=True, # Only gene expression
    # make_unique=True
)
# %%
# Cell metadata (prefiltered guide-calls for singlets)
cell_metadata = pd.read_csv("../data/cell_metadata.txt", sep="\t").set_index("cell_barcode")
cell_metadata
# %%
singlet_barcodes = cell_metadata.index
singlet_barcodes
# %%
aggr_dat.obs_names
# %%
# Prepare gene expression count data --------------------------------------

# Guide singlet cell filter (prefiltered in metadata table)

aggr_dat = aggr_dat[aggr_dat.obs_names.isin(singlet_barcodes), :]
print(aggr_dat)
aggr_dat
# Generate AnnData Object --------------------------------------------------
# %%
tcells = sc.AnnData(X=aggr_dat.X, 
                    obs=cell_metadata.loc[aggr_dat.obs_names],
                    var=aggr_dat.var)

sc.pp.filter_cells(tcells, min_genes=200)
sc.pp.filter_genes(tcells, min_cells=3)

print(tcells)
print(tcells.obs.head())
# %%
tcells
# Save Data ---------------------------------------------------------------
# %%
tcells.write_h5ad("data/tcells_all_unfilt.h5ad")
