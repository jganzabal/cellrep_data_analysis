# Donor demultiplexing metadata from souporcell

# %%
#  Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce

# %%
# Import AnnData object ---------------------------------------------------

tcells = sc.read_h5ad("data/tcells_all_filt.h5ad")

#%% 
# Import SoupOrCell Data --------------------------------------------------

import os
data_path = "../data/Perturb-seq/data/souporcell/"
souporcell = pd.DataFrame({'sample': os.listdir(data_path)})
souporcell = souporcell[souporcell['sample'] != '.DS_Store']
souporcell
# %%
souporcell['directory'] = data_path + souporcell['sample'] + "/clusters.tsv"
souporcell['data'] = souporcell['directory'].apply(lambda x: pd.read_csv(x, sep='\t'))
souporcell
# %%
souporcell['data'].iloc[0]
# %%
# Fix well identifiers
souporcell['well_id'] = [str(i) for i in range(1, 9)]

def fix_barcode(df, well_id):
    df['barcode'] = df['barcode'].str[:-1] + well_id
    return df

souporcell['data'] = [fix_barcode(df, well_id) for df, well_id in zip(souporcell['data'], souporcell['well_id'])]
# %%
souporcell
# %%
# Combine to single dataframe
souporcell = pd.concat(souporcell['data'].tolist(), ignore_index=True)
# %%
souporcell.status.value_counts()
# %%
# SoupOrCell Donor Values -------------------------------------------------

# Need to harmonize donor values across wells, as "0" or "1" assignment is 
# arbitrary each time. Donors were matched using output .vcf file

souporcell_donor_calls = pd.read_csv("../data/Perturb-seq/data/souporcell_match_vcf_res/donor_calls.txt", sep='\t')
# %%
# Add to souporcell dataframe
souporcell['Well_ID'] = souporcell['barcode'].str[-1].astype(int)
souporcell = souporcell.merge(souporcell_donor_calls, on='Well_ID')

# Add donor assignments to each droplet
souporcell['Donor'] = np.where(souporcell['assignment'] == '0', 
                               souporcell['Souporcell_call0_DonorA_or_B'], 
                               souporcell['Souporcell_call1_DonorA_or_B'])
souporcell.loc[souporcell['status'] == 'unassigned', 'Donor'] = 'unassigned'
# %%
souporcell
# %%
len(souporcell), len(tcells.obs)
# %%
# Check that all cell barcodes in AnnData object are present in souporcell df
assert all(tcells.obs.index.isin(souporcell['barcode']))
# %%
# Add to metadata
souporcell_filtered = souporcell[souporcell['barcode'].isin(tcells.obs.index)]
tcells.obs['souporcell_donor'] = souporcell_filtered.set_index('barcode')['Donor']
# %%
tcells
# %%
# Check male/female donors using RPS4Y1 (Y chromosome gene)
# From donor info sheets donor 1 was female, and donor 2 was male

sc.pl.violin(tcells, "RPS4Y1", groupby="souporcell_donor")
# %%
# Donor "A" is likely the male, so reassign as "Donor2"  

tcells.obs['donor'] = np.where(tcells.obs['souporcell_donor'] == 'A', 'Donor2', 'unassigned')
tcells.obs.loc[tcells.obs['souporcell_donor'] == 'B', 'donor'] = 'Donor1'
tcells.obs = tcells.obs.drop(columns=['souporcell_donor'])

# %%
tcells.obs.head()
# %%
# Save AnnData Object -----------------------------------------------------

tcells.write_h5ad("data/tcells_all_filt.h5ad", compression="gzip")
