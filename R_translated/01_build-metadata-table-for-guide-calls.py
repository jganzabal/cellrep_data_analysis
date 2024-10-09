# %%
import pandas as pd
import re
# %%
# Import data -------------------------------------------------------------

# Guide calls from CellRanger
guide_calls = pd.read_csv("data/Perturb-seq/data/cellranger-guidecalls-aggregated-unfiltered.txt", sep="\t")
guide_calls
# %%
guide_calls['num_features'].sum()
# %%
# %%
# Build metadata table ----------------------------------------------------
metadata = guide_calls[guide_calls['num_features'] == 1].copy()
metadata
# %%
# Filter for number of UMIs > 5
metadata['num_umis'] = metadata['num_umis'].astype(float)
metadata = metadata[metadata['num_umis'] >= 5]
metadata
# %%
# Add condition, guide target ("gene"), perturbed/control (crispr) columns
metadata['condition'] = metadata['cell_barcode'].apply(lambda x: 'Nostim' if re.search(r'-[1-4]', x) else 'Stim')
metadata['gene'] = metadata['feature_call'].str[:-2]
metadata['crispr'] = metadata['gene'].apply(lambda x: 'NT' if 'NO-TARGET' in x else 'perturbed')
metadata
# %%
# Select and rename columns
metadata = metadata[['cell_barcode', 'condition', 'crispr', 'feature_call', 'gene']]
metadata = metadata.rename(columns={'feature_call': 'guide_id'})
# %%
metadata
# %%
# Save metadata table -----------------------------------------------------

metadata.to_csv("data/cell_metadata.txt", sep="\t", index=False)
# %%
