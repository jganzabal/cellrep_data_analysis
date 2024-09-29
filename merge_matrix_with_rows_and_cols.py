# %%
from scipy.io import mmread
import pandas as pd
import re
# %%
guide_calls = pd.read_csv("data/Perturb-seq/data/cellranger-guidecalls-aggregated-unfiltered.txt", sep="\t")
features = pd.read_csv('data/GSE190604_features.tsv', sep='\t', header=None)
features.columns = ['Ensembl Gene ID', 'Gene Name', 'Method']
bar_codes_df = pd.read_csv('data/GSE190604_barcodes.tsv', sep='\t', header=None)
bar_codes_df.columns=['Cell line']
# %%
# Add condition, guide target ("gene"), perturbed/control (crispr) columns
bar_codes_df['condition'] = bar_codes_df['Cell line'].apply(lambda x: 'Nostim' if re.search(r'-[1-4]', x) else 'Stim')
features['Gene Clean Name'] = features.apply(lambda x: x['Gene Name'][:-2] if x['Method']=='CRISPR Guide Capture' else x['Gene Name'], axis=1)
features['CRISPr'] = features['Gene Name'].apply(lambda x: 'NT' if 'NO-TARGET' in x else 'perturbed')
# %%
bar_codes_df
# %%
mat = mmread('data/GSE190604_matrix.mtx')
# %%
nonzero_elements = mat.data
nonzero_indices = mat.nonzero()
# Create a DataFrame with nonzero elements and their indices
nonzero_df = pd.DataFrame({
    'row': nonzero_indices[0],
    'col': nonzero_indices[1],
    'RNA counts': nonzero_elements
})
# %%
del mat
# %%
df_rows_merged = pd.merge(features[['Gene Name', 'Method', 'Gene Clean Name', 'CRISPr']], nonzero_df, right_on='row', left_index=True, how='right').drop(columns=['row'])
df_rows_merged
# %%
df_merged = pd.merge(bar_codes_df, df_rows_merged, right_on='col', left_index=True, how='right').drop(columns=['col'])
df_merged
# %%
df_merged.to_parquet('data_out/all_merged.parquet')
# %%
del df_rows_merged
# %%
df_indexed = df_merged.set_index(['Cell line', 'Gene Name', 'Method'])
del df_merged
# %%
df_indexed_cripr = df_indexed.loc[:, :, 'CRISPR Guide Capture']
# %%
guide_calls
# %%
df_indexed.loc['GGGAGATAGACCGTTT-1'].loc['ABCB10-1']
# %%
df_indexed.loc['GTCACTCAGCAGCCTC-8'].loc['WT1-2']
# %%
# 250
df_indexed.loc['GTCACTCAGCAGCCTC-8'].loc['VAV1-2']
# %%
df_indexed.loc['GTCACTCAGCAGCCTC-8', :, 'CRISPR Guide Capture'].sort_values('RNA counts', ascending=False).head(60)
# %%
guide_calls['feature_call'] = guide_calls.feature_call.apply(lambda x: x.split('|'))
guide_calls['num_umis'] = guide_calls.num_umis.apply(lambda x: x.split('|'))
# %%
guide_calls_exploded = guide_calls.explode(column=['feature_call', 'num_umis'])
# %%
guide_calls_exploded
# %%
guide_calls_exploded['num_umis'] = guide_calls_exploded['num_umis'].astype(int)
# %%
guide_calls_exploded_indexed = guide_calls_exploded.rename(columns={'cell_barcode': 'Cell line', 'feature_call': 'Gene Name'}).set_index(['Cell line', 'Gene Name'])
# %%
merged = guide_calls_exploded_indexed.merge(
    df_indexed_cripr, how='left', left_index=True, right_index=True
)
# %%
nan_ratio = merged['RNA counts'].isna().sum()/len(merged)
print(nan_ratio)
# %%
merged_clean = merged[~merged['RNA counts'].isna()]
# %%
exact_coincidence_ratio = (merged_clean['num_umis'].astype(float) == merged_clean['RNA counts']).sum()/len(merged_clean)
print(exact_coincidence_ratio)
# %%
error = 0.1
print((((merged_clean['num_umis'] - merged_clean['RNA counts'])/merged_clean['num_umis']) > error).sum()/len(merged_clean))
# %%
