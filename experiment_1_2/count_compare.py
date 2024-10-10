# %%
import pandas as pd
# %%
all_data_df = pd.read_excel(
    '../data/GSE174292/GSE174255	/GSE174255_sgRNA-Read-Counts.xlsx', sheet_name=None
)
# %%
dfs = []
for i, df in all_data_df.items():
    if 'CalabreseSet' in i:
        dfs.append(df)
# %%
df_1 = dfs[0].set_index(['sgRNA', 'Gene'])
df_2 = dfs[1].set_index(['sgRNA', 'Gene'])
# %%
df_1.sum()
# %%
df_2.sum()
# %%

# %%
raw_df = pd.read_csv(
    '../data/GSE190846/GSE190846_supp_CD4_CRISPR_screens_read_counts.tsv', 
    sep='\t', index_col=False
).set_index(['sgRNA', 'Gene'])
raw_df.columns
# %%
raw_df.sum(axis=1).std()/raw_df.sum(axis=1).mean()
# %%

# %%
raw_df.columns = ['sgRNA', 'Gene'] + ['_'.join(c.split('_')[1:5]) for c in raw_df.columns[2:]]
# %%
raw_df_1 = raw_df.rename(columns={
    'Donor15_IFNG_high_Rep1': 'Donor1_IFNG_high_raw_1', 'Donor15_IFNG_low_Rep1': 'Donor1_IFNG_low_raw_1',
    'Donor16_IFNG_high_Rep1': 'Donor2_IFNG_high_raw_1', 'Donor16_IFNG_low_Rep1': 'Donor2_IFNG_low_raw_1'
})[
    ['sgRNA', 'Gene',  'Donor1_IFNG_high_raw_1', 'Donor1_IFNG_low_raw_1', 'Donor2_IFNG_high_raw_1', 'Donor2_IFNG_low_raw_1']
].sort_values('sgRNA').reset_index(drop=True)
# %%
raw_df_2 = raw_df.rename(columns={
    'Donor15_IFNG_high_Rep2': 'Donor1_IFNG_high_raw_2', 'Donor15_IFNG_low_Rep2': 'Donor1_IFNG_low_raw_2',
    'Donor16_IFNG_high_Rep2': 'Donor2_IFNG_high_raw_2', 'Donor16_IFNG_low_Rep2': 'Donor2_IFNG_low_raw_2'
})[
    ['sgRNA', 'Gene',  'Donor1_IFNG_high_raw_2', 'Donor1_IFNG_low_raw_2', 'Donor2_IFNG_high_raw_2', 'Donor2_IFNG_low_raw_2']
].sort_values('sgRNA').reset_index(drop=True)
del raw_df
# %%
# raw_df = pd.read_csv('input_data/IFNG_count.txt', sep='\t').sort_values('sgRNA').reset_index(drop=True).rename(columns={
#     'Donor15_IFNG_high': 'Donor1_IFNG_high_raw', 'Donor15_IFNG_low': 'Donor1_IFNG_low_raw',
#     'Donor16_IFNG_high': 'Donor2_IFNG_high_raw', 'Donor16_IFNG_low': 'Donor2_IFNG_low_raw'
# })
# %%
orig_df = pd.read_csv('input_data/IFNG_count.orig.txt', sep='\t')[
    ['sgRNA', 'Gene',  'Donor1_IFNG_high', 'Donor1_IFNG_low', 'Donor2_IFNG_high', 'Donor2_IFNG_low']
].sort_values('sgRNA').reset_index(drop=True)

# %%
orig_df.shape, raw_df_1.shape
# %%
(orig_df['sgRNA'] == raw_df_1['sgRNA']).sum()/len(raw_df_1)
# %%
(orig_df['Gene'] == raw_df_1['Gene']).sum()/len(raw_df_1)
# %%
orig_df = orig_df.set_index(['sgRNA', 'Gene'])
raw_df_1 = raw_df_1.set_index(['sgRNA', 'Gene'])
raw_df_2 = raw_df_2.set_index(['sgRNA', 'Gene'])
# %%
orig_df.sum(axis=1)
# %%
raw_df_1.sum(axis=1)
# %%
raw_df_2.sum(axis=1)
# %%
orig_df.iloc[0].values
# %%
i = 6
(raw_df_1.iloc[i].values + raw_df_2.iloc[i].values)/orig_df.iloc[i].values
# %%

# %%
df_all = pd.concat([orig_df, raw_df_1, raw_df_2], axis=1)
# %%
df_all.columns
# %%
df_all.corr()
# %%
df_all.describe()
# %%
orig_df.sum(axis=1).std()/ orig_df.sum(axis=1).mean()
# %%
raw_df_1.sum(axis=1).std()/raw_df_1.sum(axis=1).mean()
# %%
raw_df_2.sum(axis=1).std()/raw_df_2.sum(axis=1).mean() 
# %%
raw_df_2.sum()
# %%
orig_df.sum()
# %%
