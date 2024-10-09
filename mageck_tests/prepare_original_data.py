# %%
import pandas as pd
# %%
all_data_df = pd.read_excel(
    '../data/GSE174292/GSE174255	/GSE174255_sgRNA-Read-Counts.xlsx', sheet_name=None
)
# %%
dfs = []
for k, df in all_data_df.items():
    print(k)
    if 'CalabreseSet' in k:
        df = df.copy()
        df.columns = df.columns[:2].tolist() + ['_'.join(c.split('_')[2:]) for c in df.columns[2:]]
        dfs.append(df)
# %%
raw_data = pd.concat(dfs, axis=0)
# %%
# %%
input_folder = 'input_data'
cytokine = 'IFNG'
nt_filename = f'{input_folder}/{cytokine}_NO-TARGET.orig.txt'
nt_df = raw_data[raw_data['Gene'].apply(lambda x: x=='NO-TARGET')]['sgRNA']
nt_df.to_csv(f'{nt_filename}', index=False, header=False)
len(nt_df)
# %%
# %%
filename = f'{input_folder}/{cytokine}_count.orig.txt'
raw_data.to_csv(filename, sep='\t', index=False)
# %%
raw_data.columns
# %% Original Data
mageck_command = f"""
mageck test -k {filename} \\
-t Donor1_IFNG_high,Donor2_IFNG_high \\
-c Donor1_IFNG_low,Donor2_IFNG_low \\
-n data_out/CRISPRa.{cytokine} \\
--norm-method none \\
--paired \\
--control-sgrna {nt_filename}
"""
print(mageck_command)
# %%
# %% Raw Data
# Donor15_IFNG_high	Donor15_IFNG_low	Donor16_IFNG_high	Donor16_IFNG_low
mageck_command = f"""
mageck test -k input_data/IFNG_count.txt \\
-t Donor15_IFNG_high,Donor16_IFNG_high \\
-c Donor15_IFNG_low,Donor16_IFNG_low \\
-n data_out/CRISPRa.{cytokine} \\
--norm-method none \\
--paired \\
--control-sgrna input_data/IFNG_NO-TARGET.txt
"""
print(mageck_command)
# %%
# cd cellrep/data_analysis/mageck_tests
# conda activate mageckenv2