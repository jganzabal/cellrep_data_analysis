# %%
import pandas as pd

df_all = pd.read_csv('../data/GSE190846/GSE190846_supp_CD4_CRISPR_screens_read_counts.tsv', sep='\t')
df_all
# %%
base_cols = ['sgRNA', 'Gene']
df_all.columns = base_cols + ['_'.join(c.split('_')[1:5]) for c in df_all.columns[2:]]
# %%
df_all.shape
# %%
def average_reps():
    donors = [15, 16]
    hls = ['high', 'low']
    for d in donors:
        for hl in hls:
            print(d, hl)
            df_all[f'Donor{d}_IFNG_{hl}'] = (df_all[f'Donor{d}_IFNG_{hl}_Rep1'] + df_all[f'Donor{d}_IFNG_{hl}_Rep2'])/2
# %%
df_all.shape
# %%
df_all
# %%
input_folder = 'input_data'
nt_filename = f'{input_folder}/complementary_NO-TARGET.txt'
nt_df = df_all[df_all['Gene'].apply(lambda x: x=='NO-TARGET')]['sgRNA']
nt_df.to_csv(f'{nt_filename}', index=False, header=False)
nt_df
# %%

# %%
filename = f'{input_folder}/{cytokine}_count.txt'
df_filt.to_csv(filename, sep='\t', index=False)




# %%
df_cols = pd.DataFrame([c.split('_') + [c] for c in df_all.columns[2:]])
df_cols.sort_values(2)
# %%
cytokine = 'TNFa'
df_cols[df_cols[2] == cytokine]
# %%
# %%

base_cols = ['sgRNA', 'Gene']
filter_str = cytokine.lower()
df_filt = df_all[['sgRNA', 'Gene'] + [col for col in df_all.columns if filter_str in col.lower()]]
df_filt = df_filt[[col for col in df_filt.columns if 'Rep2' not in col]]
df_filt.columns = base_cols + ['_'.join(c.split('_')[1:4]) for c in df_filt.columns[2:]]
df_filt
# %%
# %%
input_folder = 'input_data'
nt_filename = f'{input_folder}/{cytokine}_NO-TARGET.txt'
nt_df = df_filt[df_filt['Gene'].apply(lambda x: x=='NO-TARGET')]['sgRNA']
nt_df.to_csv(f'{nt_filename}', index=False, header=False)
nt_df
# %%
filename = f'{input_folder}/{cytokine}_count.txt'
df_filt.to_csv(filename, sep='\t', index=False)
# %%
# %%
# cd cellrep/data_analysis/mageck_tests
# conda activate mageckenv2
# %%
mageck_command = f"""
mageck test -k {filename} \\
-t Donor15_TNFa_high,Donor16_TNFa_high \\
-c Donor15_TNFa_low,Donor16_TNFa_low \\
-n data_out/CRISPRa.{cytokine} \\
--norm-method none \\
--paired \\
--control-sgrna {nt_filename} \\
--keep-tmp
"""
print(mageck_command)
# %%
mageck_command = f"""
mageck test -k {filename} \\
-t Donor15_TNFa_high,Donor16_TNFa_high \\
-c Donor15_TNFa_low,Donor16_TNFa_low \\
-n data_out/CRISPRa.{cytokine} \\
--norm-method none \\
--paired \\
--control-sgrna {nt_filename}
"""

print(mageck_command)
# %%
mageck_command = f"""
mageck test -k {filename} \\
-t Donor15_IFNG_high,Donor16_IFNG_high \\
-c Donor15_IFNG_low,Donor16_IFNG_low \\
-n data_out/CRISPRa.{cytokine} \\
--norm-method none \\
--paired \\
--control-sgrna {nt_filename}
"""

print(mageck_command)

# %%
mageck count -k IL2.txt \
-t Sample05_Donor15_IL2_high_S6_R1_001.fastq.gz,Sample13_Donor16_IL2_high_S8_R1_001.fastq.gz \
-c Sample06_Donor15_IL2_low_S5_R1_001.fastq.gz,Sample14_Donor16_IL2_low_S7_R1_001.fastq.gz -n IL2 –trim-5 22,23,24,25,26,28,29,30
# %%
mageck test -k IFNG.txt \
-t Sample01_Donor15_IFNG_high_Rep1_S2_R1_001.fastq.gz,Sample09_Donor16_IFNG_high_Rep1_S4_R1_001.fastq.gz \
-c Sample03_Donor15_IFNG_low_Rep1_S1_R1_001.fastq.gz,Sample11_Donor16_IFNG_low_Rep1_S3_R1_001.fastq.gz -n IFNG –trim-5 22,23,24,25,26,28,29,30