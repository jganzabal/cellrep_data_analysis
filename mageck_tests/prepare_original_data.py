# %%
import pandas as pd

def read_cd8_data_for_mageck(filename= '../data/GSE174292/GSE174255/GSE174255_sgRNA-Read-Counts.xlsx', crispr_type='CRISPRa'):
    all_data_df = pd.read_excel(
       filename, sheet_name=None
    )
    if crispr_type=='CRISPRa':
        dfs = []
        for k, df in all_data_df.items():
            if 'CalabreseSet' in k:
                print(k)
                df = df.copy()
                df.columns = df.columns[:2].tolist() + ['_'.join(c.split('_')[2:]) for c in df.columns[2:]]
                dfs.append(df)
    else:
        dfs = []
        for k, df in all_data_df.items():
            if 'CalabreseSet' not in k:
                print(k)
                df = df.copy()
                df.columns = df.columns[:2].tolist() + ['_'.join(c.split('_')[2:]) for c in df.columns[2:]]
                dfs.append(df)
    return pd.concat(dfs, axis=0)

def data_to_file_for_mageck(raw_data, cell_type, crispr_type, input_folder='input_data'):
    nt_filename = f'{input_folder}/{crispr_type}.{cell_type}_NO-TARGET.orig.txt'
    nt_df = raw_data[raw_data['Gene'].apply(lambda x: x=='NO-TARGET')]['sgRNA']
    nt_df.to_csv(f'{nt_filename}', index=False, header=False)

    filename = f'{input_folder}/{crispr_type}.{cell_type}_count.orig.txt'
    raw_data.to_csv(filename, sep='\t', index=False)

# %%
# cytokine = 'IFNG'
# # %% Original Data
# mageck_command = f"""
# mageck test -k {filename} \\
# -t Donor1_IFNG_high,Donor2_IFNG_high \\
# -c Donor1_IFNG_low,Donor2_IFNG_low \\
# -n data_out/CRISPRa.{cytokine} \\
# --norm-method none \\
# --paired \\
# --control-sgrna {nt_filename}
# """
# print(mageck_command)
# # %%
# # %% Raw Data
# # Donor15_IFNG_high	Donor15_IFNG_low	Donor16_IFNG_high	Donor16_IFNG_low
# mageck_command = f"""
# mageck test -k input_data/IFNG_count.txt \\
# -t Donor15_IFNG_high,Donor16_IFNG_high \\
# -c Donor15_IFNG_low,Donor16_IFNG_low \\
# -n data_out/CRISPRa.{cytokine} \\
# --norm-method none \\
# --paired \\
# --control-sgrna input_data/IFNG_NO-TARGET.txt
# """
# print(mageck_command)
# # %%
# # cd cellrep/data_analysis/mageck_tests
# # conda activate mageckenv2