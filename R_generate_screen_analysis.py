# %%
# Import Packages ---------------------------------------------------------
import pandas as pd
import numpy as np
import os
import re
# %%
# Import Data -------------------------------------------------------------

# List files
gene_files = [os.path.join("data/Genome-wide-screens/data/gene_summaries", file) for file in os.listdir("data/Genome-wide-screens/data/gene_summaries") if file.endswith(".txt")]
sgrna_files = [os.path.join("data/Genome-wide-screens/data/sgRNA_summaries", file) for file in os.listdir("data/Genome-wide-screens/data/sgRNA_summaries") if file.endswith(".txt")]
# 
# Read data
gene_summaries = [pd.read_csv(file, sep='\t') for file in gene_files]
sgrna_summaries = [pd.read_csv(file, sep='\t') for file in sgrna_files]
# %%
gene_summaries

# %%

# Get Individual Donor Log2FoldChange Values ------------------------------

# Function to get median LFC
def get_median_lfc(df):
    # Separate sgrna into gene, sgrna, and donor
    df[['gene', 'sgrna', 'Donor']] = df['sgrna'].str.split('_', expand=True)
    df['sgrna'] = df['gene'] + '_' + df['sgrna']
    
    # Pivot donors into two columns
    df_pivot = df.pivot_table(index=['sgrna', 'Gene'], columns='Donor', values='LFC').reset_index()
    df_pivot.columns = ['sgrna', 'Gene', 'Donor1_LFC', 'Donor2_LFC']
    
    # Get the median LFCs
    median_lfcs = df_pivot.groupby('Gene').agg({'Donor1_LFC': 'median', 'Donor2_LFC': 'median'}).reset_index()
    
    return median_lfcs

donor_lfcs = [get_median_lfc(df) for df in sgrna_summaries]
# %%
donor_lfcs[4]
# %%
# Add Screen Column Identifiers --------------------------------------------

# Function to extract screen identifiers
def get_screen_identifier(files, suffix):
    return [re.sub(f".*?CRISPR(.*?)\\.{suffix}$", r"CRISPR\\1", file).replace('.', '_') for file in files]

# Add Screen column to donor LFCs
names_sgrna_dat = get_screen_identifier(sgrna_files, 'sgrna_summary.txt')
for i, df in enumerate(donor_lfcs):
    df['Screen'] = names_sgrna_dat[i]

# Add Screen column to gene summaries
names_gene_dat = get_screen_identifier(gene_files, 'gene_summary.txt')
for i, df in enumerate(gene_summaries):
    df['Screen'] = names_gene_dat[i]
# %%
df
# %%
# Get single FDR value ----------------------------------------------------

# Function to harmonize FDR values
def get_harmonized_fdr(df):
    df['FDR'] = df[['neg|fdr', 'pos|fdr']].min(axis=1)
    return df

gene_summaries = [get_harmonized_fdr(df) for df in gene_summaries]
# %%
gene_summaries[0]
# %%
# Get LFC Z-scores --------------------------------------------------------

# Function to calculate z-scores
def get_zscore(df):
    screen_id = df['Screen'].unique()[0]
    
    df['LFC'] = df['pos|lfc']
    
    # Reverse LFC for CRISPRi screens
    if 'CRISPRi' in screen_id:
        df['rev_lfc'] = -df['LFC']
        df['zscore'] = (df['rev_lfc'] - df['rev_lfc'].mean()) / df['rev_lfc'].std()
        df.drop(columns=['rev_lfc'], inplace=True)
    else:
        df['zscore'] = (df['LFC'] - df['LFC'].mean()) / df['LFC'].std()
    
    return df

gene_summaries = [get_zscore(df) for df in gene_summaries]
# %%
gene_summaries[0]
# %%
# Bind Screens into Single Dataframe --------------------------------------

# Combine gene summaries and donor LFC data
gene_df = pd.concat(gene_summaries, ignore_index=True)
donor_df = pd.concat(donor_lfcs, ignore_index=True)
# %%
gene_df.rename(columns={'id': 'Gene'}, inplace=True)
# Merge gene and donor data
full_dat = pd.merge(gene_df, donor_df, on=['Gene', 'Screen'], how='outer')

# Filter out "NO-TARGET"
full_dat = full_dat[full_dat['Gene'] != 'NO-TARGET']

# %%
full_dat
# %%
donor_df.columns
# %%
# Add Hit Call ------------------------------------------------------------

# Define hit calls based on FDR and LFC thresholds
full_dat['Hit'] = (full_dat['FDR'] < 0.05) & (full_dat['LFC'].abs() > 0.5)

# Define Hit_Type
full_dat['Hit_Type'] = np.where((full_dat['Hit']) & (full_dat['Screen'].str.contains('CRISPRi')) & (full_dat['LFC'] < 0),
                                'Positive Regulator', 'NA')
full_dat['Hit_Type'] = np.where((full_dat['Hit']) & (full_dat['Screen'].str.contains('CRISPRi')) & (full_dat['LFC'] > 0),
                                'Negative Regulator', full_dat['Hit_Type'])
full_dat['Hit_Type'] = np.where((full_dat['Hit']) & (full_dat['Screen'].str.contains('CRISPRa')) & (full_dat['LFC'] > 0),
                                'Positive Regulator', full_dat['Hit_Type'])
full_dat['Hit_Type'] = np.where((full_dat['Hit']) & (full_dat['Screen'].str.contains('CRISPRa')) & (full_dat['LFC'] < 0),
                                'Negative Regulator', full_dat['Hit_Type'])

# %%
# Clean up columns for export ---------------------------------------------

dat_for_export = full_dat[['Gene', 'Screen', 'Donor1_LFC', 'Donor2_LFC', 'LFC', 'zscore', 'FDR', 'Hit', 'Hit_Type']]

dat_for_export.head()
# %%

# %%
