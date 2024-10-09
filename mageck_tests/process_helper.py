import os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from adjustText import adjust_text

def get_median_lfc(df):
    # Each row is a single sgRNA/Donor measurement
    # sgRNA column contains both sgrna and donor ids --> split into 2 cols
    # r0 = Donor1, r1 = Donor2
    df[['gene', 'sgrna', 'Donor']] = df['sgrna'].apply(lambda x: pd.Series(x.split('_')))
    df['sgrna'] = df['gene'] + '_' + df['sgrna']
    
    # Spread donors into 2 columns
    df = df[['sgrna', 'Gene', 'Donor', 'LFC']].pivot_table(index=['sgrna', 'Gene'], columns='Donor', values='LFC', aggfunc='median').reset_index()
    df.columns = ['sgrna', 'Gene', 'Donor1_LFC', 'Donor2_LFC']
    
    # Get median
    median_lfcs = df[['Gene', 'Donor1_LFC', 'Donor2_LFC']].groupby('Gene').median().reset_index()
    return median_lfcs
# %%
def read_data(folder='data_out'):
    # Load Data
    gene_files = [f for f in os.listdir(f"{folder}") if f.endswith('.txt') and 'gene_summary' in f]
    gene_files = [os.path.join(f"{folder}/", f) for f in gene_files]
    gene_summaries = [pd.read_csv(f, sep='\t') for f in gene_files]
    print(gene_files)

    sgrna_files = [f for f in os.listdir(f"{folder}") if f.endswith('.txt') and 'sgrna_summary' in f]
    sgrna_files = [os.path.join(f"{folder}/", f) for f in sgrna_files]
    sgrna_summaries = [pd.read_csv(f, sep='\t') for f in sgrna_files]
    print(sgrna_files)
    names_sgrna_dat = [(f.split('.')[1], f.split('.')[0].split('/')[1]) for f in sgrna_files]
    names_gene_dat = [(f.split('.')[1], f.split('.')[0].split('/')[1]) for f in gene_files]
    return gene_summaries, sgrna_summaries, names_gene_dat, names_sgrna_dat
# %%
def get_harmonized_fdr(df):
    df['FDR'] = df.apply(lambda x: min(x['neg|fdr'], x['pos|fdr']), axis=1)
    return df

def get_zscore(df):
    screen_id = df['Screen'].iloc[0]
    df['LFC'] = df['pos|lfc']
    
    # Reverse LFC for CRISPRi screens so positive regulators have positive 
    # z-scores (e.g. CD28), which will be useful for comparing screens
    if 'CRISPRi' in screen_id:
        df['rev_lfc'] = -1 * df['LFC']
        df['zscore'] = (df['rev_lfc'] - df['rev_lfc'].mean()) / df['rev_lfc'].std()
    elif 'CRISPRa' in screen_id:
        df['zscore'] = (df['LFC'] - df['LFC'].mean()) / df['LFC'].std()
    
    if 'rev_lfc' in df.columns:
        df.drop(columns=['rev_lfc'], inplace=True)
    return df

def generate_screen_analysis_table(folder='data_out', hit_FDR_thres=0.05, hit_LFC_thres=0.5, save_data=True):
    # Load data
    gene_summaries, sgrna_summaries, names_gene_dat, names_sgrna_dat = read_data(folder=folder)

    # Gert median
    donor_lfcs = [get_median_lfc(df) for df in sgrna_summaries]
    
    # Add columns (cytokine and scrren)
    donor_lfcs = [df.assign(Screen=screen, Cytokine=cytokine) for df, (cytokine, screen) in zip(donor_lfcs, names_sgrna_dat)]
    gene_summaries = [df.assign(Screen=screen, Cytokine=cytokine) for df, (cytokine, screen) in zip(gene_summaries, names_gene_dat)]

    print(f'Both pos and neg fdr minor than 0.05')
    filtered_gene_summaries = [df[(df['neg|fdr'] < 0.05) & (df['pos|fdr'] < 0.05)] for df in gene_summaries]
    print([len(d) for d in filtered_gene_summaries])

    # Get single FDR value
    gene_summaries = [get_harmonized_fdr(df) for df in gene_summaries]

    # Get LFC Z-scores
    gene_summaries = [get_zscore(df) for df in gene_summaries]
    
    # Bind Screens into Single Dataframe
    gene_df = pd.concat(gene_summaries, ignore_index=True).rename(columns={'id': 'Gene'})
    donor_df = pd.concat(donor_lfcs, ignore_index=True)
    
    # Join with Donor LFC data
    full_dat = pd.merge(gene_df, donor_df, on=['Gene', 'Screen', 'Cytokine'], how='outer')
    
    # Filter out no-taget measurements (not useful at gene level and returns NAs)
    full_dat = full_dat[full_dat['Gene'] != 'NO-TARGET']
    
    # Add Hit Call
    full_dat['Hit'] = (full_dat['FDR'] < hit_FDR_thres) & (abs(full_dat['LFC']) > hit_LFC_thres)
    
    full_dat['Hit_Type'] = full_dat.apply(lambda x: 'Positive Regulator' if x['Hit'] and 'CRISPRi' in x['Screen'] and x['LFC'] < 0 else 'NA', axis=1)
    full_dat['Hit_Type'] = full_dat.apply(lambda x: 'Negative Regulator' if x['Hit'] and 'CRISPRi' in x['Screen'] and x['LFC'] > 0 else x['Hit_Type'], axis=1)
    full_dat['Hit_Type'] = full_dat.apply(lambda x: 'Positive Regulator' if x['Hit'] and 'CRISPRa' in x['Screen'] and x['LFC'] > 0 else x['Hit_Type'], axis=1)
    full_dat['Hit_Type'] = full_dat.apply(lambda x: 'Negative Regulator' if x['Hit'] and 'CRISPRa' in x['Screen'] and x['LFC'] < 0 else x['Hit_Type'], axis=1)
    
    # Clean up columns for export
    dat_for_export = full_dat[['Gene', 'Screen', 'Donor1_LFC', 'Donor2_LFC', 'LFC', 'zscore', 'FDR', 'Hit', 'Hit_Type', 'Cytokine']]

    if save_data:
        dat_for_export.to_csv(f'{folder}/screen_analysis_table.csv')
    return dat_for_export


def plot_donors_scatter(screen_analysis_table):

    not_a_hit_label = 'Not a Hit'
    screen_analysis_table['Hit_Type'] = screen_analysis_table['Hit_Type'].fillna(not_a_hit_label)

    color_dict = {
        not_a_hit_label: 'gray',
        'Negative Regulator': 'blue',
        'Positive Regulator': 'red',
    }

    screen_analysis_table['color'] = screen_analysis_table['Hit_Type'].apply(color_dict.get)

    fig, ax = plt.subplots()
    ax.scatter(
        screen_analysis_table['Donor1_LFC'], 
        screen_analysis_table['Donor2_LFC'],
        c=screen_analysis_table['color'],
        s=2
    )
    texts = []
    for gene, row in screen_analysis_table.iterrows():
        x = row['Donor1_LFC']
        y = row['Donor2_LFC']
        if (x**2 + y**2)**(1/2) > 2:
            
            texts.append(
                ax.text(
                    x, y, gene,
                    # color=GREY30,
                    fontsize=8,
                    fontname="Poppins"
                )
            )
            # ax.annotate(gene, (x, y))
    adjust_text(
        texts, 
        expand_points=(2, 4),
        arrowprops=dict(
            arrowstyle="->", 
            # color=GREY50, 
            lw=0.5
        ),
        ax=fig.axes[0]
    )
    ax.set_xlabel('Donor 15, LFC')
    ax.set_ylabel('Donor 16, LFC')
    return fig

def plot_vulcano(screen_analysis_table, x_label='', y_label='', eps=1e-6, top_10_low=[], top_10_high=[]):
    genes_to_show = top_10_low + top_10_high
    screen_analysis_table['color_volcano'] = screen_analysis_table.index.map(lambda x: 'red' if x in top_10_high else 'blue' if x in top_10_low else 'gray')


    screen_analysis_table['log_10_FDR'] = -np.log10(eps + screen_analysis_table['FDR'])
    # To plot them after
    screen_analysis_table_filt = screen_analysis_table[
        screen_analysis_table.index.isin(genes_to_show)
    ]
    
    fig, ax = plt.subplots()

    ax.scatter(
        screen_analysis_table['LFC'],
        screen_analysis_table['log_10_FDR'],
        c=screen_analysis_table['color_volcano'],
        s=1
    )

    ax.scatter(
        screen_analysis_table_filt['LFC'],
        screen_analysis_table_filt['log_10_FDR'], 
        c=screen_analysis_table_filt['color_volcano'],
        s=2
    )
    plt.axvline(-2,color="grey",linestyle="--")
    plt.axvline(2,color="grey",linestyle="--")
    plt.axhline(2,color="grey",linestyle="--")


    texts = []
    for gene, row in screen_analysis_table_filt.iterrows():
        if gene in genes_to_show:
            x = row['LFC']
            y = row['log_10_FDR']
            # ax.annotate(gene, (x, y))
            texts.append(
                ax.text(
                    x, y, gene,
                    # color=GREY30,
                    fontsize=8,
                    fontname="Poppins"
                )
            )
    adjust_text(
        texts, 
        expand_points=(2, 4),
        arrowprops=dict(
            arrowstyle="->", 
            # color=GREY50, 
            lw=0.5
        ),
        ax=fig.axes[0]
    )
    ax.set_xlim([-3, 3])
    ax.set_ylim([-0.1, 4.5])
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    return fig

