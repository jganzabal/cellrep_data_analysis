# %%
%load_ext autoreload
%autoreload 2
# %%
import pandas as pd
from matplotlib import pyplot as plt
from process_helper import plot_volcano, generate_screen_analysis_table, plot_donors_scatter, read_data
# %%
a, b, c, d = read_data('data_out')
# %%
b[0]
# %%

# %%
_ = generate_screen_analysis_table(
    'data_out',
    hit_FDR_thres=0.05,
    hit_LFC_thres=0.5
)
# %%
screen_analysis_table = pd.read_csv('data_out/screen_analysis_table.csv')
# screen_analysis_table = pd.read_csv('../data/Genome-wide-screens/tableS2.csv')
screen_analysis_table
# %%
cytokine = 'TNFa'
out_folder = 'data_out'
screen = 'CRISPRa'
# %%
screen_analysis_table = screen_analysis_table[
    (screen_analysis_table['Screen'] == screen) & 
    (screen_analysis_table['Cytokine'] == cytokine) 
].set_index('Gene')
# %%
top_10_low = pd.read_csv('data_out/CRISPRa.TNFa.gene.low.txt', sep='\t')['group_id'][:10].to_list()
top_10_high = pd.read_csv('data_out/CRISPRa.TNFa.gene.high.txt', sep='\t')['group_id'][:10].to_list()
# %%
fig_1 = plot_donors_scatter(screen_analysis_table)
# %%

fig_2 = plot_vulcano(
    screen_analysis_table,
    x_label='LFC ' + '$TNRa^{hi}/TNRa^{low}$',
    y_label='$-Log_{10}(FDR)$',
    eps=1e-6,
    top_10_low=top_10_low,
    top_10_high=top_10_high
)
# ax.set_ylim([ymin, ymax])
# %%
