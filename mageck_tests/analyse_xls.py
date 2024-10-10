# %%
import pandas as pd
# %%
df = pd.read_excel(
    '../data/Genome-wide-screens/data/science.abj4008_table_s2.xlsx', 
)
# %%
df.groupby(['Screen_Version', 'CRISPRa_or_i', 'CD4_or_CD8', 'Cytokine']).size()
# %%
df.columns
# %%
