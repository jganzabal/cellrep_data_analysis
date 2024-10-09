#%%
# Generate log2(CD4/CD8) expression score for each cell

# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %%
# Import AnnData Object ----------------------------------------------------

tcells = sc.read_h5ad("data/tcells_all_filt.h5ad")

# %% Get CD4/CD8 score for each cell -----------------------------------------

# Extract dataframe of CD4 and CD8A/CD8B gene expression
CD4_8_df = tcells[:, ['CD4', 'CD8A', 'CD8B']].to_df()
display(CD4_8_df.head())

# For CD8 average CD8A and CD8B expression for each cell
CD4_8_df['CD8'] = CD4_8_df[['CD8A', 'CD8B']].mean(axis=1)
display(CD4_8_df.head())

min_cd4 = CD4_8_df.loc[CD4_8_df['CD4'] != 0, 'CD4'].min()
print(f"Min CD4 value != 0 is {min_cd4:.2f}")

min_cd8 = CD4_8_df.loc[CD4_8_df['CD8'] != 0, 'CD8'].min()
print(f"Min CD8 value != 0 is {min_cd8:.2f}")

# %%
# Calculate CD4/CD8 score (log2FC). Use minimum value/2 as pseudocount
pseudocount = CD4_8_df.loc[CD4_8_df['CD8'] != 0, 'CD8'].min() / 2

CD4_8_df['CD4'] += pseudocount
CD4_8_df['CD8'] += pseudocount
CD4_8_df['CD4.CD8.Score'] = np.log2(CD4_8_df['CD4'] / CD4_8_df['CD8'])
display(CD4_8_df.head())

# Add CD4.CD8.Score to AnnData object obs
tcells.obs['CD4.CD8.Score'] = CD4_8_df['CD4.CD8.Score']
display(tcells.obs.head())
# %%
# %% Plot CD4/CD8 Scores --------------------------------------------------

# Distribution of scores
fig, axs = plt.subplots(2, 3, figsize=(12, 10))
for i, condition in enumerate(tcells.obs['condition'].unique()):
    for j, donor in enumerate(tcells.obs['donor'].unique()):
        # print(i, j)
        sns.kdeplot(data=tcells.obs[(tcells.obs['condition'] == condition) & 
                                    (tcells.obs['donor'] == donor)], 
                    x='CD4.CD8.Score', ax=axs[i, j])
        axs[i, j].set_title(f"{condition} - {donor}")
plt.tight_layout()
plt.show()
# %%
# Set cutoffs for CD4/CD8 assignment
fig, ax = plt.subplots(figsize=(10, 6))
for condition in tcells.obs['condition'].unique():
    sns.kdeplot(data=tcells.obs[tcells.obs['condition'] == condition], 
                x='CD4.CD8.Score', ax=ax, label=condition)
ax.axvline(-0.9, linestyle='--', color='black')
ax.axvline(1.4, linestyle='--', color='black')
ax.legend()
ax.set_title('CD4/CD8 Score Distribution')
plt.show()
# %%
# Add CD4/CD8 assignments to obs -------------------------------------

tcells.obs['CD4_or_CD8'] = pd.cut(tcells.obs['CD4.CD8.Score'], 
                                  bins=[-np.inf, -0.9, 1.4, np.inf], 
                                  labels=['CD8', 'unassigned', 'CD4'])

print(tcells.obs.head())

# Save AnnData Object ------------------------------------------------------

tcells.write_h5ad("data/tcells_all_filt.h5ad", compression="gzip")

# %%
tcells