# The goal is to assign an unbiased T cell activation score to each cell, 
# based on the expression of differentially expressed genes upon T cell activation.
# The No-Target control cells will be used to identify these differentially expressed genes.

# Importing required packages
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
import seaborn as sns
import matplotlib.pyplot as plt

# Load the Seurat object with both restim and resting conditions ----
tcells = sc.read_h5ad("data/tcells_all_filt.h5ad")  # Loading the AnnData object

# Differential expression of No-Target control cells: stim vs nostim ----
# Filtering for NT (No-Target) control cells
ctrl_cells = tcells[tcells.obs['crispr'] == 'NT']

# Ensure 'condition' is a categorical variable and reorder it for differential expression analysis
ctrl_cells.obs['condition'] = pd.Categorical(ctrl_cells.obs['condition'], categories=['Stim', 'Nostim'])

# Performing differential expression analysis using Wilcoxon rank-sum test
sc.tl.rank_genes_groups(ctrl_cells, 'condition', method='wilcoxon', min_in_group_fraction=0.25)
stim_de = sc.get.rank_genes_groups_df(ctrl_cells, group='Stim')

# Display the top upregulated and downregulated genes
print(stim_de.sort_values('logfoldchanges').head())  # Downregulated genes
print(stim_de.sort_values('logfoldchanges', ascending=False).head())  # Upregulated genes

# Remove tcells object to free memory
del tcells

# Load re-stimulated cells dataset ----
stim = sc.read_h5ad("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad")

# Calculate activation scores ----

# Using Log2FoldChange values as gene weights for the activation score
gene_weights = stim_de[stim_de['pvals_adj'] < 0.001].sort_values(by='logfoldchanges', ascending=False)
gene_weights = gene_weights[['names', 'logfoldchanges']].rename(columns={'names': 'mRNA_gene', 'logfoldchanges': 'weight'})

# Get average expression values for these genes in NT control cells
ctrl_avg_sct = sc.get.obs_df(stim, keys=gene_weights['mRNA_gene'], layer="SCT")
ctrl_avg_sct['mRNA_gene'] = gene_weights['mRNA_gene']
ctrl_avg_sct['NT.AVG'] = np.log1p(ctrl_avg_sct.mean(axis=1))  # Log transformation
ctrl_avg_sct = ctrl_avg_sct[['mRNA_gene', 'NT.AVG']]

# Calculate activation scores for each cell
# Extract SCT expression data for the relevant genes and compute scores
expr_data = pd.DataFrame(stim[:, gene_weights['mRNA_gene']].X.T, columns=stim.obs_names, index=gene_weights['mRNA_gene'])
expr_data = expr_data.reset_index().melt(id_vars='index', var_name='cell', value_name='value')
expr_data = expr_data.rename(columns={'index': 'mRNA_gene'})

# Merge gene weights and control average data
expr_data = expr_data.merge(gene_weights, on='mRNA_gene').merge(ctrl_avg_sct, on='mRNA_gene')
expr_data['score'] = (expr_data['value'] * expr_data['weight']) / expr_data['NT.AVG']

# Group by cell and sum to calculate total activation score for each cell
cell_scores = expr_data.groupby('cell').agg({'score': 'sum'}).reset_index().rename(columns={'score': 'activation.score'})

# Add cell scores to the object metadata
stim.obs = stim.obs.merge(cell_scores.set_index('cell'), left_index=True, right_index=True)

# Plotting violin plot of activation scores
sc.pl.violin(stim, keys='activation.score', groupby='crispr', show=False)
plt.show()

# Save the Seurat object (converted to AnnData format) ----
stim.write("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad", compression='gzip')

# UMAP plot ----
umap_df = pd.DataFrame(stim.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'])
umap_df['cell'] = stim.obs_names
umap_df = umap_df.merge(stim.obs[['activation.score']], left_index=True, right_index=True)

# Plot UMAP with activation scores
plt.figure(figsize=(10, 8))
sns.scatterplot(data=umap_df, x='UMAP_1', y='UMAP_2', hue='activation.score', palette='viridis', s=10)
plt.title("UMAP plot colored by activation scores")
plt.show()

# Boxplot for each sgRNA target vs NT control ----
activation_summary = stim.obs.groupby('crispr').agg({'activation.score': ['median', 'mean']}).reset_index()
activation_summary.columns = ['crispr', 'median_activation_score', 'mean_activation_score']

# Wilcoxon rank-sum test for each gene vs NT control
nt_scores = stim.obs.loc[stim.obs['crispr'] == 'NO-TARGET', 'activation.score']

wilcox_results = []
for gene in stim.obs['crispr'].unique():
    if gene != 'NO-TARGET':
        gene_scores = stim.obs.loc[stim.obs['crispr'] == gene, 'activation.score']
        res = wilcoxon(gene_scores, nt_scores)
        wilcox_results.append({'gene': gene, 'p_value': res.pvalue})

wilcox_df = pd.DataFrame(wilcox_results)
wilcox_df['stars'] = wilcox_df['p_value'].apply(lambda x: "***" if x < 0.001/69 else ("**" if x < 0.01/69 else ("*" if x < 0.05/69 else "")))

# Boxplot with significance stars
plt.figure(figsize=(10, 8))
sns.boxplot(x='crispr', y='activation.score', data=stim.obs, order=activation_summary['crispr'])
plt.xticks(rotation=90)
plt.show()