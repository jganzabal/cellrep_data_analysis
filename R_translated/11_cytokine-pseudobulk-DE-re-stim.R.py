# Pseudobulk differential expression analysis of sgRNA targets versus NT control
# cells for each cytokine (or other effector) gene - Re-stimulated condition

# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData
from typing import List, Dict
import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

# Import AnnData Object ---------------------------------------------------

stim = sc.read_h5ad("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad")

# Gene Lists --------------------------------------------------------------

cytokine_genes = pd.read_csv("data/GO_0005125_Cytokines.txt", sep="\t", header=None)[2].unique().sort_values().tolist()
print(cytokine_genes)

granzyme_genes = ["GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "THBS1"]
# Why TSP1: https://science.sciencemag.org/content/368/6493/897

detected_features = stim.var_names.tolist()

cytokines_in_data = [gene for gene in cytokine_genes if gene in detected_features]
print(cytokines_in_data)
print([gene for gene in cytokine_genes if gene not in detected_features])

granzyme_genes_in_data = [gene for gene in granzyme_genes if gene in detected_features]
print(granzyme_genes_in_data)

test_genes = cytokines_in_data + granzyme_genes_in_data

# Find Differentially Expressed Genes -------------------------------------

# Versus non-targeting control cells, genes in library set

# Get sgRNA target genes from library
targetgenes = stim.obs[stim.obs['gene'] != "NO-TARGET"]['gene'].unique().tolist()
print(targetgenes)

# Function to get differential expression for target gene of interest versus
# non-targeting control cells (in a given condition)
def target_gene_overexpression(adata: AnnData, 
                               Gene: str, 
                               GeneSet: List[str] = test_genes) -> pd.DataFrame:
    adata_subset = adata[adata.obs['gene'].isin(["NO-TARGET", Gene])]
    adata_subset.obs['crispr'] = adata_subset.obs['gene'].map(lambda x: "perturbed" if x != "NO-TARGET" else "NT")
    sc.tl.rank_genes_groups(adata_subset, 'crispr', groups=['perturbed'], reference='NT', method='wilcoxon', gene_symbols='index')
    de_results = sc.get.rank_genes_groups_df(adata_subset, group='perturbed')
    de_results = de_results[de_results['names'].isin(GeneSet)]
    return de_results

# Dataframe for differential expression testing
de_df = pd.DataFrame({'gene': targetgenes})
print(de_df)

# Apply function
de_stim_all = []
for gene in de_df['gene']:
    de_result = target_gene_overexpression(stim, gene)
    de_result['crispr_gene'] = gene
    de_stim_all.append(de_result)

de_stim_all = pd.concat(de_stim_all)
de_stim_all = de_stim_all.rename(columns={'names': 'mRNA_gene'})
print(de_stim_all)

# Save table --------------------------------------------------------------

#de_stim_all.to_csv("DE_NT-vs-sgRNA_cytokine-genes_re-stim.txt", sep="\t", index=False)

# Heatmap plot ------------------------------------------------------------

# Add guide and cytokine categories 

guide_categories = pd.read_csv("data/guide-target-categories.txt", sep="\t")

# Switch to factors and rearrange levels
category_order = ["Negative Regulator", "Control", "TCR Stim/Co-Stim", 
                  "GTPase signaling", "NF-kB Signaling", "TNF receptor superfamily", 
                  "Other Cytokine Receptor", "Other Signal Transduction",
                  "Transcription Factor", "Other / Unknown"]
guide_categories['gene_functional_category'] = pd.Categorical(
    guide_categories['gene_functional_category'],
    categories=category_order,
    ordered=True
)

# Add some line-breaks to categories for more compact plotting
category_recode = {
    "Negative Regulator": "Negative\nRegulator",
    "TCR Stim/Co-Stim": "TCR\nStim/Co-Stim",
    "GTPase signaling": "GTPase\nsignaling",
    "NF-kB Signaling": "NF-kB\nSignaling",
    "TNF receptor superfamily": "TNF\nreceptor\nsuperfamily",
    "Other Cytokine Receptor": "Other\nCytokine\nReceptor",
    "Other Signal Transduction": "Other\nSignal\nTransduction",
    "Transcription Factor": "Transcription\nFactor",
    "Other / Unknown": "Other /\nUnknown"
}
guide_categories['gene_functional_category'] = guide_categories['gene_functional_category'].cat.rename_categories(category_recode)

# Do the same for cytokine categories
cytokine_categories = pd.read_csv("data/cytokine_categories.txt", sep="\t")
cytokine_categories = cytokine_categories.dropna(subset=['Category'])

category_order = ["Th1", "Th2", "Chemokine", "Growth Factor", "Apoptosis", "Other", "Granzyme/Cytolytic"]
cytokine_categories['Category'] = pd.Categorical(
    cytokine_categories['Category'],
    categories=category_order,
    ordered=True
)

category_recode = {
    "Growth Factor": "Growth\nFactor",
    "Granzyme/Cytolytic": "Granzyme/\nCytolytic"
}
cytokine_categories['Category'] = cytokine_categories['Category'].cat.rename_categories(category_recode)

# Join with data
plot_df = cytokine_categories.rename(columns={'Cytokine': 'mRNA_gene'}).merge(
    de_stim_all, on='mRNA_gene'
).merge(
    guide_categories.rename(columns={'gene': 'crispr_gene'}), on='crispr_gene'
).rename(columns={'Category': 'cytokine_category'})

# Remove crispra target genes that do not have any significant cytokine genes
genes2keep = plot_df[plot_df['pvals_adj'] < 0.05]['crispr_gene'].unique()

# Plot
max_lfc = 1
min_lfc = -1
min_padj = 1e-10

plot_data = plot_df[plot_df['crispr_gene'].isin(genes2keep)].copy()
plot_data['logfoldchanges'] = plot_data['logfoldchanges'].clip(min_lfc, max_lfc)
plot_data['pvals_adj'] = plot_data['pvals_adj'].clip(min_padj)
plot_data['significant'] = plot_data['pvals_adj'] < 0.05

g = sns.FacetGrid(plot_data, 
                  row='gene_functional_category', 
                  col='cytokine_category', 
                  margin_titles=True, 
                  sharex=False, 
                  sharey=False, 
                  height=4, 
                  aspect=1)

g.map(sns.scatterplot, 'crispr_gene', 'mRNA_gene', 
      'logfoldchanges', '-log10(pvals_adj)', 
      edgecolor='black', linewidth=0.5)

g.add_legend(title='Log2 Fold Change', bbox_to_anchor=(1.05, 1), loc='upper left')
g.set_axis_labels('CRISPRa Target', 'Cytokine Gene')
g.set_titles(col_template='{col_name}', row_template='{row_name}')

for ax in g.axes.flat:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
plt.tight_layout()
plt.show()

# The rest of the code (CD4/CD8 T cells analysis) follows a similar pattern.
# You would need to subset the data for CD4 and CD8 T cells, perform differential expression,
# and create similar plots as above for each subset.

