# Pseudobulk differential expression analysis of sgRNA targets versus NT control
# cells for each cytokine (or other effector) gene - Resting condition

# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Import AnnData Object ---------------------------------------------------

resting = sc.read_h5ad("data/HuTcellsCRISPRaPerturbSeq_Resting.h5ad")

# Gene Lists --------------------------------------------------------------

cytokine_genes = pd.read_csv("data/GO_0005125_Cytokines.txt", sep="\t", header=None)
cytokine_genes = sorted(cytokine_genes[2].unique())

granzyme_genes = ["GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "THBS1"]
# Why TSP1: https://science.sciencemag.org/content/368/6493/897

detected_features = resting.var_names

cytokines_in_data = [gene for gene in cytokine_genes if gene in detected_features]
print(cytokines_in_data)
print([gene for gene in cytokine_genes if gene not in detected_features])

granzyme_genes_in_data = [gene for gene in granzyme_genes if gene in detected_features]
print(granzyme_genes_in_data)

test_genes = cytokines_in_data + granzyme_genes_in_data

# Find Differentially Expressed Genes -------------------------------------

# Versus non-targeting control cells, genes in library set

# Get sgRNA target genes from library
targetgenes = resting.obs[resting.obs['gene'] != "NO-TARGET"]['gene'].unique()
print(targetgenes)

# Function to get differential expression for target gene of interest versus
# non-targeting control cells (in a given condition)
def target_gene_overexpression(adata, Gene, GeneSet=test_genes):
    adata_subset = adata[adata.obs['gene'].isin(["NO-TARGET", Gene])]
    sc.tl.rank_genes_groups(adata_subset, 'crispr', groups=['perturbed'], reference='NT', method='wilcoxon')
    markers = sc.get.rank_genes_groups_df(adata_subset, group='perturbed')
    markers = markers[markers['names'].isin(GeneSet)]
    markers = markers.rename(columns={'names': 'gene', 'scores': 'avg_log2FC', 'pvals_adj': 'p_val_adj'})
    return markers[['gene', 'avg_log2FC', 'p_val_adj']]

# Dataframe for differential expression testing
de_df = pd.DataFrame({'gene': targetgenes})

# CD8 DE 
CD8 = resting[resting.obs['CD4.or.CD8'] == "CD8"]

de_resting_CD8 = de_df.copy()
de_resting_CD8['DiffExpr'] = de_resting_CD8['gene'].apply(lambda x: target_gene_overexpression(CD8, x))
de_resting_CD8 = de_resting_CD8.explode('DiffExpr').reset_index(drop=True)
de_resting_CD8 = pd.concat([de_resting_CD8.drop('DiffExpr', axis=1), pd.DataFrame(de_resting_CD8['DiffExpr'].tolist())], axis=1)
de_resting_CD8 = de_resting_CD8.rename(columns={'gene': 'crispr_gene', 'gene': 'mRNA_gene'})

# CD4 DE 
CD4 = resting[resting.obs['CD4.or.CD8'] == "CD4"]

de_resting_CD4 = de_df.copy()
de_resting_CD4['DiffExpr'] = de_resting_CD4['gene'].apply(lambda x: target_gene_overexpression(CD4, x))
de_resting_CD4 = de_resting_CD4.explode('DiffExpr').reset_index(drop=True)
de_resting_CD4 = pd.concat([de_resting_CD4.drop('DiffExpr', axis=1), pd.DataFrame(de_resting_CD4['DiffExpr'].tolist())], axis=1)
de_resting_CD4 = de_resting_CD4.rename(columns={'gene': 'crispr_gene', 'gene': 'mRNA_gene'})

# Get Cytokine / CRISPRa categories ---------------------------------------

# Add guide and cytokine categories 

guide_categories = pd.read_csv("data/guide-target-categories.txt", sep="\t")

# Switch to factors and rearrange levels
guide_categories['gene_functional_category'] = pd.Categorical(
    guide_categories['gene_functional_category'],
    categories=["Negative Regulator", "Control", "TCR Stim/Co-Stim", 
                "GTPase signaling", "NF-kB Signaling", "TNF receptor superfamily", 
                "Other Cytokine Receptor", "Other Signal Transduction",
                "Transcription Factor", "Other / Unknown"],
    ordered=True
)

# Add some line-breaks to categories for more compact plotting
guide_categories['gene_functional_category'] = guide_categories['gene_functional_category'].cat.rename_categories({
    "Negative Regulator": "Negative\nRegulator",
    "TCR Stim/Co-Stim": "TCR\nStim/Co-Stim",
    "GTPase signaling": "GTPase\nsignaling",
    "NF-kB Signaling": "NF-kB\nSignaling",
    "TNF receptor superfamily": "TNF\nreceptor\nsuperfamily",
    "Other Cytokine Receptor": "Other\nCytokine\nReceptor",
    "Other Signal Transduction": "Other\nSignal\nTransduction",
    "Transcription Factor": "Transcription\nFactor",
    "Other / Unknown": "Other /\nUnknown"
})

# Do the same for cytokine categories
cytokine_categories = pd.read_csv("data/cytokine_categories.txt", sep="\t")
cytokine_categories = cytokine_categories.dropna(subset=['Category'])

cytokine_categories['Category'] = pd.Categorical(
    cytokine_categories['Category'],
    categories=["Th1","Th2","Chemokine","Growth Factor","Apoptosis","Other",
                "Granzyme/Cytolytic"],
    ordered=True
)

cytokine_categories['Category'] = cytokine_categories['Category'].cat.rename_categories({
    "Growth Factor": "Growth\nFactor",
    "Granzyme/Cytolytic": "Granzyme/\nCytolytic"
})

# Plots -------------------------------------------------------------------

max_lfc = 1
min_lfc = -1
min_padj = 1e-10

# Add categories to dataframes for plotting
CD4_plot_df = pd.merge(cytokine_categories.rename(columns={'Cytokine': 'mRNA_gene'}),
                       de_resting_CD4,
                       on='mRNA_gene')
CD4_plot_df = pd.merge(CD4_plot_df,
                       guide_categories.rename(columns={'gene': 'crispr_gene'}),
                       on='crispr_gene')
CD4_plot_df = CD4_plot_df.rename(columns={'Category': 'cytokine_category'})

CD8_plot_df = pd.merge(cytokine_categories.rename(columns={'Cytokine': 'mRNA_gene'}),
                       de_resting_CD8,
                       on='mRNA_gene')
CD8_plot_df = pd.merge(CD8_plot_df,
                       guide_categories.rename(columns={'gene': 'crispr_gene'}),
                       on='crispr_gene')
CD8_plot_df = CD8_plot_df.rename(columns={'Category': 'cytokine_category'})

# Filter out targets with no significant DE expressed genes 
genes2keep_CD4CD8 = pd.concat([CD4_plot_df, CD8_plot_df])
genes2keep_CD4CD8 = genes2keep_CD4CD8[genes2keep_CD4CD8['p_val_adj'] < 0.05]['crispr_gene'].unique()

# Plots
def create_plot(plot_df, title):
    plot_df = plot_df[plot_df['crispr_gene'].isin(genes2keep_CD4CD8)]
    plot_df['avg_log2FC'] = plot_df['avg_log2FC'].clip(lower=min_lfc, upper=max_lfc)
    plot_df['p_val_adj'] = plot_df['p_val_adj'].clip(lower=min_padj)
    plot_df['significant'] = plot_df['p_val_adj'] < 0.05

    g = sns.FacetGrid(plot_df, col="cytokine_category", row="gene_functional_category", 
                      margin_titles=True, sharex=False, sharey=False, height=3, aspect=1)
    g.map(sns.scatterplot, "crispr_gene", "mRNA_gene", "avg_log2FC", "-log10(p_val_adj)", 
          edgecolor="black", linewidth=0.5)
    g.add_legend(title="avg_log2FC", label_order=sorted(plot_df['avg_log2FC'].unique()))
    g.fig.suptitle(title)
    g.set_axis_labels("CRISPRa Target", "Cytokine Gene")
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    plt.tight_layout()
    return g

CD4_plot = create_plot(CD4_plot_df, "CD4+ T cells")
CD8_plot = create_plot(CD8_plot_df, "CD8+ T cells")

plt.show()

# plt.savefig("cytokines-resting-DE.pdf", bbox_inches="tight")
