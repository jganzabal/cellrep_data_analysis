# For each cluster, test for guide/target-gene overrepresentation

# For each cluster do the following:

## Generate contingency table of num cells for each gene (sgRNA target):

## e.g.
##                 gene.not.interest gene.of.interest
## In_cluster                   2613               28
## not_in_cluster              15310               29

## Use Fisher's Exact Test to test for overrepresentation


# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Import AnnData Objects ---------------------------------------------------

dat = pd.DataFrame({
    'condition': ['Resting', 'Re-stimulated'],
    'object': [
        sc.read_h5ad("data/HuTcellsCRISPRaPerturbSeq_Resting.h5ad"),
        sc.read_h5ad("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad")
    ]
})
print(dat)


# Test for guide target overrepresentation in each cluster ----------------

dat['cluster_ids'] = dat['object'].apply(lambda x: x.obs['seurat_clusters'].cat.categories.tolist())
print(dat['cluster_ids'])

dat['genes'] = dat['object'].apply(lambda x: x.obs['gene'].unique().tolist())
print(dat['genes'])

print(dat)


# Function to build contingency tables for a given gene/cluster combo 
def get_contingency_table(test_gene, cluster, data):
    
    # Initialize contingency table dataframe
    df = pd.DataFrame({'gene_test': [0, 0], 'gene_other': [0, 0]},
                      index=['In_Cluster', 'Not_In_Cluster'])
    
    # Get total cells inside and outside cluster of interest
    cells_in_cluster = data[data['seurat_clusters'] == cluster].shape[0]
    cells_outside_cluster = data[data['seurat_clusters'] != cluster].shape[0]
    
    # Build Contingency table
    ## Num cells for test gene in cluster
    df.loc['In_Cluster', 'gene_test'] = data[(data['gene'] == test_gene) & (data['seurat_clusters'] == cluster)].shape[0]
    ## Num cells of other genes in cluster
    df.loc['In_Cluster', 'gene_other'] = cells_in_cluster - df.loc['In_Cluster', 'gene_test']
    ## Num cells for test gene outside cluster
    df.loc['Not_In_Cluster', 'gene_test'] = data[(data['gene'] == test_gene) & (data['seurat_clusters'] != cluster)].shape[0]
    ## Num cells for other genes outside cluster
    df.loc['Not_In_Cluster', 'gene_other'] = cells_outside_cluster - df.loc['Not_In_Cluster', 'gene_test']
    
    return df

# Example: LAT2 in cluster 2 (stim)
print(get_contingency_table("LAT2", "2", dat['object'].iloc[1].obs))

# Dataframe for testing
dat['fisher_test_df'] = dat.apply(lambda row: pd.DataFrame({
    'cluster_id': row['cluster_ids'],
    'gene': row['genes']
}).explode('gene'), axis=1)

print(dat['fisher_test_df'])

# Apply fishers exact test
def apply_fisher_test(row, data):
    contingency_table = get_contingency_table(row['gene'], row['cluster_id'], data)
    odds_ratio, p_value = stats.fisher_exact(contingency_table)
    return pd.Series({'odds_ratio': odds_ratio, 'p_value': p_value})

for i, row in dat.iterrows():
    dat.loc[i, 'fisher_test_df'] = dat.loc[i, 'fisher_test_df'].join(
        dat.loc[i, 'fisher_test_df'].apply(apply_fisher_test, data=row['object'].obs, axis=1)
    )

print(dat['fisher_test_df'])

# FDR correction
for i, row in dat.iterrows():
    _, q_values, _, _ = multipletests(row['fisher_test_df']['p_value'], method='fdr_bh')
    dat.loc[i, 'fisher_test_df']['q_value'] = q_values

fisher_test_df = dat[['condition', 'fisher_test_df']].explode('fisher_test_df').reset_index(drop=True)
fisher_test_df = pd.concat([fisher_test_df['condition'], fisher_test_df['fisher_test_df'].apply(pd.Series)], axis=1)
print(fisher_test_df)


# Plots -------------------------------------------------------------------

# Plot p.value distributions
g = sns.FacetGrid(fisher_test_df, col="condition")
g.map(plt.hist, "p_value", bins=20)
plt.show()

# Get log odds ratio for plotting purposes
fisher_test_df['log_or'] = np.log2(fisher_test_df['odds_ratio'] + 0.01)

# Cluster Enrichment Heatmap
fisher_test_df['log_or'] = fisher_test_df['log_or'].clip(-3, 3)
fisher_test_df['q_value'] = fisher_test_df['q_value'].clip(1e-100, None)

g = sns.FacetGrid(fisher_test_df, col="condition", height=10, aspect=1.5)
g.map(sns.scatterplot, "gene", "cluster_id", "log_or", "-np.log10(q_value)", 
      palette="RdBu", hue_norm=(-3, 3), size_norm=(0, 100))
g.add_legend(title="log2(Odds Ratio)", label_order=sorted(fisher_test_df['log_or'].unique()))
g.set_xticklabels(rotation=45, ha='right')
plt.tight_layout()
plt.show()

# Save table --------------------------------------------------------------

fisher_test_df.to_csv("data/sgRNA-enrichment-in-clusters.txt", sep="\t", index=False)

