# Generate 2D contour plots of sgRNA distributions in UMAP space

# Load Packages -----------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Optional
import anndata

# Import AnnData Objects ---------------------------------------------------

dat = pd.DataFrame({
    'condition': ['Resting', 'Re-stimulated'],
    'object': [
        sc.read_h5ad("data/HuTcellsCRISPRaPerturbSeq_Resting.h5ad"),
        sc.read_h5ad("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad")
    ]
})
print(dat)

# Re-stim Contour plots ----------------------------------------------------

n_bins = 12
plot_limits_x = (-8.5, 9.5)
plot_limits_y = (-5.5, 6.5)

# Contour plot of all cells (re-stimulated)
umap_coords = dat['object'].iloc[1].obsm['X_umap']
plt.figure(figsize=(10, 8))
sns.kdeplot(x=umap_coords[:, 0], y=umap_coords[:, 1], cmap="viridis", shade=True, cbar=True, n_levels=n_bins)
plt.xlim(plot_limits_x)
plt.ylim(plot_limits_y)
plt.title("All Cells (Re-stimulated)")
plt.show()

# Function to generate a contour plot for a given gene or category
def make_contour_plot(adata: anndata.AnnData, 
                      Gene: Optional[str] = None, 
                      category: Optional[str] = None, 
                      control_backdrop: Optional[plt.Figure] = None):
    
    # Subset data by gene or gene category
    if Gene is not None:
        dat = adata[adata.obs['gene'] == Gene]
    elif category is not None:
        dat = adata[adata.obs['gene_category'] == category]
    else:
        dat = adata
    
    umap_coords = dat.obsm['X_umap']
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Build plot w control layer
    if control_backdrop is not None:
        ax = control_backdrop.gca()
        sns.kdeplot(x=umap_coords[:, 0], y=umap_coords[:, 1], cmap="viridis", ax=ax, n_levels=n_bins)
    else:
        sns.kdeplot(x=umap_coords[:, 0], y=umap_coords[:, 1], cmap="viridis", shade=True, cbar=True, ax=ax, n_levels=n_bins)
    
    ax.set_xlim(plot_limits_x)
    ax.set_ylim(plot_limits_y)
    ax.set_title(f"{Gene or category}")
    
    return fig

# Generate contour of only NT control cells for backdrop comparison
control_adata = dat['object'].iloc[1][dat['object'].iloc[1].obs['gene'] == "NO-TARGET"]
control_umap = control_adata.obsm['X_umap']

control_fig, control_ax = plt.subplots(figsize=(10, 8))
sns.kdeplot(x=control_umap[:, 0], y=control_umap[:, 1], cmap="Greys", shade=True, cbar=True, ax=control_ax, n_levels=n_bins)
control_ax.set_xlim(plot_limits_x)
control_ax.set_ylim(plot_limits_y)
control_ax.set_title("Non-Targeting Controls")
plt.show()

# Test function
make_contour_plot(dat['object'].iloc[1], Gene="MAP4K1")
make_contour_plot(dat['object'].iloc[1], Gene="MAP4K1", control_backdrop=control_fig)

# Apply function to all sgRNA targets
gene_contour_plots = []
for gene in dat['object'].iloc[1].obs['gene'].unique():
    gene_contour_plots.append(make_contour_plot(dat['object'].iloc[1], Gene=gene, control_backdrop=control_fig))

# Plots for manuscript main figure

genes_to_plot = ["VAV1", "MAP4K1", "FOXQ1", "GATA3", "IL1R1", "TNFRSF1A", "TBX21"]

fig, axes = plt.subplots(3, 3, figsize=(20, 20))
axes = axes.flatten()

# All Cells
sns.kdeplot(x=umap_coords[:, 0], y=umap_coords[:, 1], cmap="viridis", shade=True, ax=axes[0], n_levels=n_bins)
axes[0].set_title("All Cells")

# Non-Targeting Controls
sns.kdeplot(x=control_umap[:, 0], y=control_umap[:, 1], cmap="Greys", shade=True, ax=axes[1], n_levels=n_bins)
axes[1].set_title("Non-Targeting Controls")

# Perturbed Cells
perturbed_adata = dat['object'].iloc[1][dat['object'].iloc[1].obs['gene'] != "NO-TARGET"]
perturbed_umap = perturbed_adata.obsm['X_umap']
sns.kdeplot(x=perturbed_umap[:, 0], y=perturbed_umap[:, 1], cmap="viridis", shade=True, ax=axes[2], n_levels=n_bins)
axes[2].set_title("Perturbed Cells")

# Individual genes
for i, gene in enumerate(genes_to_plot):
    gene_adata = dat['object'].iloc[1][dat['object'].iloc[1].obs['gene'] == gene]
    gene_umap = gene_adata.obsm['X_umap']
    sns.kdeplot(x=gene_umap[:, 0], y=gene_umap[:, 1], cmap="viridis", shade=True, ax=axes[i+3], n_levels=n_bins)
    axes[i+3].set_title(gene)

for ax in axes:
    ax.set_xlim(plot_limits_x)
    ax.set_ylim(plot_limits_y)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])

plt.tight_layout()
plt.show()

# Resting cells -----------------------------------------------------------

# Contour plots -----------------------------------------------------------

n_bins = 12
plot_limits_x = (-7.3, 10.5)
plot_limits_y = (-7.5, 7.8)

# Contour plot of all cells
umap_coords = dat['object'].iloc[0].obsm['X_umap']
plt.figure(figsize=(10, 8))
sns.kdeplot(x=umap_coords[:, 0], y=umap_coords[:, 1], cmap="viridis", shade=True, cbar=True, n_levels=n_bins)
plt.xlim(plot_limits_x)
plt.ylim(plot_limits_y)
plt.title("All Cells (Resting)")
plt.show()

# Backdrop control contour
control_adata = dat['object'].iloc[0][dat['object'].iloc[0].obs['gene'] == "NO-TARGET"]
control_umap = control_adata.obsm['X_umap']

control_fig, control_ax = plt.subplots(figsize=(10, 8))
sns.kdeplot(x=control_umap[:, 0], y=control_umap[:, 1], cmap="Greys", shade=True, cbar=True, ax=control_ax, n_levels=n_bins)
control_ax.set_xlim(plot_limits_x)
control_ax.set_ylim(plot_limits_y)
control_ax.set_title("Non-Targeting Controls")
plt.show()

# Apply function to guide genes
gene_contour_plots = []
for gene in dat['object'].iloc[0].obs['gene'].unique():
    gene_contour_plots.append(make_contour_plot(dat['object'].iloc[0], Gene=gene, control_backdrop=control_fig))

# Plots for manuscript main figure

fig, axes = plt.subplots(3, 3, figsize=(20, 20))
axes = axes.flatten()

# All Cells
sns.kdeplot(x=umap_coords[:, 0], y=umap_coords[:, 1], cmap="viridis", shade=True, ax=axes[0], n_levels=n_bins)
axes[0].set_title("All Cells")

# Non-Targeting Controls
sns.kdeplot(x=control_umap[:, 0], y=control_umap[:, 1], cmap="Greys", shade=True, ax=axes[1], n_levels=n_bins)
axes[1].set_title("Non-Targeting Controls")

# Perturbed Cells
perturbed_adata = dat['object'].iloc[0][dat['object'].iloc[0].obs['gene'] != "NO-TARGET"]
perturbed_umap = perturbed_adata.obsm['X_umap']
sns.kdeplot(x=perturbed_umap[:, 0], y=perturbed_umap[:, 1], cmap="viridis", shade=True, ax=axes[2], n_levels=n_bins)
axes[2].set_title("Perturbed Cells")

# Individual genes
for i, gene in enumerate(genes_to_plot):
    gene_adata = dat['object'].iloc[0][dat['object'].iloc[0].obs['gene'] == gene]
    gene_umap = gene_adata.obsm['X_umap']
    sns.kdeplot(x=gene_umap[:, 0], y=gene_umap[:, 1], cmap="viridis", shade=True, ax=axes[i+3], n_levels=n_bins)
    axes[i+3].set_title(gene)

for ax in axes:
    ax.set_xlim(plot_limits_x)
    ax.set_ylim(plot_limits_y)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])

plt.tight_layout()
plt.show()
