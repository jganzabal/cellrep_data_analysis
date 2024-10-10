# Idea is to assign an unbiased T cell activation score to each cell, based on 
# the expression of differentially expressed genes upon T cell activation -->
# Will use NO-TARGET control cells to determine these differentially expressed
# genes

# %%Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(broom)


# Differential expression of No-Target control cells stim vs nostim -------

# Import Seurat object with both restim and resting conditions
tcells <- LoadH5Seurat("data_out/tcells_all_filt_06.h5Seurat")

# Differential expression of No-Target control cells stim vs nostim -------
# Keep No target
ctrl_cells <- tcells %>%
  subset(crispr == "NT")

ctrl_cells <- ctrl_cells %>%
  RegroupIdents("condition")

# Does Differential expression Analysis
# It uses features to do the differential expression Analysis
# i.e. gene columns
stim_de <- ctrl_cells %>%
  FindMarkers(
    ident.1 = "Stim",
    ident.2 = "Nostim",
    min.pct = 0.25  # only test genes that are detected in a minimum
  )                 # fraction of min.pct cells in either of the two populations

# Show less expressed genes
stim_de %>%
  arrange(avg_log2FC) %>%
  head()

# Show most expressed genes
stim_de %>%
  arrange(desc(avg_log2FC)) %>%
  head()

# Remove tcells object to free up memory
rm(tcells)


# Import Re-stim Only object -----------------------------------------------

stim <- LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5Seurat")



# Get activation scores ---------------------------------------------------

# Use Log2FoldChange values for gene weights

gene_weights <- stim_de %>%
  dplyr::filter(p_val_adj < 0.001) %>%    # Relevant genes
  rownames_to_column("mRNA_gene") %>%     # Names de index Messenger RNA
  tibble() %>%                            # Convert to dataframe
  arrange(desc(avg_log2FC)) %>%           # Order on columnavg_log2FC
  select(mRNA_gene, weight = avg_log2FC)  # Keep mRNA_gene and avg_log2FC
gene_weights

# Get average expression values for these genes in NT control cells
ctrl_avg_sct <- AverageExpression(
  stim,
  features = gene_weights$mRNA_gene,
  group.by = "crispr", # SCT.perturbed, SCT.NT
  assays = "SCT"
) %>%
  as.data.frame() %>%
  rownames_to_column("mRNA_gene") %>%
  tibble() %>%
  select(mRNA_gene, NT.AVG = SCT.NT) %>%
  mutate(NT.AVG = log(NT.AVG + 1)) # Re-log data
ctrl_avg_sct

# Calculate activation scores for eachc cell
cell_scores <- stim@assays$SCT[gene_weights$mRNA_gene,] %>%
  as.data.frame() %>%
  rownames_to_column("mRNA_gene") %>%
  tibble() %>%
  gather(cell, value, -mRNA_gene) %>%
  full_join(gene_weights, by = "mRNA_gene") %>%
  full_join(ctrl_avg_sct, by = "mRNA_gene") %>%
  mutate(score = (value * weight) / NT.AVG) %>%
  group_by(cell) %>%
  summarise(activation.score = sum(score))

cell_scores

# Add cell scores to object metadata
cell_scores <- cell_scores %>%
  as.data.frame() %>%
  column_to_rownames("cell")
stim <- stim %>%
  AddMetaData(cell_scores)

head(stim)

VlnPlot(stim, "activation.score", pt.size = 0) +
  theme(legend.position = "NULL")


# Save Seurat Object ------------------------------------------------------

SaveH5Seurat(stim, 
             "data_out/HuTcellsCRISPRaPerturbSeq_Re-stimulated_10.h5Seurat",
             overwrite = TRUE)


# UMAP plot ---------------------------------------------------------------

umap_df_stim <- stim@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  tibble()

umap_df_stim <- full_join(umap_df_stim,
                          stim@meta.data %>%
                            as.data.frame() %>%
                            rownames_to_column("cell"),
                          by = "cell")
umap_df_stim

umap_activation.score <- umap_df_stim %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, z = activation.score)) +
  stat_summary_hex(bins = 50) +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom") 
umap_activation.score


# Boxplot -----------------------------------------------------------------

# Stats - Each sgRNA target vs NT control
# Cell metadata
metadata <- stim@meta.data %>%
  rownames_to_column("cell") %>%
  tibble()

activation_summary <- metadata %>%
  group_by(gene) %>%
  summarise(median.activation.score = median(activation.score),
            mean.activation.score = mean(activation.score)) %>%
  arrange(desc(median.activation.score))
activation_summary

# Crispr genes
# The test will be done with rows, not columns
# Can't use FindMarkers?
wilcox_df <- metadata %>%
  select(gene, activation.score) %>%
  dplyr::filter(gene != "NO-TARGET") %>%
  group_by(gene) %>%
  nest() %>%
  mutate(data = map(data, pull))
wilcox_df

NT_activation.scores <- metadata %>%
  dplyr::filter(gene == "NO-TARGET") %>%
  pull(activation.score)
NT_activation.scores

wilcox_df <- wilcox_df %>%
  mutate(wilcox.res = map(data, wilcox.test, y = NT_activation.scores))

wilcox_df <- wilcox_df %>%
  mutate(wilcox.tidy = map(wilcox.res, broom::tidy)) %>%
  select(gene, wilcox.tidy) %>%
  unnest(wilcox.tidy)

wilcox_df


# Add significance stars with Bonferroni correction
wilcox_df <- wilcox_df %>%
  mutate(stars = if_else(p.value < (0.05/69), "*", ""),
         stars = if_else(p.value < (0.01/69), "**", stars),
         stars = if_else(p.value < (0.001/69), "***", stars))
wilcox_df

metadata %>%
  mutate(gene = factor(gene, levels = rev(activation_summary$gene))) %>%
  mutate(color = if_else(gene_category == "Negative_Both", "Negative Regulator", "Positive Regulator"),
         color = if_else(gene %in% c("NO-TARGET", "IFNG", "IL2"), "Control", color),
         color = factor(color, levels = c("Control", "Negative Regulator", "Positive Regulator"))) %>%
  ggplot(aes(x = gene, y = activation.score)) +
  geom_boxplot(aes(fill = color), outlier.shape = NA) +
  geom_hline(yintercept = activation_summary %>%
               dplyr::filter(gene == "NO-TARGET") %>%
               pull(median.activation.score),
             linetype = "dashed") +
  geom_text(aes(y = 605, label = stars),
            data = wilcox_df,
            size = 4.5,
            lineheight = 0.25,
            vjust = 0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  coord_flip(ylim = c(-300, 600)) +
  scale_fill_manual(values = c("Grey", "#67a9cf", "#ef8a62")) +
  xlab("CRISPRa Target Gene") +
  ylab("Activation Score") 
