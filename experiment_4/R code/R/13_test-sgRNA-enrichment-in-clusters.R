# For each Seurat cluster, test for guide/target-gene overrespresentation

# For each cluster do the following:

## Generate contingency table of num cells for each gene (sgRNA target):

## e.g.
##                 gene.not.interest gene.of.interest
## In_cluster                   2613               28
## not_in_cluster              15310               29

## Use Fisher's Exact Test to test for overrepresentatino


# %%Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(broom)
library(qvalue)

# Import Seurat Objects ---------------------------------------------------
  
dat <- tibble(
  condition = c("Resting", "Re-stimulated"),
  object = list(
    LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Resting.h5Seurat"),
    LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5Seurat")
  )
)
dat


# Test for guide target overrepresentation in each cluster ----------------

dat <- dat %>%
  mutate(cluster_ids = map(object, function(x) levels(x@meta.data$seurat_clusters)))
dat$cluster_ids

dat <- dat %>%
  mutate(genes = map(object, function(x) x@meta.data$gene %>%
                       as.character() %>%
                       unique() %>%
                       list()))
dat$genes

dat


# Function to build contingency tables for a given gene/cluster combo 
get_contingency_table <- function(test_gene, cluster, data = tcells@meta.data) {
  
  # Initialize contingency table dataframe
  df <- data.frame(gene_test = c(0, 0), gene_other = c(0, 0))
  rownames(df) <- c("In_Cluster", "Not_In_Cluster")
  
  # Get total cells inside and outside cluster of interest
  cells_in_cluster <- data %>%
    dplyr::filter(seurat_clusters == cluster) %>%
    nrow()
  cells_outside_cluster <- data %>%
    dplyr::filter(seurat_clusters != cluster) %>%
    nrow()
  
  # Build Contingency table
  ## Num cells for test gene in cluster
  df$gene_test[1] <- data %>%
    dplyr::filter(gene == test_gene & seurat_clusters == cluster) %>%
    nrow()
  ## Num cells of other genes in cluster
  df$gene_other[1] <- cells_in_cluster - df$gene_test[1]
  ## Num cells for test gene outside cluster
  df$gene_test[2] <- data %>%
    dplyr::filter(gene == test_gene & seurat_clusters != cluster) %>%
    nrow()
  ## Num cells for other genes outside cluster
  df$gene_other[2] <- cells_outside_cluster - df$gene_test[2]
  
  return(df)
}

# Example: LAT2 in cluster 2 (stim)
get_contingency_table("LAT2", 2, dat$object[[2]]@meta.data)

# Dataframe for testing
dat <- dat %>%
  mutate(
    fisher_test_df = map2(cluster_ids, genes,
                          function(x,y) tibble(cluster_id = x, gene = y) %>%
                            unnest(gene)))

dat$fisher_test_df

# Apply fishers exact test
#   Resting
dat$fisher_test_df[[1]] <- dat$fisher_test_df[[1]] %>%
  mutate(test.res = map2(
    gene, cluster_id, 
    function(x,y) get_contingency_table(x,y, data = dat$object[[1]]@meta.data) %>%
      fisher.test() %>%
      broom::tidy()))
# Re-stim
dat$fisher_test_df[[2]] <- dat$fisher_test_df[[2]] %>%
  mutate(test.res = map2(
    gene, cluster_id, 
    function(x,y) get_contingency_table(x,y, data = dat$object[[2]]@meta.data) %>%
      fisher.test() %>%
      broom::tidy()))
dat$fisher_test_df

# Extract pvalue and odds ratios
dat$fisher_test_df[[1]] <- dat$fisher_test_df[[1]] %>%
  mutate(odds.ratio = map(test.res, function(x) pull(x, estimate)) %>% unlist(),
         p.value = map(test.res, function(x) pull(x, p.value)) %>% unlist()) 

dat$fisher_test_df[[2]] <- dat$fisher_test_df[[2]] %>%
  mutate(odds.ratio = map(test.res, function(x) pull(x, estimate)) %>% unlist(),
         p.value = map(test.res, function(x) pull(x, p.value)) %>% unlist()) 

dat$fisher_test_df

# FDR correction
#fisher_test_df$q.value <- qvalue(fisher_test_df$p.value)

#dat$fisher_test_df[[1]]$q.value <- qvalue(dat$fisher_test_df[[1]]$p.value)

dat <- dat %>%
  mutate(q.values = map(fisher_test_df, function(x) qvalue(x$p.value)),
         stats_df = map2(fisher_test_df, q.values, function(x,y) mutate(x, q.value = y$qvalues)))
dat

fisher_test_df <- dat %>%
  select(condition, stats_df) %>%
  unnest(stats_df) %>%
  select(-test.res)
fisher_test_df


# Plots -------------------------------------------------------------------

# Plot p.value distributions
fisher_test_df %>%
  ggplot(aes(p.value)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(~condition)


# Get log odds ratio for plotting purposes
fisher_test_df <- fisher_test_df %>%
  mutate(log_or = log2(odds.ratio+0.01))
fisher_test_df

# Cluster Enrichment Heatmap
crispr_enrichment_plot <- fisher_test_df %>%
  mutate(cluster_id = as_factor(cluster_id)) %>%
  # Cap log odds ratio and p.values
  mutate(log_or = if_else(log_or > 3, 3, log_or),
         log_or = if_else(log_or < -3, -3, log_or),
         q.value = if_else(q.value < 1e-100, 1e-100, q.value)) %>% 
  ggplot(aes(y = cluster_id, x = gene, fill = log_or, size = -log10(q.value))) +
  geom_point(shape = 21, color = "black") +
  # Only show points with -log10(p.value)>3 (p=0.001)
  scale_size_area(max_size = 4, limits = c(-log10(0.01),100)) +
  scale_fill_distiller(palette = "RdBu", limits = c(-3, 3)) +
  #scale_color_scico(palette = "berlin", limits = c(-3, 3)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  facet_grid(condition ~ .)
crispr_enrichment_plot



# Save table --------------------------------------------------------------

fisher_test_df %>%
  write_tsv("data/sgRNA-enrichment-in-clusters.txt")

