# Generate a heatmap of top markers in each cluster, with top markers, guides,
# cytokine genes, and cluster names indicated to the right (Fig. 4H)

# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)
library(scico)

# %% Import Data -------------------------------------------------------------

dat <- tibble(
  condition = c("Resting", "Re-stimulated"),
  object = list(
    LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Resting.h5Seurat"),
    LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5Seurat")
  )
)
dat

marker_df <- read_tsv("data/clusterMarkers.txt")
marker_df

guide_df <- read_tsv("data/sgRNA-enrichment-in-clusters.txt")
guide_df


# Get genes to plot for each cluster --------------------------------------

# Function to pull the top genes N for each cluster, such that each gene only
# appears once. If a gene is a signficantly enriched marker more than once, the
# cluster with the greater logFoldChange enrichment is prioritized

get_unique_cluster_markers <- function(df,
                                       top_genes = 50,
                                       p_adj_filt = 0.25) {
  df <- df %>%
    dplyr::filter(p_val_adj < p_adj_filt & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  clusters <- c(min(df$cluster):max(df$cluster))
  
  topGenesList <- vector("list", length(clusters))
  names(topGenesList) <- clusters
  
  for (i in 1:nrow(df)) {
    gene <- df$gene[i]
    cluster <- as.character(df$cluster[i]) 
    genes_used <- reduce(topGenesList, c)
    if(!gene %in% genes_used) {
      if(length(topGenesList[[cluster]]) < top_genes) {
        topGenesList[[cluster]] <- c(topGenesList[[cluster]], gene)
      }
    }
  }
  return(topGenesList)
}

markers <- marker_df %>%
  group_by(condition) %>%
  nest() %>%
  mutate(genes_for_heatmap = map(data, get_unique_cluster_markers))
markers

markers <- markers %>%
  mutate(genes_for_heatmap = map(genes_for_heatmap, 
                                 stack),
         genes_for_heatmap = map(genes_for_heatmap, rename, 
                                 gene = values, 
                                 cluster = ind))
markers
markers$genes_for_heatmap[[1]] %>% head()
markers$genes_for_heatmap[[2]] %>% head()


# # Re-level clusters for plotting in correct order
# markers$genes_for_heatmap[[1]]$cluster <- factor(markers$genes_for_heatmap[[1]]$cluster,
#                                                  levels = levels(dat$object[[1]]))
# 
# markers$genes_for_heatmap[[2]]$cluster <- factor(markers$genes_for_heatmap[[2]]$cluster,
#                                                  levels = levels(dat$object[[2]]))
# 
# # Arrange data by cluster levels
# markers$genes_for_heatmap[[1]] <- markers$genes_for_heatmap[[1]] %>%
#   arrange(cluster)
# markers$genes_for_heatmap[[2]] <- markers$genes_for_heatmap[[2]] %>%
#   arrange(cluster)
# 
# markers

# Join with dat dataframe
dat <- left_join(dat,
                 rename(markers, 
                        marker_df = data, 
                        marker_cluster_df = genes_for_heatmap),
                 by = "condition")
dat

rm(markers) # Clean env


# Get average expression for these genes ----------------------------------

dat <- dat %>%
  mutate(
    average_expression = map2(
      object, marker_cluster_df,
      function(x,y) AverageExpression(x, 
                                      features = y$gene,
                                      assays = "SCT")
    )
  )
dat

dat <- dat %>%
  mutate(average_expression = map(average_expression,
                                  function(x) x$SCT %>%
                                    as.data.frame() %>%
                                    rownames_to_column("gene") %>%
                                    gather("cluster", "value", -gene)))
dat
dat$average_expression[[1]] %>% head()
dat$average_expression[[2]] %>% head()


# Build dataframes for heatmaps -------------------------------------------

pseudocount <- 0.001
limit_high <- 2
limit_low <- -2

heatmap_df_nostim <- dat$average_expression[[1]] %>%
  # Factors to preserve ordering
  mutate(gene = factor(gene, levels = rev(dat$marker_cluster_df[[1]]$gene)),
         cluster = factor(cluster, levels = levels(dat$object[[1]]))) %>%
  # Scale each gene across clusters
  mutate(log2value = log2(value + pseudocount)) %>%
  group_by(gene) %>%
  nest() %>%
  mutate(data = map(data, function(x) mutate(x, zscore = scale(log2value)))) %>%
  unnest(data) %>%
  ungroup() %>%
  # Set data to max/min limits
  mutate(zscore = if_else(zscore > limit_high, limit_high, zscore),
         zscore = if_else(zscore < limit_low, limit_low, zscore))
heatmap_df_nostim

heatmap_df_stim <- dat$average_expression[[2]] %>%
  # Factors to preserve ordering
  mutate(gene = factor(gene, levels = rev(dat$marker_cluster_df[[2]]$gene)),
         cluster = factor(cluster, levels = levels(dat$object[[2]]))) %>%
  # Scale each gene across clusters
  mutate(log2value = log2(value + pseudocount)) %>%
  group_by(gene) %>%
  nest() %>%
  mutate(data = map(data, function(x) mutate(x, zscore = scale(log2value)))) %>%
  unnest(data) %>%
  ungroup() %>%
  # Set data to max/min limits
  mutate(zscore = if_else(zscore > limit_high, limit_high, zscore),
         zscore = if_else(zscore < limit_low, limit_low, zscore))
heatmap_df_stim

dat$heatmap_df <- list(heatmap_df_nostim, heatmap_df_stim)
dat

rm(heatmap_df_nostim, heatmap_df_stim)


# Heatmaps ----------------------------------------------------------------

heatmap_nostim <- dat$heatmap_df[[1]] %>%
  ggplot(aes(x = cluster, y = gene, fill = zscore)) +
  geom_tile() +
  scale_fill_distiller(palette = "PiYG") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom") +
  ggtitle("Resting")

heatmap_stim <- dat$heatmap_df[[2]] %>%
  ggplot(aes(x = cluster, y = gene, fill = zscore)) +
  geom_tile() +
  scale_fill_distiller(palette = "PiYG") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom") +
  ggtitle("Re-stimulated")

dat$heatmap <- list(heatmap_nostim, heatmap_stim)

dat$heatmap[[1]] + dat$heatmap[[2]]

rm(heatmap_nostim, heatmap_stim) # Clean env


# Add horizontal lines between clusters in heatmaps -----------------------

# Need vector of positions for the breaks --> derive from # unique markers
# in each cluster

# First get vector of number of unique gene markers for each cluster
dat <- dat %>%
  mutate(unique_genes_vec = map(marker_cluster_df, 
                                function(x) x %>%
                                  group_by(cluster) %>%
                                  summarise(n_genes = n()) %>%
                                  pull(n_genes) %>%
                                  rev() # Reverse, as ggplot goes from bottom to top on Y axis))
  ))

dat$unique_genes_vec[[1]]

# Convert to vector of positions for line breaks, where each element is added
# to the previous element to get the position

get_line_breaks <- function(x) {
  for (i in seq_along(x)) {
    if(i != 1) {
      x[i] <- x[i] + x[i-1]
    }
  }
  x <- x + 0.5 # Add 0.5 so line break ends up between the sections
  return(x)
}


dat <- dat %>%
  mutate(line_breaks = map(unique_genes_vec, get_line_breaks))
dat
dat$line_breaks[[1]]


# Add line breaks to heatmaps
dat$heatmap[[1]] <- dat$heatmap[[1]] +
  geom_hline(yintercept = dat$line_breaks[[1]])
dat$heatmap[[2]] <- dat$heatmap[[2]] +
  geom_hline(yintercept = dat$line_breaks[[2]])

dat$heatmap[[1]] + dat$heatmap[[2]]



# Get top enriched genes in each cluster ----------------------------------

# For adding beside heatmap

## Add column to id the unique markers previously identified 
dat <- dat %>%
  mutate(unique_markers = map(marker_cluster_df, 
                              function(x) paste0(x$gene,"_",x$cluster)))
dat$unique_markers[[2]] %>% head()

dat <- dat %>%
  mutate(marker_df = map2(
    marker_df, 
    unique_markers,
    function(x,y) x %>%
      mutate(gene_cluster = paste0(gene,"_",cluster),
             unique_marker = if_else(gene_cluster %in% y, TRUE, FALSE))
  ))

dat$marker_df[[2]]

dat <- dat %>%
  mutate(top_markers = map(
    marker_df, 
    function(x) x %>%
      dplyr::filter(unique_marker) %>%
      group_by(cluster) %>%
      nest() %>%
      mutate(top_genes = map(
        data, 
        function(x) x %>%
          dplyr::filter(p_val_adj < 0.25, avg_log2FC > 0) %>%
          arrange(desc(avg_log2FC)) %>%
          slice(1:4) %>%
          pull(gene))) %>%
      select(-data) %>%
      unnest(top_genes)
  ))

dat$top_markers


# Build plots

#   Order by cluster factor levels to match heatmap (by dendogram)
dat <- dat %>%
  mutate(
    top_markers = map2(
      top_markers, object,
      function(x,y) x %>%
        mutate(cluster = factor(cluster, levels = levels(y))) %>%
        arrange(cluster)
    )
  )
dat$top_markers

#   Build coordinate positioning for text labels
dat <- dat %>%
  mutate(
    top_markers = map(
      top_markers,
      function(x) x %>%
        group_by(cluster) %>%
        nest() %>%
        mutate(data = map(
          data, 
          function(x) 
            mutate(x, nudge = seq(-0.3, 0.3, length.out = nrow(x))))) %>%
        unnest(data) %>%
        # Note: as.integer() returns the factor's level positions
        mutate(clust.position = as.integer(cluster) + nudge - 1) 
    ))
dat$top_markers[[2]]

# Plot markers
dat <- dat %>%
  mutate(
    top_marker_plot = 
      map(top_markers,
          function(x) x %>%
            ggplot(aes(x = 1, y = clust.position, label = top_genes)) +
            geom_text(size = 2) +
            geom_hline(yintercept = seq(0.5, 14.5, by = 1), color = "black") +
            scale_y_reverse(breaks = c(15:0)) +
            theme_bw() +
            theme(panel.grid = element_blank(),
                  axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank(),
                  plot.title = element_text(size=10)) +
            xlab("") +
            ggtitle("Markers") +
            coord_cartesian(ylim = c(14.8,0.2))
      ))

dat$top_marker_plot[[1]] + dat$top_marker_plot[[2]]

dat$heatmap[[1]] + dat$top_marker_plot[[1]] + 
  dat$heatmap[[2]] + dat$top_marker_plot[[2]] +
  plot_layout(widths = c(4, 1, 4, 1), nrow = 1)



# Get top guides for each cluster -----------------------------------------

# Get top enriched guides in each cluster
top_guides <- guide_df %>%
  group_by(condition, cluster_id) %>%
  nest() %>%
  mutate(top_genes = map(data, function(x) 
    dplyr::filter(x, q.value < 0.05, log_or > 1) %>%
      arrange(desc(log_or)) %>%
      slice(1:6) %>%
      pull(gene))) %>%
  select(-data) %>%
  unnest(top_genes) %>%
  ungroup() %>%
  rename(cluster = cluster_id)
top_guides

# Add to main dataframe
dat <- dat %>%
  left_join(
    top_guides %>%
      group_by(condition) %>%
      nest() %>%
      rename(top_guides = data),
    by = "condition"
  )

rm(top_guides)

# Re-level and arrange clusters
dat <- dat %>%
  mutate(
    top_guides = map2(
      top_guides, object,
      function(x,y) x %>%
        mutate(cluster = factor(cluster, levels = levels(y))) %>%
        arrange(cluster)
    )
  )
dat$top_guides 


#   Build coordinate positioning for text labels
dat <- dat %>%
  mutate(
    top_guides = map(
      top_guides,
      function(x) x %>%
        group_by(cluster) %>%
        nest() %>%
        mutate(data = map(
          data, 
          function(x) 
            mutate(x, nudge = seq(-0.4, by=0.16, length.out = nrow(x))))) %>%
        unnest(data) %>%
        # Note: as.integer() returns the factor's level positions
        mutate(clust.position = as.integer(cluster) + nudge - 1) 
    ))
dat$top_guides[[2]]

# Plot guides
dat <- dat %>%
  mutate(
    top_guide_plot = 
      map(top_guides,
          function(x) x %>%
            ggplot(aes(x = 1, y = clust.position, label = top_genes)) +
            geom_text(size = 2) +
            geom_hline(yintercept = seq(0.5, 14.5, by = 1), color = "black") +
            scale_y_reverse(breaks = c(15:0)) +
            theme_bw() +
            theme(panel.grid = element_blank(),
                  axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank(),
                  plot.title = element_text(size=10)) +
            xlab("") +
            ggtitle("Guides") +
            coord_cartesian(ylim = c(14.8,0.2))
      ))

dat$top_guide_plot[[1]] + dat$top_guide_plot[[2]]


dat$heatmap[[1]] + dat$top_marker_plot[[1]] + dat$top_guide_plot[[1]] +
  dat$heatmap[[2]] + dat$top_marker_plot[[2]] + dat$top_guide_plot[[2]] +
  plot_layout(widths = c(4, 1, 1, 4, 1, 1), nrow = 1)


# Get top cytokine genes from each cluster --------------------------------

cytokine_genes <- read_tsv("data/GO_0005125_Cytokines.txt", col_names = F)
cytokine_genes <- cytokine_genes$X3 %>% unique() %>% sort()
cytokine_genes

# Get top cyotkine genes (by logFoldChange) that are markers in each cluster
dat <- dat %>%
  mutate(
    top_cytokines_in_cluster = map(
      marker_df,
      function(x) x %>%
        dplyr::filter(gene %in% cytokine_genes & p_val_adj < 0.05 & avg_log2FC > 0) %>%
        group_by(cluster) %>%
        nest() %>%
        mutate(data = map(data, function(x) arrange(x, desc(avg_log2FC)) %>%
                            slice(1:4))) %>%
        unnest(cols = data) %>%
        mutate(top_cytokines_in_cluster = paste0(gene,"_",cluster)) %>%
        pull(top_cytokines_in_cluster)
    )
  )
dat$top_cytokines_in_cluster

dat <- dat %>%
  mutate(
    top_cytokines = map2(
      heatmap_df,
      top_cytokines_in_cluster,
      function(x,y) x %>%
        mutate(gene_cluster = paste0(gene,"_",cluster)) %>%
        dplyr::filter(gene_cluster %in% y) %>%
        select(cluster, gene) %>%
        arrange(cluster)
    )
  )
dat$top_cytokines

#   Build coordinate positioning for text labels
dat <- dat %>%
  mutate(
    top_cytokines = map(
      top_cytokines,
      function(x) x %>%
        group_by(cluster) %>%
        nest() %>%
        mutate(data = map(
          data, 
          function(x) 
            mutate(x, nudge = seq(-0.3, 0.3, length.out = nrow(x))))) %>%
        unnest(data) %>%
        # Note: as.integer() returns the factor's level positions
        mutate(clust.position = as.integer(cluster) + nudge - 1) 
    ))
dat$top_cytokines[[2]]
dat$top_cytokines[[1]]


# Plot cytokines
dat <- dat %>%
  mutate(
    top_cytokine_plot = 
      map(top_cytokines,
          function(x) x %>%
            ggplot(aes(x = 1, y = clust.position, label = gene)) +
            geom_text(size = 2) +
            geom_hline(yintercept = seq(0.5, 14.5, by = 1), color = "black") +
            scale_y_reverse(breaks = c(15:0), position = "right") +
            theme_bw() +
            theme(panel.grid = element_blank(),
                  axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank(),
                  plot.title = element_text(size=10)) +
            xlab("") +
            ggtitle("Cytokines") +
            coord_cartesian(ylim = c(14.8,0.2))
      ))

dat$top_cytokine_plot[[1]] + dat$top_cytokine_plot[[2]] 


dat$heatmap[[1]] + dat$top_marker_plot[[1]] + dat$top_guide_plot[[1]] + dat$top_cytokine_plot[[1]] +
  dat$heatmap[[2]] + dat$top_marker_plot[[2]] + dat$top_guide_plot[[2]] + dat$top_cytokine_plot[[2]] +
  plot_layout(widths = c(4, 1, 1, 1, 4, 1, 1, 1), nrow = 1)


# Label clusters on right -------------------------------------------------


stim_cluster_names <- tibble(
  cluster = levels(dat$object[[2]]),
  cluster_name = c("IFNG High 1",
                   "Negative Regulators",
                   "Th2",
                   "IL2 High",
                   "Proliferative (S)",
                   "Proliferative (G2/M)",
                   "TNF Locus High",
                   "GNLY High",
                   "CCL3/4 High, \nIFNG Low",
                   "CD8 Common",
                   "EMP1 Guides",
                   "IFNG High 2",
                   "IL2 High 1",
                   "CD4 Common",
                   "IL2 High 2")) %>%
  mutate(cluster_name_full = paste0(cluster,": ",cluster_name)) %>%
  pull(cluster_name_full) %>%
  as_factor()
stim_cluster_names

nostim_cluster_names <- tibble(
  cluster = levels(dat$object[[1]]),
  cluster_name = c("Mitochondrial High",
                   "LTB High",
                   "CD8 CCL5 High",
                   "TCR High",
                   "NKG7/GZMA High",
                   "CD4 Common",
                   "EMP1 Guides",
                   "CD8 Common",
                   "CD8 S100B High",
                   "Proliferating (S)",
                   "SERPINB1 High / \nBICDL2/CEACAM1 Guides",
                   "CD4 Effector Memory?",
                   "Proliferating (G2/M)",
                   "MHC Class II High",
                   "GNLY High",
                   "FOXL2NB Guides")) %>%  
  mutate(cluster_name_full = paste0(cluster,": ",cluster_name)) %>%
  pull(cluster_name_full) %>%
  as_factor()
nostim_cluster_names

dat$cluster_names <- list(nostim_cluster_names, stim_cluster_names)
rm(nostim_cluster_names, stim_cluster_names)


dat <- dat %>%
  mutate(
    cluster_id_plot = 
      map(cluster_names,
          function(x) x %>%
            tibble() %>%
            rename("cluster" = ".") %>%
            mutate(clust.position = row_number() - 1) %>%
            ggplot(aes(x = 1, y = clust.position, label = cluster)) +
            geom_text(size = 2) +
            geom_hline(yintercept = seq(0.5, 14.5, by = 1), color = "black") +
            scale_y_reverse(breaks = c(15:0), position = "right") +
            theme_bw() +
            theme(panel.grid = element_blank(),
                  axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank(),
                  plot.title = element_text(size=10)) +
            xlab("") +
            ggtitle("Cluster") +
            coord_cartesian(ylim = c(14.8,0.2))
      ))

dat$cluster_id_plot[[1]] + dat$cluster_id_plot[[2]] 


dat$heatmap[[1]] + dat$top_marker_plot[[1]] + dat$top_guide_plot[[1]] + 
  dat$top_cytokine_plot[[1]] + dat$cluster_id_plot[[1]] +
  dat$heatmap[[2]] + dat$top_marker_plot[[2]] + dat$top_guide_plot[[2]] + 
  dat$top_cytokine_plot[[2]] + dat$cluster_id_plot[[2]] + 
  plot_layout(widths = c(4, 1, 1, 1, 1.5, 4, 1, 1, 1, 1.5), nrow = 1)



# CD4/CD8 scores for each cluster -----------------------------------------

# Re-stimulated 
stim_summary <- dat$object[[2]]@meta.data %>%
  as.data.frame() %>%
  tibble() %>%
  select(seurat_clusters, CD4.CD8.Score) %>%
  group_by(seurat_clusters) %>%
  summarise(CD4.CD8.Score = mean(CD4.CD8.Score))
stim_summary

stim_max <- stim_summary$CD4.CD8.Score %>%
  abs() %>%
  max()

stim_heatmap <- stim_summary %>%
  mutate(seurat_clusters = fct_rev(seurat_clusters)) %>%
  ggplot(aes(x = 1, y = seurat_clusters, fill = CD4.CD8.Score)) +
  geom_tile(color = "black") +
  scale_fill_scico(palette = "cork",
                   limits = c(-1*stim_max, stim_max)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "top") +
  ggtitle("Stim")
stim_heatmap

# Resting 
nostim_summary <- dat$object[[1]]@meta.data %>%
  as.data.frame() %>%
  tibble() %>%
  select(seurat_clusters, CD4.CD8.Score) %>%
  group_by(seurat_clusters) %>%
  summarise(CD4.CD8.Score = mean(CD4.CD8.Score))
nostim_summary

nostim_max <- nostim_summary$CD4.CD8.Score %>%
  abs() %>%
  max()

nostim_heatmap <- nostim_summary %>%
  mutate(seurat_clusters = fct_rev(seurat_clusters)) %>%
  ggplot(aes(x = 1, y = seurat_clusters, fill = CD4.CD8.Score)) +
  geom_tile(color = "black") +
  scale_fill_scico(palette = "cork",
                   limits = c(-1*nostim_max, nostim_max)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "top") +
  ggtitle("NoStim")
nostim_heatmap

nostim_heatmap + stim_heatmap

