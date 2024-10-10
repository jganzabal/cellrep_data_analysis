# Reorder/name clusters based on dendogram (tree)

# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)


# Import Data -------------------------------------------------------------

dat <- tibble(
  condition = c("NoStim", "Stim"),
  object = list(
    LoadH5Seurat("data_out/Nostim.h5Seurat"),
    LoadH5Seurat("data_out/Stim.h5Seurat")
  )
)
dat


# Get cluster order for plotting ------------------------------------------

# Use the order they appear in the cluster dendogram
#   First build dendograms
dat <- dat %>%
  mutate(object = map(object, BuildClusterTree))
# Note: for some reason the arg reorder = TRUE is not working for BuildClusterTree
# so using another solution

# Now get tree order

# solution from:
# https://stackoverflow.com/questions/34364660/how-to-get-correct-order-of-tip-labels-in-ape-after-calling-ladderize-function

get_tree_order <- function(obj) {
  tree <- obj@tools$.f
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  ordered_tips <- ordered_tips - 1 # Something weird with the clusters starting at 0
  return(ordered_tips)
}

dat <- dat %>%
  mutate(tree_order = map(object, get_tree_order))
dat

# NoStim
plot(dat$object[[1]]@tools$.f)
dat$tree_order[[1]]

# Stim
plot(dat$object[[2]]@tools$.f)
dat$tree_order[[2]]


# Relevel cluster metadata

dat$object[[1]]@meta.data$seurat_clusters <- factor(dat$object[[1]]@meta.data$seurat_clusters,
                                                    levels = dat$tree_order[[1]])
dat$object[[1]]@meta.data$seurat_clusters %>% head()

dat$object[[2]]@meta.data$seurat_clusters <- factor(dat$object[[2]]@meta.data$seurat_clusters,
                                                    levels = dat$tree_order[[2]])
dat$object[[2]]@meta.data$seurat_clusters %>% head()

# Cluster conversion tables
nostim_cluster_conversion <- tibble(
  original_id = levels(dat$object[[1]]@meta.data$seurat_clusters),
  new_id = c(1:length(levels(dat$object[[1]]@meta.data$seurat_clusters))) %>%
    as.character())
nostim_cluster_conversion

stim_cluster_conversion <- tibble(
  original_id = levels(dat$object[[2]]@meta.data$seurat_clusters),
  new_id = c(1:length(levels(dat$object[[2]]@meta.data$seurat_clusters))) %>%
    as.character())
stim_cluster_conversion

cluster_conversion <- bind_rows(
  mutate(nostim_cluster_conversion, condition = "NoStim"),
  mutate(stim_cluster_conversion, condition = "Stim")
)
cluster_conversion

# Convert to new cluster ids based on tree order
nostim.new.cluster.ids <- nostim_cluster_conversion$new_id
names(nostim.new.cluster.ids) <- nostim_cluster_conversion$original_id
dat$object[[1]] <- RenameIdents(dat$object[[1]], nostim.new.cluster.ids)
DimPlot(dat$object[[1]], reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
levels(dat$object[[1]]@meta.data$seurat_clusters) <- nostim_cluster_conversion$new_id

stim.new.cluster.ids <- stim_cluster_conversion$new_id
names(stim.new.cluster.ids) <- stim_cluster_conversion$original_id
dat$object[[2]] <- RenameIdents(dat$object[[2]], stim.new.cluster.ids)
DimPlot(dat$object[[2]], reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
levels(dat$object[[2]]@meta.data$seurat_clusters) <- stim_cluster_conversion$new_id

